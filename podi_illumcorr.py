#! /usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

"""

This module contains all functionality related to fringin in ODI frames, from
creating the fringe templates, to finding the ideal fringe scaling factor to
actually subtracting the fringe frame.

Standalone routines
----------------------

* **-make_template**

  Create a fringe template from collectcells-reduced frames

  ``./podi_fringing.py -make_template (-op=nanmedian.bn) output_template.fits file1.fits file2.fits``

* **-esomethod**

  Determine the optimal fringe scaling and perform fringe removal

  ``./podi_fringing.py -esomethod fringe_template.fits input.fits output.fits``


Methods
---------------------------

"""



import sys
import os
import pyfits
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import scipy
import scipy.stats
import scipy.optimize
import bottleneck
import scipy.signal

import Queue
import threading
import multiprocessing
import ctypes
import time

from podi_definitions import *
from podi_commandline import *
import podi_imcombine
import podi_fitskybackground
import time



if __name__ == "__main__":

    filelist = get_clean_cmdline()[2:]
    outfile = get_clean_cmdline()[1]

    masked_list = []

    for fitsfile in filelist:

        # get some info so we know what to call the output frame
        hdulist = pyfits.open(fitsfile)
        obsid = hdulist[0].header['OBSID']
        
        # Run Sextractor
        segmask = "%s_segmentation.fits" % (obsid)

        if (not os.path.isfile(segmask)):
            sex_cmd = """%(sex)s -c %(conf)s \
                         -PARAMETERS_NAME %(params)s \
                         -CHECKIMAGE_TYPE SEGMENTATION \
                         -CHECKIMAGE_NAME %(segfile)s \
                         -FILTER_NAME %(filtername)s \
                         %(image)s
                """ % {
                'sex': sitesetup.sextractor,
                'conf': "%s/.config/illumcorr.conf" % (sitesetup.exec_dir),
                'params': "%s/.config/illumcorr.param" % (sitesetup.exec_dir),
                'filtername': "%s/.config/gauss_5.0_9x9.conv" % (sitesetup.exec_dir),
                'segfile': segmask,
                'image': fitsfile,
                }
            print " ".join(sex_cmd.split())

            start_time = time.time()
            try:
                ret = subprocess.Popen(sex_cmd.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
                (sex_stdout, sex_stderr) = ret.communicate()
                if (not ret.returncode == 0):
                    print sex_stdout
                    print sex_stderr
            except OSError as e:
                print >>sys.stderr, "Execution failed:", e
            end_time = time.time()
            print("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
        else:
            print "segmentation mask (%s) exist, re-using it" % (segmask)

        #
        # Now use the mask and the input frame to mask out 
        # all sources in the frame. At the same time, rescale the frame by
        # dividing the intensity by the global median skylevel
        #
        hdu_out = []
        hdu_out.append(hdulist[0])

        masked_frame = "%s_masked.fits" % (obsid)
        if (not os.path.isfile(masked_frame)):
            mask_hdu = pyfits.open(segmask)

            for ext in hdulist:
                if (not is_image_extension(ext)):
                    continue

                if (ext.header['CELLMODE'].find("V") >= 0):
                    # This is a video extension, do not use it
                    continue

                # Now search for the right extension in the mask frame
                found_mask = False
                for mask_ext in mask_hdu:
                    if (ext.name == mask_ext.name):
                        # found it
                        found_mask = True
                        print "found the mask for extension",ext.name

                        # Set all detected pixels to NaN to ignore them during the
                        # final imcombine
                        ext.data[mask_ext.data > 0] = numpy.NaN

                        # Rescale with the global sky-level
                        # maybe better to re-compute based on the segmentation mask
                        ext.data /= hdulist[0].header['SKYLEVEL']

                        hdu_out.append(ext)

                if (not found_mask):
                    print "Can't find extension ",ext.name," in mask"


            print "writing masked frame to", masked_frame

            clobberfile(masked_frame)
            hdulist_out = pyfits.HDUList(hdu_out)
            hdulist_out.writeto(masked_frame, clobber=True)

            mask_hdu.close()

        masked_list.append(masked_frame)

        print "done with this one, taking next frame"

    #
    # Stack all files with all stars masked out
    #
    print "Stacking the following frames:"
    print "\n".join(["  ** %s" % a for a in masked_list])
    presmoothed_file = outfile[:-5]+".raw.fits"
    if (not os.path.isfile(presmoothed_file)):
        combined = podi_imcombine.imcombine(input_filelist=masked_list,
                                            outputfile=None,
                                            operation="sigmaclipmedian",
                                            return_hdu=True,
                                            subtract=None,
                                            scale=None)

        combined.writeto(outfile[:-5]+".raw.fits")
    else:
        print "reading pre-smoothed file from file"
        combined = pyfits.open(presmoothed_file)

    # Now go through all extensions and apply some smoothing or 
    # low-pass filtering to further increase signal-to-noise
    # combined = pyfits.open("illumcorr.fits")

    # smooth_filter_1d = scipy.signal.flattop(10, True).reshape(1,10)
    # smooth_filter_2d = numpy.sqrt(smooth_filter_1d * smooth_filter_1d.T)
    # print numpy.sum(smooth_filter_1d)
    # print numpy.sum(smooth_filter_2d)

    smooth_filter_2d = numpy.ones((5,5))
    smooth_filter_2d /= numpy.sum(smooth_filter_2d)

    for ext in combined:
        if (not is_image_extension(ext)):
            continue

        print "smoothing ", ext.name

        data = ext.data.copy()
        data[numpy.isnan(data)] = 1.0
        data_out = scipy.signal.convolve2d(data, 
                                           smooth_filter_2d,
                                           mode='same',
                                           boundary='fill',
                                           fillvalue=1)

        data_out[numpy.isnan(ext.data)] = numpy.NaN
        ext.data = data_out

    combined.writeto(outfile, clobber=True)
