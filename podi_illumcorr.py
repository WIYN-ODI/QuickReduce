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

This module contains all functionality related to illumincation corrections 
for ODI frames, from creating the templates, to applying the correction to 
individual files or a bunch of files.

Standalone routines
----------------------

* **-create**

  Create an illumination correction from QR-reduced frames

  ``./podi_illumcorr.py -create output.fits file1.fits file2.fits``

* **-apply1**

  Apply the fringe correction to a single frame

  ``./podi_illumcorr.py -apply1 illumcorr_template.fits input.fits output.fits``


* **-applyall**

  Apply the finge correction to a bunch of files in one execution

  ``./podi_illumcorr.py -applyall illumcorr_template.fits INSERT file1.fits file2.fits``

  In this mode, the output file is created from the input file by inserting 
  INSERT between the filename and the .fits extension. 

  Example: file1.fits will be corrected ans saved as file1.INSERT.fits


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
import podi_imarith
import podi_logging


def get_illumination_filename(calibdir, filtername, binning):

    if (os.path.isfile(calibdir)):
        return calibdir

    elif (os.path.isdir(calibdir)):
        fn = "%s/illumination__%s_bin%d.fits" % (calibdir, filtername, binning)
        return fn

    elif (calibdir == None):
        fn = "illumination__%s_bin%d.fits" % (filtername, binning)
        return fn

    return calibdir
    

def apply_illumination_correction(input_hdu, illumcorr_hdu, output_file=None):

    logger = logging.getLogger("ApplyIllumCorr")

    input_is_hdu = True
    if (not type(input_hdu) == pyfits.hdu.hdulist.HDUList):
        input_is_hdu = False
        try:
            if (os.path.isfile(input_hdu)):
                hdu = pyfits.open(input_hdu)
                input_hdu = hdu
            else:
                logger.error("The specified illumination correction frame (%s) does not exist" % (
                    str(input_hdu)))
                return None
        except:
            logger.error("Can't figure out what format input_hdu is in")
            return None

    illum_is_hdu = True
    if (not type(illumcorr_hdu) == pyfits.hdu.hdulist.HDUList):
        illum_is_hdu = False
        try:
            hdu = pyfits.open(illumcorr_hdu)
            illumcorr_hdu = hdu
        except:
            logger.error("Can't figure out what format illumcorr_hdu is in")
            return None

    for ext in input_hdu:
        if (not is_image_extension(ext)):
            continue

        # Now search for the right illumination frame
        for il in illumcorr_hdu:
            if (ext.name == il.name):
                ext.data /= il.data
                break
    
    if (not illum_is_hdu):
        illumcorr_hdu.close()

    if (not input_is_hdu and not output_file == None):
        clobberfile(output_file)
        input_hdu.writeto(output_file, clobber=True)
        return output_file

    
    return input_hdu

def compute_illumination_frame(queue, return_queue, tmp_dir=".", redo=False):

    root = logging.getLogger("CompIllumination")

    while (True):
        cmd = queue.get()
        if (cmd == None):
            root.debug("Received shutdown command")
            queue.task_done()
            return

        fitsfile = cmd
        root.debug("Received new work: %s" % (fitsfile))

        tempfiles = []

        # get some info so we know what to call the output frame
        hdulist = pyfits.open(fitsfile)
        obsid = hdulist[0].header['OBSID']
        logger = logging.getLogger("CompIllum(%s)" % (obsid))

        # Run Sextractor
        segmask = "%s/%s_segmentation.fits" % (tmp_dir, obsid)
        masked_frame = "%s/%s_masked.fits" % (tmp_dir, obsid)

        if (not (os.path.isfile(segmask) and os.path.isfile(masked_frame)) or redo):
            logger.info("Starting work (source detection, masking, normalization) ...")
        else:
            logger.info("No work necessary, re-using existing data ...")
            
        if (not os.path.isfile(segmask) or redo):
            logger.debug("Creating segmentation mask: %s" % (segmask))
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
            
            logger.debug("Starting Sextractor:\n%s" % (" ".join(sex_cmd.split())))

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
            logger.debug("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
        else:
            logger.debug("segmentation mask (%s) exist, re-using it" % (segmask))

        #
        # Now use the mask and the input frame to mask out 
        # all sources in the frame. At the same time, rescale the frame by
        # dividing the intensity by the global median skylevel
        #
        hdu_out = []
        hdu_out.append(hdulist[0])

        if (not os.path.isfile(masked_frame) or redo):
            logger.debug("Preparing masked frame: %s" % (masked_frame))

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
                        logger.debug("found the mask for extension %s" % (ext.name))
                        
                        mask_grown = scipy.ndimage.filters.convolve(
                            input=mask_ext.data, 
                            weights=numpy.ones((10,10)), 
                            output=None, 
                            mode='constant', cval=0.0)

                        # Set all detected pixels to NaN to ignore them during the
                        # final imcombine
#                        ext.data[mask_ext.data > 0] = numpy.NaN
                        ext.data[mask_grown > 0] = numpy.NaN

                        # Rescale with the global sky-level
                        # maybe better to re-compute based on the segmentation mask
                        ext.data /= hdulist[0].header['SKYLEVEL']

                        hdu_out.append(ext)

                if (not found_mask):
                    logger.debug("Can't find extension %s in mask" % (ext.name))


            logger.debug("writing masked frame to %s" % (masked_frame))

            clobberfile(masked_frame)
            hdulist_out = pyfits.HDUList(hdu_out)
            hdulist_out.writeto(masked_frame, clobber=True)

            mask_hdu.close()

        return_queue.put(masked_frame)
        queue.task_done()

        logger.debug("done with this one, taking next frame")

    root.debug("Terminating process")
    return


def prepare_illumination_correction(filelist, outfile, tmpdir=".", redo=False):

    logger = logging.getLogger("CreateIllumCorr")

    number_files_sent_off = 0
    queue = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.Queue()

    logger.info("Preparing all individual illumination reference frames")
    for fitsfile in filelist:
        logger.debug("Queuing %s" % (fitsfile))
        queue.put((fitsfile))
        number_files_sent_off += 1

    processes = []
    for i in range(sitesetup.number_cpus):
        logger.debug("Starting process #%d" % (i+1))
        p = multiprocessing.Process(target=compute_illumination_frame,
                                    kwargs = {'queue': queue,
                                              'return_queue': return_queue,
                                              'tmp_dir': tmpdir,
                                              'redo': redo,
                                          },
                                    # args=(queue, return_queue),
        )
        queue.put(None)
        p.start()
        processes.append(p)

    masked_list = []
    for i in range(number_files_sent_off):
        masked_frame = return_queue.get()
        if (not masked_frame == None):
            masked_list.append(masked_frame)

    for p in processes:
        p.join()

    logger.info("All files prepared, combining ...")

    # Pick the first file, and extract some information about filter and binning
    filtername = "unknown"
    binning = 0
    for infile in filelist:
        try:
            hdulist = pyfits.open(infile)
            filtername = hdulist[0].header['FILTER']
            binning = hdulist[0].header['BINNING']
            hdulist.close()
            break
        except:
            pass

    illumcorrfile = get_illumination_filename(outfile, filtername, binning)
    #print "\n"*5,illumcorrfile,"\n",outfile,"\n"*5

    presmoothed_file = illumcorrfile[:-5]+".raw.fits"

    #
    # Stack all files with all stars masked out
    #
    logger.debug("Stacking the following frames:\n%s" % (
        "\n".join(["  ** %s" % a for a in masked_list])))

    if (not os.path.isfile(presmoothed_file) or redo):
        logger.debug("Starting imcombine!")
        combined = podi_imcombine.imcombine(input_filelist=masked_list,
                                            outputfile=None,
                                            operation="sigmaclipmedian",
                                            return_hdu=True,
                                            subtract=None,
                                            scale=None)

        clobberfile(presmoothed_file)
        logger.debug("Writing stacked, pre-smoothed frame to %s" % (presmoothed_file))
        combined.writeto(presmoothed_file)
    else:
        logger.debug("Found raw (pre-smoothed) combined frame, reusing it")
        combined = pyfits.open(presmoothed_file)

    # Now go through all extensions and apply some smoothing or 
    # low-pass filtering to further increase signal-to-noise
    # combined = pyfits.open("illumcorr.fits")

    # smooth_filter_1d = scipy.signal.flattop(10, True).reshape(1,10)
    # smooth_filter_2d = numpy.sqrt(smooth_filter_1d * smooth_filter_1d.T)
    # print numpy.sum(smooth_filter_1d)
    # print numpy.sum(smooth_filter_2d)

    logger.info("Starting smoothing to increase signal-to-noise")

    smooth_filter_2d = numpy.ones((5,5))
    smooth_filter_2d /= numpy.sum(smooth_filter_2d)

    for ext in combined:
        if (not is_image_extension(ext)):
            continue

        logger.debug("smoothing %s" % (ext.name))

        data = ext.data.copy()
        data[numpy.isnan(data)] = 1.0
        # data_out = scipy.signal.convolve2d(data, 
        #                                    smooth_filter_2d,
        #                                    mode='same',
        #                                    boundary='fill',
        #                                    fillvalue=1)
        data_out = scipy.ndimage.filters.median_filter(input=data, 
                                                       size=3,
                                                       mode='constant', cval=1.0)

        data_out[numpy.isnan(ext.data)] = numpy.NaN
        ext.data = data_out

    combined.writeto(illumcorrfile, clobber=True)
    logger.info("Work done, output written to %s" % (outfile))

    return

if __name__ == "__main__":

    # Setup everything we need for logging
    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)

    if (cmdline_arg_isset("-create")):

        filelist = get_clean_cmdline()[2:]
        outfile = get_clean_cmdline()[1]
        tmpdir = cmdline_arg_set_or_default("-tmp", sitesetup.swarp_singledir)
        redo = cmdline_arg_isset("-redo")

        print outfile, "\n"*5

        prepare_illumination_correction(filelist, outfile, tmpdir, redo)

    elif (cmdline_arg_isset("-apply1")):

        illumcorr = get_clean_cmdline()[1]
        inputfile = get_clean_cmdline()[2]
        outputfile = get_clean_cmdline()[3]

        podi_imarith.imarith(inputfile, "/", illumcorr, outputfile)

    elif (cmdline_arg_isset("-applyall")):

        illumcorr = get_clean_cmdline()[1]
        insert = get_clean_cmdline()[2]

        for inputfile in get_clean_cmdline()[3:]:
            outputfile = "%s.%s.fits" % (inputfile[:-5], insert)
            podi_imarith.imarith(inputfile, "/", illumcorr, outputfile)



        pass

    else:

        print """\
Options are:

  -create output.fits input*.fits

  -apply1 illumcorr.fits input.fits output.fits

  -applyall illumcorr.fits insert_string input*.fits

"""
        
    # Shutdown logging to shutdown cleanly
    podi_logging.shutdown_logging(options)
