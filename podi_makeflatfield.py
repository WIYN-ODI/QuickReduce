#! /usr/bin/env python3
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

This module handles the normalization of flat-fields.

"""

import sys
import os
import astropy.io.fits as pyfits
import numpy
import scipy

from podi_definitions import *
from podi_commandline import *
import podi_focalplanelayout

import podi_logging
import logging

def normalize_flatfield(filename, outputfile, 
                        binning_x=8, binning_y=8, 
                        repeats=3, 
                        batchmode_hdu=None,
                        normalize_otas=None):

    logger = logging.getLogger("NormFlatField")
    logger.debug("Starting to normalize %s" % (str(batchmode_hdu)))

    if (batchmode_hdu is not None):
        hdulist = batchmode_hdu
    else:
        hdulist = pyfits.open(filename)

    filter = hdulist[0].header['FILTER']

    fpl = podi_focalplanelayout.FocalPlaneLayout(hdulist)

    list_of_otas_to_normalize = fpl.get_science_area_otas(filter, include_vignetted=False)
    if (normalize_otas is not None):
        list_of_otas_to_normalize = normalize_otas

    logger.info("Using these OTAs to normalize overall flux:\n%s" % (", ".join(["%02d" % ota for ota in list_of_otas_to_normalize])))

    flatfield_data = numpy.zeros(
        shape=(len(list_of_otas_to_normalize)*4096*4096//(binning_x*binning_y)),
        dtype=numpy.float32)
    flatfield_data[:] = numpy.nan

    # also prepare to store the global gain value
    gain_sum = 0
    gain_count = 0

    datapos = 0
    for extension in range(1, len(hdulist)): #hdulist[1:]:
        if (not is_image_extension(hdulist[extension])):
            continue

        fppos = int(hdulist[extension].header['FPPOS'][2:4])

        #print list_of_otas_to_normalize
        try:
            index = list_of_otas_to_normalize.index(fppos)
        except:
            # We didn't find this OTA in the list, so skip it
            hdulist[extension].header["FF_NORM"] = (False, "Used in normalization")
            extension += 1
            continue

        hdulist[extension].header["FF_NORM"] = (True, "Used in normalization")

        gain_ota = hdulist[extension].header['GAIN']
        gain_ota_count = hdulist[extension].header['NGAIN']
        
        gain_sum += gain_ota * gain_ota_count
        gain_count += gain_ota_count

        # We now know that we should include this OTA in the
        # calculation of the flat-field normalization
        logger.debug("Adding OTA %02d to flat-field ..." % fppos)
        #flatfield_data = numpy.concatenate((flatfield_data, extension.data.flatten()))
        #flatfield_data[extension,:,:] = extension.data

        if (binning_x>1 or binning_y>1):
            sx, sy = hdulist[extension].data.shape[0], hdulist[extension].data.shape[1]
            bx, by = sx//binning_x, sy//binning_y
            one_d = numpy.reshape(hdulist[extension].data, (by,binning_y,bx,binning_x)).mean(axis=-1).mean(axis=1).flatten()
        else:
            one_d = hdulist[extension].data.flatten()
            
        flatfield_data[datapos:datapos+one_d.shape[0]] = one_d

        datapos += one_d.shape[0]
        #print datapos
        del one_d

    # Remove all remaining NaN values and atruncate the array to the values actually used
    finite = numpy.isfinite(flatfield_data[:datapos])
    flatfield_data = flatfield_data[:datapos][finite]

    # Now we are through all flatfields, compute the median value
    logger.debug(" computing median ...")
    sigma_min, sigma_max = -1e5, 1e6
    for i in range(repeats):

        valid = (flatfield_data > sigma_min) & (flatfield_data < sigma_max)

        ff_median_level = numpy.median(flatfield_data[valid])
        ff_std = numpy.std(flatfield_data[valid])

        sigma_min = ff_median_level - 2 * ff_std
        sigma_max = ff_median_level + 3 * ff_std
        #print i, numpy.sum(valid), datapos, ff_median_level, ff_std, sigma_min, sigma_max
        
    if (ff_median_level <= 0):
        logger.error("Something went wrong or this is no flatfield frame")
        ff_median_level = 1.0

    #stdout_write("\b\b\b(% 7.1f) ..." % (ff_median_level))
    logger.debug("Found median level % 7.1f ADU, normalizing ..." % (ff_median_level))

    # Now normalize all OTAs with the median flatfield level
    #stdout_write(" normalizing ...")

    # Create a new HDU list for the normalized output
    # hdu_out = [] #pyfits.PrimaryHDU(header=hdulist[0].header)]

    # for extension in range(0, len(hdulist)):
    #     if (not is_image_extension(hdulist[extension])):
    #         hdu_out.append(hdulist[extension])
    #         continue

    #     data = hdulist[extension].data.copy()
    #     data /= ff_median_level
    #     data[data < 0.1] = numpy.nan
    #     new_hdu = pyfits.ImageHDU(data=data, header=hdulist[extension].header)

    #     #hdulist[extension].data /= ff_median_level
    #     #hdulist[extension].data[hdulist[extension].data < 0.1] = numpy.nan
    #     new_hdu.header.add_history("FF-level: %.1f" % (ff_median_level))
    #     hdu_out.append(new_hdu)
        
    # hdulist = pyfits.HDUList(hdu_out)

    hdu_out = [] #pyfits.PrimaryHDU(header=hdulist[0].header)]

    for extension in hdulist:
        if (not is_image_extension(extension)):
            continue

        # data = hdulist[extension].data.copy()
        # data /= ff_median_level
        # data[data < 0.1] = numpy.nan
        # new_hdu = pyfits.ImageHDU(data=data, header=hdulist[extension].header)
        extension.data /= ff_median_level
        extension.data[extension.data < 0.1] = numpy.nan

        #hdulist[extension].data /= ff_median_level
        #hdulist[extension].data[hdulist[extension].data < 0.1] = numpy.nan
        # new_hdu.header.add_history("FF-level: %.1f" % (ff_median_level))
        # hdu_out.append(new_hdu)
        
    # hdulist = pyfits.HDUList(hdu_out)

    #
    # compute the global gain value and store it in primary header
    #
    logger.debug("Computing global gain value (sum=%.1f, #=%d)" % (gain_sum, gain_count))
    global_gain = gain_sum / gain_count if (gain_count > 0) else -1
    hdulist[0].header['GAIN'] = global_gain if (gain_count > 0) else -1.
    hdulist[0].header['NGAIN'] = gain_count

    logger.debug("writing results to file (%s) ..." % (outputfile))
    clobberfile(outputfile)
    hdulist.writeto(outputfile, overwrite=True)
    logger.info("done!")
       
if __name__ == "__main__":

    binning_x = int(cmdline_arg_set_or_default("-binx", 8))
    binning_y = int(cmdline_arg_set_or_default("-biny", 8))
    repeats = int(cmdline_arg_set_or_default("-reps", 3))
    
    clean_list = get_clean_cmdline()
    if (cmdline_arg_isset("-multi")):
        #
        # Ok, we should work on a number of files
        #
        add_filter_to_output_filename = cmdline_arg_isset("-addfilter")
        for filename in clean_list[1:]:
            directory, basename = os.path.split(filename)

            # If -keeppath is specified, put the normalized output file in
            # the same directory where we got the files from
            if (cmdline_arg_isset("-keeppath")):
                if (directory == ""):
                    out_directory = "."
                else:
                    out_directory = directory
            #
            # Or maybe the user specified the output directory explicitely?
            #
            elif (cmdline_arg_isset("-outdir")):
                out_directory = get_cmdline_arg("-outdir")
            #
            # by default, put all output files in the current directory
            #
            else:
                out_directory = "."

            #
            # Now construct the output filename
            #
            if (add_filter_to_output_filename):
                hdulist = pyfits.open(filename)
                filter = hdulist[1].header['FILTER']
                outputfile = "%s/%s.norm.%s.fits" % (out_directory, basename[:-5], filter)
                hdulist.close()
                del hdulist
            else:
                outputfile = "%s/%s.norm.fits" % (out_directory, basename[:-5])
            print(filename, outputfile)

            # And finally, do the actual work
            normalize_flatfield(filename, outputfile, binning_x=binning_x, binning_y=binning_y, repeats=repeats)
    else:
        filename = clean_list[1]
        outputfile = clean_list[2]
        normalize_flatfield(filename, outputfile, binning_x=binning_x, binning_y=binning_y, repeats=repeats)

