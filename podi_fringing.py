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

import sys
import os
import pyfits
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import scipy
import scipy.stats


import Queue
import threading
import multiprocessing
import ctypes

from podi_definitions import *
import podi_imcombine

avg_sky_countrates = {
    "odi_i": 3.8,
    "odi_z": 4.0,
    "823_v2": 0.1,
}


def make_fringing_template(input_filelist, outputfile, return_hdu=False, skymode='local'):

    # First loop over all filenames and make sure all files exist
    hdu_filelist = []
    for file in input_filelist:
        if (os.path.isfile(file)):
            hdu_filelist.append(pyfits.open(file))

    if (len(hdu_filelist) <= 0):
        stdout_write("No existing files found in input list, hence nothing to do!\n")
        return
    
    # Read the input parameters
    # Note that file headers are copied from the first file

    # Create the primary extension of the output file
    ref_hdulist = hdu_filelist[0]
    primhdu = pyfits.PrimaryHDU(header=ref_hdulist[0].header)

    # Add PrimaryHDU to list of OTAs that go into the output file
    out_hdulist = [primhdu]

    filtername = primhdu.header['FILTER']
    if (outputfile == "auto"):
        outputfile = "fringes__%s.fits" % (filtername)

    print "Output file=",outputfile


    #
    # Now loop over all extensions and compute the mean
    #
    for cur_ext in range(1, len(ref_hdulist)):
        
        data_blocks = []
        # Check what OTA we are dealing with
        if (not is_image_extension(ref_hdulist[cur_ext].header)):
            continue
        extname = ref_hdulist[cur_ext].header['EXTNAME']

        if (filtername in otas_for_photometry):
            useful_otas = otas_for_photometry[filtername]
            ota_id = int(extname[3:5])
            if (not ota_id in useful_otas):
                continue
        
        stdout_write("\rCombining frames for OTA %s (#% 2d/% 2d) ..." % (extname, cur_ext, len(ref_hdulist)-1))

        # Now open all the other files, look for the right extension, and copy their image data to buffer
        for file_number in range(0, len(filelist)):

            try:
                this_hdu = hdu_filelist[file_number][extname]
            except:
                continue

            # Skip all OTAs that are marked as video/guide OTAs
            cellmode = this_hdu.header['CELLMODE']
            if (cellmode.find("V") >= 0):
                continue
            
            skylevel = this_hdu.header['SKY_MEDI']
            if (skymode == 'global'):
                skylevel = hdu_filelist[file_number][0].header['SKYLEVEL']
            if ("EXPTIME" in hdu_filelist[file_number][0].header):
                exptime = hdu_filelist[file_number][0].header['EXPTIME']
                filter = hdu_filelist[file_number][0].header['FILTER']
                if (filter in avg_sky_countrates):
                    max_skylevel = 2 * avg_sky_countrates[filter] * exptime
                    if (skylevel > max_skylevel):
                        stdout_write(" (%.1f)" % (skylevel/exptime))
                        continue

            fringing = (this_hdu.data - skylevel) / skylevel
            stdout_write(" %.1f" % (skylevel/exptime))
            data_blocks.append(fringing)

            # delete the data block to free some memory, since we won't need it anymore
            del this_hdu.data

        stdout_write(" combining ...")
        #combined = podi_imcombine.imcombine_data(data_blocks, "nanmedian")
        combined = podi_imcombine.imcombine_data(data_blocks, "nanmedian.bn")

        # Create new ImageHDU
        # Insert the imcombine'd frame into the output HDU
        # Copy all headers from the reference HDU
        hdu = pyfits.ImageHDU(header=ref_hdulist[cur_ext].header, data=combined)

        # Append the new HDU to the list of result HDUs
        out_hdulist.append(hdu)
        stdout_write(" done!\n")

        del hdu

    # Add the assumed skylevel countrate to primary header so we can use it
    # when it comes to correcting actual data
    out_hdulist[0].header.update("SKYCNTRT", avg_sky_countrates[filter], "average countrate of dark sky")

    return_hdu = False
    out_hdu = pyfits.HDUList(out_hdulist)
    if (not return_hdu and outputfile != None):
        stdout_write(" writing results to file %s ..." % (outputfile))
        clobberfile(outputfile)
        out_hdu.writeto(outputfile, clobber=True)
        out_hdu.close()
        del out_hdu
        del out_hdulist
        stdout_write(" done!\n")
    elif (return_hdu):
        stdout_write(" returning HDU for further processing ...\n")
        return out_hdu
    else:
        stdout_write(" couldn't write output file, no filename given!\n")

    return

if __name__ == "__main__":

    if (cmdline_arg_isset("-singles")):
        for filename in get_clean_cmdline()[1:]:
            outputfile = filename[:-5]+".fringe.fits"
            if (cmdline_arg_isset("-noclobber") and os.path.isfile(outputfile)):
                stdout_write("\n%s already exists, skipping!\n\n" % (outputfile))
                continue

            stdout_write("Converting %s to fringemask ...\n" % (filename))
            hdulist = pyfits.open(filename)

            out_hdu = [pyfits.PrimaryHDU(header=hdulist[0].header)]

            for ext in range(len(hdulist)):
                if (not is_image_extension(hdulist[ext].header)):
                    continue
                
                # Skip all OTAs that are marked as video/guide OTAs
                cellmode = hdulist[ext].header['CELLMODE']
                if (cellmode.find("V") >= 0):
                    continue

                skylevel = hdulist[ext].header['SKY_MEDI']
                
                fringing = (hdulist[ext].data - skylevel) / skylevel
                stdout_write("   %s = %.1f\n" % (hdulist[ext].header['EXTNAME'], skylevel))

                out_hdu.append(pyfits.ImageHDU(header=hdulist[ext].header,
                                               data=fringing))

            stdout_write("writing (%s)..." % (outputfile))
            out_hdulist = pyfits.HDUList(out_hdu)
            out_hdulist.writeto(outputfile, clobber=True)
            stdout_write(" done!\n\n")

        sys.exit(0)
                               
    
    if (cmdline_arg_isset("-make_template")):
        outputfile = get_clean_cmdline()[1]
        filelist = get_clean_cmdline()[2:]
        operation = cmdline_arg_set_or_default("-op", "mean")
        operation = cmdline_arg_set_or_default("-op", "nanmedian.bn")
        make_fringing_template(filelist, outputfile, operation)

        sys.exit(0)

    if (cmdline_arg_isset("-samplesky")):
        import podi_fitskybackground
        boxsize = int(cmdline_arg_set_or_default("-boxsize", "15"))
        count = int(cmdline_arg_set_or_default("-count", "12500"))

        print "Using %d boxes with +/- %d pixel width" % (count, boxsize)

        filename = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]

        print filename,"-->",outputfile
        hdulist = pyfits.open(filename)
        for i in range(len(hdulist)):
            if (is_image_extension(hdulist[i].header)):
                extname = hdulist[i].header['EXTNAME']
                stdout_write("\rWorking on %s ..." % (extname))
                regions = numpy.array(podi_fitskybackground.sample_background(hdulist[i].data, None, None, min_found=count, boxwidth=boxsize))
                
                median = numpy.median(regions[:,4])
                std = numpy.std(regions[:,4])

                histogram, binedges = numpy.histogram(regions[:,4], bins=500, range=(median-20*std,median+20*std))
                nice_histogram = numpy.empty(shape=(histogram.shape[0], 3))
                nice_histogram[:,0] = binedges[:-1]
                nice_histogram[:,1] = binedges[1:]
                nice_histogram[:,2] = histogram[:]
                
                dumpfile = "%s.%s" % (outputfile, extname)
                df = open(dumpfile, "w")
                numpy.savetxt(df, regions)
                print >>df, "\n\n\n\n\n\n\n"
                numpy.savetxt(df, nice_histogram)
                print >>df, "\n\n\n\n\n\n\n"

                validpixels = hdulist[i].data
                fullframe, edges = numpy.histogram(validpixels, bins=500, range=(median-20*std,median+20*std))
                nice_histogram[:,2] = fullframe[:]
                numpy.savetxt(df, nice_histogram)

                df.close()
            

    if (cmdline_arg_isset("-correct")):

        inputframe = get_clean_cmdline()[1]
        template = get_clean_cmdline()[2]
        scaling = int(get_clean_cmdline()[3])
        output = get_clean_cmdline()[4]

        inputhdu = pyfits.open(inputframe)
        templatehdu = pyfits.open(template)

        for i in range(len(inputhdu)):
            if (not is_image_extension(inputhdu[i].header)):
                continue

            extname = inputhdu[i].header['EXTNAME']
            stdout_write("\rWorking on %s ... " % (extname))

            try:
                fringemap = templatehdu[extname].data * scaling
                inputhdu[i].data -= fringemap
            except:
                print "Problem with extension",extname
                pass
                continue

        stdout_write(" writing ...")

        inputhdu.writeto(output, clobber=True)
        
        stdout_write(" done!\n")

            
        sys.exit(0)
