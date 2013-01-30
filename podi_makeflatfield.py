#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy

from podi_definitions import *

def normalize_flatfield(filename, outputfile, binning_x=8, binning_y=8, repeats=3, batchmode_hdu=None):

    if (not batchmode_hdu == None):
        hdulist = batchmode_hdu
    else:
        hdulist = pyfits.open(filename)
        
    filter = hdulist[0].header['FILTER']
    # print "This is filter",filter    # print "Using binning %d,%d" % (binning_x, binning_y)

    list_of_otas_to_normalize = otas_to_normalize_ff[filter]

    flatfield_data = numpy.zeros(shape=(13*4096*4096/(binning_x*binning_y)), dtype=numpy.float32)
    flatfield_data[:] = numpy.NaN

    datapos = 0
    for extension in range(1, len(hdulist)): #hdulist[1:]:
        fppos = int(hdulist[extension].header['FPPOS'][2:4])

        #print list_of_otas_to_normalize
        try:
            index = list_of_otas_to_normalize.index(fppos)
        except:
            # We didn't find this OTA in the list, so skip it
            hdulist[extension].header.update("FF_NORM", False, "Used in normalization")
            extension += 1
            continue

        hdulist[extension].header.update("FF_NORM", True, "Used in normalization")
        
        # We now know that we should include this OTA in the
        # calculation of the flat-field normalization
        stdout_write("\rAdding OTA %02d to flat-field ..." % fppos)
        #flatfield_data = numpy.concatenate((flatfield_data, extension.data.flatten()))
        #flatfield_data[extension,:,:] = extension.data

        if (binning_x>1 or binning_y>1):
            sx, sy = hdulist[extension].data.shape[0], hdulist[extension].data.shape[1]
            bx, by = sx/binning_x, sy/binning_y
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
    stdout_write(" computing median ...")
    sigma_min, sigma_max = -1e5, 1e6
    for i in range(repeats):

        valid = (flatfield_data > sigma_min) & (flatfield_data < sigma_max)

        ff_median_level = numpy.median(flatfield_data[valid])
        ff_std = numpy.std(flatfield_data[valid])

        sigma_min = ff_median_level - 2 * ff_std
        sigma_max = ff_median_level + 3 * ff_std
        #print i, numpy.sum(valid), datapos, ff_median_level, ff_std, sigma_min, sigma_max
        
    if (ff_median_level <= 0):
        print "Something went wrong or this is no flatfield frame"
        ff_median_level = 1.0

    stdout_write("\b\b\b(% 7.1f) ..." % (ff_median_level))
    
    # Now normalize all OTAs with the median flatfield level
    stdout_write(" normalizing ...")
    for extension in range(1, len(hdulist)):
        hdulist[extension].data /= ff_median_level
        hdulist[extension].data[hdulist[extension].data < 0.1] = numpy.NaN
        hdulist[extension].header.add_history("FF-level: %.1f" % (ff_median_level))
        
    stdout_write(" writing results ...")
    if (os.path.isfile(outputfile)):
	os.remove(outputfile)
    hdulist.writeto(outputfile, clobber=True)
    stdout_write(" done!\n")
       
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
            print filename, outputfile            

            # And finally, do the actual work
            normalize_flatfield(filename, outputfile, binning_x=binning_x, binning_y=binning_y, repeats=repeats)
    else:
        filename = clean_list[1]
        outputfile = clean_list[2]
        normalize_flatfield(filename, outputfile, binning_x=binning_x, binning_y=binning_y, repeats=repeats)

