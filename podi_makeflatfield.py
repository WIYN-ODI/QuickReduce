#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy

from podi_definitions import *

def normalize_flatfield(filename, outputfile):
    hdulist = pyfits.open(filename)
    filter = hdulist[1].header['FILTER']
    print "This is filter",filter

    list_of_otas_to_normalize = which_otas_to_use[filter]

    flatfield_data = numpy.zeros(shape=(13*4100*4100), dtype=numpy.float32)
    flatfield_data[:] = numpy.NaN

    datapos = 0
    for extension in range(1, len(hdulist)): #hdulist[1:]:
        fppos = int(hdulist[extension].header['FPPOS'][2:4])
 	#print fppos

 	try:
	    index = list_of_otas_to_normalize.index(fppos)
	except:
	    # We didn't find this OTA in the list, so skip it
            hdulist[extension].header.update("FF_NORM", 0, "Used in normalization")
            extension += 1
	    continue

        hdulist[extension].header.update("FF_NORM", 1, "Used in normalization")
        
        # We now know that we should include this OTA in the
	# calculation of the flat-field normalization
	stdout_write("\rAdding OTA %02d to flat-field ..." % fppos)
	#flatfield_data = numpy.concatenate((flatfield_data, extension.data.flatten()))
	#flatfield_data[extension,:,:] = extension.data
        
        one_d = hdulist[extension].data.flatten()
        flatfield_data[datapos:datapos+one_d.shape[0]] = one_d

        datapos += one_d.shape[0]
        del one_d

    # Now we are through all flatfields, compute the median value
    stdout_write(" computing median ...")
    ff_median_level = numpy.median(flatfield_data[:datapos])
    if (ff_median_level <= 0):
        print "Something went wrong or this is no flatfield frame"
        ff_median_level = 1.0

    # Now normalize all OTAs with the median flatfield level
    stdout_write(" normalizing ...")
    for extension in range(1, len(hdulist)):
        hdulist[extension].data /= ff_median_level
        hdulist[extension].data[hdulist[extension].data < 0.2] = numpy.NaN

    stdout_write(" writing results ...")
    if (os.path.isfile(outputfile)):
	os.remove(outputfile)
    hdulist.writeto(outputfile, clobber=True)
    stdout_write(" done!\n")
       
if __name__ == "__main__":

    if (sys.argv[1] == "-multi"):
        for filename in sys.argv[2:]:
            outputfile = filename[:-5]+".norm.fits"
            print filename, outputfile            
            normalize_flatfield(filename, outputfile)
    else:
        filename = sys.argv[1]
        outputfile = sys.argv[2]
        normalize_flatfield(filename, outputfile)

