#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy
import scipy.stats


from podi_definitions import *


def imcombine(filelist, outputfile):

    # For now assume we have 13 extensions, and all extensions are in the same order
    reference_filename = filelist[0]
    ref_hdulist = pyfits.open(reference_filename)
    filter = ref_hdulist[1].header['FILTER']

    primhdu = pyfits.PrimaryHDU()
    out_hdulist = pyfits.HDUList([primhdu])
    out_hdulist.writeto(outputfile, clobber=True)
    out_hdulist.close()
    del out_hdulist
    
    for cur_ext in range(1, len(ref_hdulist)):
        # Check what OTA we are dealing with
        ref_fppos = ref_hdulist[cur_ext].header['FPPOS']

        stdout_write("Working on OTA %s (#%d)" % (ref_fppos, cur_ext))

        # Allocate enough memory to load a single OTA from all files
        buffer = numpy.zeros(shape=(ref_hdulist[cur_ext].data.shape[0], ref_hdulist[cur_ext].data.shape[1], len(filelist)))
        # Set the full buffer to NaN
        buffer[:,:,:] = numpy.NaN
        print "buffer-size",buffer.shape
        
        # Copy the reference data
        buffer[:,:,0] = ref_hdulist[cur_ext].data[:,:]
        del ref_hdulist[cur_ext].data

        # Now open all the other files, look for the right extension, and copy their image data to buffer
        for file_number in range(1, len(filelist)):
            filename = filelist[file_number]
            hdulist = pyfits.open(filename)
            for i in range(1, len(hdulist)):
                fppos = hdulist[i].header['FPPOS']
                if (fppos == ref_fppos):
                    buffer[:,:,file_number] = hdulist[i].data[:,:]
                    break
            hdulist.close()
            del hdulist

        # Now compute median/average/sum/etc
        avg = numpy.mean(buffer, axis=2)
        print avg.shape
            
        # Create new ImageHDU
        hdu = pyfits.ImageHDU()

        # Insert the imcombine'd frame into the output HDU
        hdu.data = avg

        # Copy all headers from the reference HDU
        cards = ref_hdulist[cur_ext].header.ascardlist()
        for c in cards:
            hdu.header.update(c.key, c.value, c.comment)

        # Open the output file and append the new HDU
        out_hdulist = pyfits.open(outputfile, mode="update")
        out_hdulist.append(hdu)
        print "#OTAs in output:",len(out_hdulist)
        out_hdulist.flush()
        out_hdulist.close()
        del out_hdulist

        del hdu
        del buffer
            
    stdout_write(" done!\n")

if __name__ == "__main__":

    outputfile = sys.argv[1]

    filelist = sys.argv[2:]

    imcombine(filelist, outputfile)
