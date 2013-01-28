#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import ephem

from podi_definitions import *


def add_circle(buffer, center_x, center_y, radius, amplitude, fuzzy):

    x, y = numpy.indices(buffer.shape)

    dx = x - center_x - buffer.shape[0]/2
    dy = y - center_y - buffer.shape[1]/2
    d2 = dx*dx + dy*dy

    print dx[0:10,0:10]
    print dy[0:10,0:10]
    print d2[0:10,0:10]

    print d2[995:1005, 995:1005]
    
    tmp_buffer = numpy.zeros(shape=buffer.shape)
    tmp_buffer[d2 < radius*radius] = amplitude

    buffer += tmp_buffer

    return

    
if __name__ == "__main__":

    # Read in the input parameters
    command_file = sys.argv[1]
    outputfile = sys.argv[2]

    # Load command file 
    cmdfile = open(command_file, "r")
    cmds = cmdfile.readlines()
    
    pixelscale = 0.5
    dimension_arcsec = 1000

    # Create the pixel array
    buffer = numpy.zeros(shape=(dimension_arcsec/pixelscale, dimension_arcsec/pixelscale))
    print buffer.shape


    for i in range(len(cmds)):
        line = cmds[i]
        
        items = line.strip().split()
        print items

        shape = items[0]
        if (shape == "filledcircle"):
            coord_j2000 = ephem.Equatorial(items[1], items[2], epoch=ephem.J2000)
            center_x = numpy.degrees(coord_j2000.ra)*3600./pixelscale
            center_y = numpy.degrees(coord_j2000.ra)*3600./pixelscale
            radius = float(items[3])/pixelscale
            amplitude = float(items[4])
            fuzzy = float(items[5])

            if (i==0):
                refx = numpy.degrees(coord_j2000.ra)
                refy = numpy.degrees(coord_j2000.dec)

            center_x -= refx*3600/pixelscale
            center_y -= refy*3600/pixelscale
            
            add_circle(buffer, center_x, center_y, radius, amplitude, fuzzy)

    primhdu = pyfits.PrimaryHDU(data=buffer)
    primhdu.header.update("CTYPE1", "RA---TAN")
    primhdu.header.update("CTYPE2", "DEC--TAN")
    primhdu.header.update("CRVAL1", refx)
    primhdu.header.update("CRVAL2", refy)
    primhdu.header.update("CRPIX1", buffer.shape[0]/2)
    primhdu.header.update("CRPIX2", buffer.shape[1]/2)
    primhdu.header.update("CD1_1", pixelscale/3600)
    primhdu.header.update("CD2_2", pixelscale/3600)

    primhdu.writeto(outputfile, clobber=True)
    
    sys.exit(0)    
    
    
    

    
    op = sys.argv[2]
    input_2 = sys.argv[3]
    output = sys.argv[4]

    stdout_write("\nOpening input files ...")
    # Open both input fits files
    hdu_1 = pyfits.open(input_1)
    hdu_2 = pyfits.open(input_2)
    stdout_write(" done!\n")

    rebin_fac = int(cmdline_arg_set_or_default("-bin", 1))
    
    # Now go though each extension and perform the operation
    for idx_img1 in range(1, len(hdu_1)):
        img1 = hdu_1[idx_img1]
        
        fppos1 = img1.header['FPPOS']
        stdout_write("\rComputing extension %s (%2d of %2d) ..." % (img1.header['EXTNAME'], idx_img1, len(hdu_1)-1))
        for img2 in hdu_2[1:]:
            fppos2 = img2.header['FPPOS']
            if (fppos2 == fppos1):
                # This is the one

                if (op == "+"):
                    img1.data += img2.data

                elif (op == "-"):
                    img1.data -= img2.data

                elif (op == "/"):
                    img1.data /= img2.data

                elif (op == "x"):
                    img1 *= img2.data

                else:
                    stdout_write("Unkwnon operation %s\n" % (op))

        if (rebin_fac > 1):
            img1.data = rebin_image(img1.data, rebin_fac)

    # Now all the work is done, all final data is stored in img2, write results to new file                    
    stdout_write(" writing output ...")
    clobberfile(output)
    hdu_1.writeto(output)
    stdout_write(" done!\n\n")

