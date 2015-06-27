#!/usr/bin/env python
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
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
import scipy
import ephem

gain_correct_frames = False
from podi_definitions import *

def scamp_header_to_minifits(filename, minifits_outputname, reference_fits, 
                             reference_extension="OTA33.SCI"):

    if (not os.path.isfile(filename)):
        return None

    headfile = open(filename, "r")
    lines = headfile.readlines()

    values = {}
    values_list = []
    comments = {}
    comments_list = []

    extlist = []
    primhdu = pyfits.PrimaryHDU()
    extlist.append(primhdu)
    
    for line in lines:
        #print line

        key, value, comment = line[0:8].strip(), line[9:30].strip(), line[32:].strip()
        if (key in ("HISTORY", "COMMENT",
                    ) ):
            # Don't know what to do with those, so skip'em
            continue

        elif (key in ("FGROUPNO", "FLXSCALE", "MAGZEROP", 
                      "ASTINST",
                      "PHOTIRMS", "PHOTINST", "PHOTLINK",
                      ) ):
            # These are some scamp-specific headers, let's not copy them
            continue

        elif (key in ("CRVAL1", "CRVAL2",
                      "CRPIX1", "CRPIX2", 
                      "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                      "EQUINOX",
                      ) or (key[0:2] == "PV" and key[3] == "_")):
            value = float(value)

        elif (key in ("RADECSYS", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2") ):
            # Strip the unnecessary quotes and spaces from the end
            value = value[1:-1].strip()

        elif (key == "END"):
            # This concludes one extension, add it to list and start new 
            # list for the next OTA

            values_list.append(values)
            comments_list.append(comments)
            values = {}
            comments = {}
            continue

        values[key] = value
        comments[key] = comment

    refhdu = pyfits.open(reference_fits)
    # Count how many image extensions exist in the reference frame
    n_imageext = 0;
    for i in range(len(refhdu)):
        if (is_image_extension(refhdu[i])):
            n_imageext += 1
        # try:
        #     extname = refhdu[i].header['EXTNAME']
        #     if (extname[0:3] == "OTA" and extname[-3:] == "SCI"):
        #         n_imageext += 1
        # except:
        #     pass

    if (len(values_list) != n_imageext):
        stdout_write("Illegal scamp solution or wrong reference fits (%d vs %d)!\n" % (len(values_list), len(refhdu)-1))
        return -1

    # Now copy all scamp headers into the minifits
    # With this data at hand, work out the shift we need to apply to the scamp solution
    print "Creating the minfits headers"
    for ota in range(len(values_list)):

        # Create a new Image extension to hold the header
        hdu = pyfits.ImageHDU()

        # Copy all headers to the new HDU
        for key in values_list[ota]:
            hdu.header[key] = (values_list[ota][key], comments_list[ota][key])

        # And add the HDU to the list that will form the minifits
        extlist.append(hdu)


    # Now go through both the minifits and the reference and assign the extname keywords
    print "Correcting CRVALs for reference position"
    ota33_found = False
    for ext in range(1, len(refhdu)):
        if (not is_image_extension(refhdu[ext])):
            continue
        extname = refhdu[ext].header['EXTNAME']
        if (extname == reference_extension): ota33_found = True
        extlist[ext].name = extname

    # Create a proper HDUList from the list of new HDUs
    minifits_hdulist = pyfits.HDUList(extlist)


    # Next (and almost finally) we need to change the CRVAL keywords to account for the
    # offsets in pointing center and the center assumed by scamp

    # Adjust the solution for the known shift in CRPIX position
    # In pODI: True pointing of telescope is ~ at pixel 4200,4200 of OTA 3,3
    # From scamp: Typically somewhere around 2000,2000 in OTA 3,3
    # NB: 0.11 / 3600. is the pixel scale in degrees/pixel, ignoring rotation and distortion
    if (ota33_found):
        dx_crpix = 4200 - minifits_hdulist[reference_extension].header["CRPIX1"]
        dy_crpix = 4200 - minifits_hdulist[reference_extension].header["CRPIX2"]
        d_ra  = dx_crpix * 0.11 / 3600.
        d_dec = dy_crpix * 0.11 / 3600.
    else:
        d_ra, d_dec = 0, 0

    # fix potential problems with d_ra > 360 or d_ra < 0
    d_ra = d_ra - math.floor(d_ra/360.)*360

    for ext in range(1, len(minifits_hdulist)):
        minifits_hdulist[ext].header['CRVAL1'] = d_ra
        minifits_hdulist[ext].header['CRVAL2'] = d_dec

    # Last step: write the minifits to file
    minifits_hdulist.writeto(minifits_outputname, clobber=True)
    
    # That's it, all done!
    return 0


if __name__ == "__main__":

    scampfile = sys.argv[1]

    outputfile = sys.argv[2]

    reference = sys.argv[3]

    try:
        ref_ext = sys.argv[4]
    except:
        ref_ext = "OTA33.SCI"
    
    scamp_header_to_minifits(scampfile, outputfile, reference, ref_ext)
    
