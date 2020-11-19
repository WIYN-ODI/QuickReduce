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
import astropy.io.fits as pyfits
import numpy
import scipy

available_otas = [00, 16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]
available_ota_coords = [(0,0), (1,6), 
                        (2,2), (2,3), (2,4),
                        (3,2), (3,3), (3,4),
                        (4,2), (4,3), (4,4),
                        (5,5), (6,1),
                        ]

broken_ota_cells = [ [(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7)],  #00
                     [],#1.6
                     [],#2.2
                     [(6,6)],#2,3
                     [],#2,4
                     [],#3,2
                     [],#3,3
                     [(1,7),(3,1),],#3,4
                     [],#4,2
                     [],#4,3
                     [],#4,4
                     [],#55
                     [(1,3),(3,1)],#61
                     ]

#available_otas = [22, 23, 24, 32, 33, 34, 42, 43, 44]
#available_ota_coords = [(2,2), (2,3), (2,4),
#                        (3,2), (3,3), (3,4),
#                        (4,2), (4,3), (4,4),
#                        ]


#available_otas = [33]
#available_ota_coords = [(3,3)]

#available_otas = [33]
   
wcs_headers = [
    "WCSASTRM",
    "WCSAXES",
    "CTYPE1",
    "CTYPE2",
    "CRVAL1",
    "CRVAL2",
    "CRPIX1",
    "CRPIX2",
    "CD1_1",
    "CD2_1",
    "CD1_2",
    "CD2_2",
    ]

def break_region_string(str_region):
    reg = str_region[1:-1]
    x,dummy,y = reg.partition(",")
    x1,dummy,x2 = x.partition(":")
    y1,dummy,y2 = y.partition(":")
    return int(x1)-1, int(x2)-1, int(y1)-1, int(y2)-1

def extract_region(data, str_region):
    x1,x2,y1,y2 = break_region_string(str_region)
    return data[y1:y2+1, x1:x2+1]


def insert_into_array(data, from_region, target, target_region):

    fx1, fx2, fy1, fy2 = break_region_string(from_region)
    tx1, tx2, ty1, ty2 = break_region_string(target_region)

    if (fx2-fx1 != tx2-tx1 or fy2-fy1 != ty2-ty1):
        print("Dimensions do not match, doing nothing")
    else:
        target[ty1:ty2+1, tx1:tx2+1] = data[fy1:fy2+1, fx1:fx2+1]

    return 0



if __name__ == "__main__":


    maskname = sys.argv[1]
    mask_hdulist = pyfits.open(maskname)

    fppos_list = []
    print("Loading mask ...",)
    for hdu in mask_hdulist[1:]:
        #print hdu.header['FPPOS']
        fppos_list.append(hdu.header['FPPOS'])
    print("done!")

    for filename in sys.argv[2:]:
        print("Opening file",filename,)
        hdulist = pyfits.open(filename)

        fppos = hdulist[0].header['FPPOS']
        
        try:
            index = fppos_list.index(fppos)
        except:
            continue

        print("found index",index,)

        # Open file
        data = hdulist[0].data

        # Also open the original mask
        mask = mask_hdulist[index+1].data

        mask_match = mask[0:data.shape[0], 0:data.shape[1]]

        data[mask_match == 0] = numpy.NaN

        output_filename = filename[0:-5]+"masked.fits"
        hdulist[0].data = data
        hdulist.writeto(output_filename, overwrite=True)

        hdulist.close()
        del hdulist
        print("done!")

    sys.exit(0)
