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

This module contains all data and most functionality to correct raw pODI
frames for the effects of cross-talk amongst cells.

For now, the crosstalk coefficient is defined globally, but the data structure
supports a more sophisticated system specifying the cross-talk coefficients
based on source and target-cell.

"""


import os
import sys
from podi_definitions import *
import pyfits
import scipy
import scipy.stats

x0 =  5.6E-5 #1.2e-4

xtalk_saturated_correction = 8
xtalk_saturation_limit = 65535



xtalk_coeffs = {

"OTA33.SCI": 
    ( 
        # Row 0:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 1:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 2:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 3:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 4:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 5:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 6:
            ((1,x0,x0,x0,x0,x0,x0,x0),
             (x0,1,x0,x0,x0,x0,x0,x0),
             (x0,x0,1,x0,x0,x0,x0,x0),
             (x0,x0,x0,1,x0,x0,x0,x0),
             (x0,x0,x0,x0,1,x0,x0,x0),
             (x0,x0,x0,x0,x0,1,x0,x0),
             (x0,x0,x0,x0,x0,x0,1,x0),
             (x0,x0,x0,x0,x0,x0,x0,1),
             ),
        # Row 7:
             ((1,x0,x0,x0,x0,x0,x0,x0),
              (x0,1,x0,x0,x0,x0,x0,x0),
              (x0,x0,1,x0,x0,x0,x0,x0),
              (x0,x0,x0,1,x0,x0,x0,x0),
              (x0,x0,x0,x0,1,x0,x0,x0),
              (x0,x0,x0,x0,x0,1,x0,x0),
              (x0,x0,x0,x0,x0,x0,1,x0),
              (x0,x0,x0,x0,x0,x0,x0,1),
              ),
     ), # end OTA33
}

xtalk_coeffs['OTA34.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA32.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA44.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA43.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA42.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA24.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA23.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA22.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA55.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA61.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA16.SCI'] = xtalk_coeffs['OTA33.SCI']
xtalk_coeffs['OTA00.SCI'] = xtalk_coeffs['OTA33.SCI']


xtalk_matrix = {}

import numpy

def invert_all_xtalk():
    """

    By default, the crosstalk coefficients are specified as the amplitude of the
    target cell. However, fto correct for cross-talk, we need to onvert the
    linear algebra system. This function handles all required steps.

    """
    
    global xtalk_matrix

    for ota in xtalk_coeffs:
        # print ota
        ota_matrices = []
        
        for row_matrix in xtalk_coeffs[ota]:
            #print row_matrix,"\n\n"
            
            inverted = numpy.linalg.inv(row_matrix)
            ota_matrices.append(inverted)
            
        xtalk_matrix[ota] = ota_matrices

invert_all_xtalk()


if __name__ == "__main__":
    print "Hello!"

    filename = sys.argv[1]
    hdulist = pyfits.open(filename)

    hdulist.info()

    overscan = numpy.ones((8,8)) * -1e9
    bglevel = numpy.ones((8,8)) * -1e9

    for row in range(8):
        for col in range(8):

            cell_name = "xy%d%d" % (col, row)

            data = extract_datasec_from_cell(hdulist[cell_name].data, binning=1)
            saturated = data >= 65535
            n_saturated = numpy.sum(saturated)

            if (n_saturated <= 10):
                # Nothing saturated in this cell, let's move on
                continue

            # We found at least a few saturated pixels
            for othercol in range(8):
                if (othercol == col):
                    # No crosstalk to the source cell
                    continue

                other_cellname = "xy%d%d" % (othercol, row)
                data = hdulist[other_cellname].data
                
                if (overscan[othercol,row] < 0):
                    # Compute overscan subtracted image
                    overscan_level = numpy.median(extract_biassec_from_cell(data, binning=1))
                    overscan[othercol,row] = overscan_level
                else:
                    overscan_level = overscan[othercol,row]

                image = extract_datasec_from_cell(data, binning=1) - overscan_level

                if (bglevel[othercol,row] < 0):
                    # Estimate the background level
                    bg_pixels = three_sigma_clip(image)
                    bgmedian = numpy.median(bg_pixels)
                    bgmode = 3*numpy.median(bg_pixels) - 2*numpy.mean(bg_pixels)
                    bglevel[othercol,row] = bgmode
                else:
                    bgmode = bglevel[othercol,row]

                image -= bgmode

                # Average the flux in the saturated pixels
                xtalk_affected = image[saturated]
                xtalk_cleaned, xtalk_mask = three_sigma_clip(xtalk_affected, return_mask=True)

                xtalk_effect = numpy.mean(xtalk_cleaned)

                print cell_name, other_cellname, overscan_level, bgmedian, numpy.mean(bg_pixels), bgmode, bgmode+overscan_level, xtalk_effect, numpy.sum(xtalk_mask), n_saturated

            

