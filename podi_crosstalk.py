#!/usr/bin/env python3
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
import astropy.io.fits as pyfits
import scipy
import scipy.stats
import itertools

x0 =  5.6E-5 #1.2e-4

xtalk_saturated_correction = 8
xtalk_saturation_limit = 65535

import podi_focalplanelayout
import podi_reductionlog

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



def apply_old_crosstalk_correction(hdulist, fpl, extname2id):

    logger = logging.getLogger("CrossTalk")
    logger.debug("Using traditional crosstalk method")

    ota = int(hdulist[0].header['FPPOS'][2:])
    extname = "OTA%02d.SCI" % ota
    xtalk_matrix = fpl.get_crosstalk_matrix(extname)

    # Allocate some memory for the cells in one row
    xtalk_corr = [None] * 8

    for row in range(8):
        for column in range(8):

            # Allocate some more memory to hold the output of the cross-talk 
            # correction of this frame
            xtalk_corr[column] = numpy.zeros(hdulist[1].data.shape, dtype=numpy.float)

            for xtalk_column in range(8):
                # Construct the name of each cell that's causing the crosstalk
                xy_name = "xy%d%d" % (xtalk_column, row)

                # Now go through the list of all cells in this row and add them to 
                # the corrected cell content output
                #print "Adding ",xy_name,"to ",extname, column, row, "(scaling",podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column],")"

                correction = hdulist[extname2id[xy_name]].data * xtalk_matrix[row][xtalk_column][column]
                if (column != xtalk_column):
                    saturated = hdulist[extname2id[xy_name]].data >= fpl.crosstalk_saturation_limit(extname)
                    correction[saturated] = -1 * fpl.crosstalk_saturation_correction(extname)

                xtalk_corr[column] += correction #hdulist[xy_name].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                #print xtalk_corr[column][100,100]

        for column in range(8):
            # Now all cells in this row have been corrected, let's write them 
            # back into the hdulist so can can continue with the overscan subtraction etc.
            xy_name = "xy%d%d" % (column, row)
            hdulist[extname2id[xy_name]].data = xtalk_corr[column]

    return hdulist


def apply_crosstalk_correction(hdulist, xtalk_file, fpl, extname2id, options, reduction_log):

    logger = logging.getLogger("CrossTalk")

    if (xtalk_file is not None and not os.path.isfile(xtalk_file)):
        reduction_log.fail('crosstalk')
        return hdulist

    reduction_log.attempt('crosstalk')

    # Find the correct crosstalk coefficient file for this detector
    # ota = int(hdulist[0].header['FPPOS'][2:])
    # extname = "OTA%02d.SCI" % ota
    # xtalk_file = fpl.get_crosstalk_file(ota, options)
    logger.debug("XTalk-File: %s" % (xtalk_file))

    if (xtalk_file is None):
        res = apply_old_crosstalk_correction(hdulist, fpl, extname2id)
        reduction_log.success('crosstalk')
        return res

    xtalk_hdu = pyfits.open(xtalk_file)
    xtalk_inv = xtalk_hdu['XTALK.INV'].data
    #print xtalk_inv.shape
    logger.debug("Using new XTALK method")

    # Allocate some memory for the cells in one row
    xtalk_corr = [None] * 8

    for row in range(8):
        for column in range(8):

            # Allocate some more memory to hold the output of the cross-talk 
            # correction of this frame
            xtalk_corr[column] = numpy.zeros(hdulist[1].data.shape,
                                             dtype=numpy.float)

            for xtalk_column in range(8):
                # Construct the name of each cell that's causing the crosstalk
                xy_name = "xy%d%d" % (xtalk_column, row)

                # Now go through the list of all cells in this row and add them to 
                # the corrected cell content output
                #print "Adding ",xy_name,"to ",extname, column, row, "(scaling",podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column],")"

                correction = hdulist[extname2id[xy_name]].data.astype(numpy.float) * xtalk_inv[row, xtalk_column, column] #xtalk_matrix[row][xtalk_column][column]
                #if (column != xtalk_column):
                #    saturated = hdulist[extname2id[xy_name]].data >= fpl.crosstalk_saturation_limit(extname)
                #        correction[saturated] = -1 * fpl.crosstalk_saturation_correction(extname)

                xtalk_corr[column] += correction.astype(numpy.float) #hdulist[xy_name].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                    #print xtalk_corr[column][100,100]

        for column in range(8):
            # Now all cells in this row have been corrected, let's write them 
            # back into the hdulist so can can continue with the overscan subtraction etc.
            xy_name = "xy%d%d" % (column, row)
            hdulist[extname2id[xy_name]].data = xtalk_corr[column]

    reduction_log.success('crosstalk')
    return hdulist



def measure_crosstalk(filename):

    hdulist = pyfits.open(filename)

    #hdulist.info()

    overscan = numpy.ones((8,8)) * -1e9
    bglevel = numpy.ones((8,8)) * -1e9

    all_entries = []

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
                data = hdulist[other_cellname].data.astype(numpy.float32)
                
                if (overscan[othercol,row] < 0):
                    # Compute overscan subtracted image
                    overscan_level = numpy.median(extract_biassec_from_cell(data, binning=1))
                    overscan[othercol,row] = overscan_level
                else:
                    overscan_level = overscan[othercol,row]

                image = extract_datasec_from_cell(data, binning=1) - overscan_level

                if (bglevel[othercol,row] < 0 or True):
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
                xtalk_raw = numpy.mean(xtalk_affected)

                print(cell_name, other_cellname, overscan_level, bgmedian, \
                    numpy.mean(bg_pixels), bgmode, bgmode+overscan_level, \
                    xtalk_effect, numpy.sum(xtalk_mask), n_saturated, xtalk_ra)

                entry = [int(cell_name[2:]), int(other_cellname[2:]), xtalk_effect, numpy.sum(xtalk_mask), n_saturated]

                all_entries.append(entry)

    all_entries = numpy.array(all_entries)

    return all_entries


def complete_xtalk_matrix(data):

    # format is:
    # row, src_column, target_column
    matrix = numpy.empty((8,8,8))
    matrix[:,:,:] = numpy.NaN
    
    for src_x, row_y in itertools.product(range(8), repeat=2):
        src_cell = src_x * 10 + row_y

        # select all measurements for this source cell
        src = data[ data[:,0] == src_cell ]
        print(src)

        # Now fill in the details
        matrix[row_y, src_x, src_x] = 1.0  # every cell contributes 100% flux to itself

        for target_x in range(8):
            _xy = target_x * 10 + row_y
            target_cell = src[:,1] == _xy
            if (numpy.sum(target_cell) > 0):
                # computed weighted average in case there's more than one 
                # measurements for this cell combination
                print(numpy.sum(src[target_cell,2]), numpy.sum(src[target_cell,3]))
                                
                avg_xtalk_signal = \
                    numpy.sum(src[target_cell,2] * src[target_cell,3]) \
                    / numpy.sum(src[target_cell,3])
                matrix[row_y, src_x, target_x] = avg_xtalk_signal / 65535.

    #
    # Now we are likely left with a lot of holes
    # 
    for src_x, target_x in itertools.product(range(8), repeat=2):
        # find the ones we are missing
        no_data = numpy.isnan(matrix[:, src_x, target_x])
        fill_in = numpy.nanmean(matrix[~no_data, src_x, target_x])
        matrix[no_data, src_x, target_x] = fill_in

        
    #
    # Now that it's complete, invert it so we can easily use it later on
    #
    matrix_inverted = numpy.empty((8,8,8))
    for row_y in range(8):
        inverted = numpy.linalg.inv(matrix[row_y,:,:])
        #print inverted
        matrix_inverted[row_y,:,:] = inverted

    imghdu = pyfits.ImageHDU(data=matrix)
    imghdu.name  = "XTALK.FWD"
    imghdu2 = pyfits.ImageHDU(data=matrix_inverted)
    imghdu2.name  = "XTALK.INV"
    p_hdu = pyfits.PrimaryHDU()
    hdulist = pyfits.HDUList([p_hdu, imghdu, imghdu2])
    hdulist.writeto("matrix_inv.fits", overwrite=True)
        
    return matrix


if __name__ == "__main__":

    if (cmdline_arg_isset("-compute")):

        filename = get_clean_cmdline()[1]
        data = numpy.loadtxt(filename)

        matrix = complete_xtalk_matrix(data)
        numpy.savetxt("xtalk.matrix", matrix.reshape((-1,8)), fmt="%10.3e")
        pass

    elif (cmdline_arg_isset("-podi")):

        x = numpy.array(xtalk_coeffs['OTA33.SCI'])
        for y, sx, tx in itertools.product(range(8), repeat=3):
            src_cell = sx*10+y
            target_cell = tx*10+y
            if (src_cell == target_cell):
                continue
            print(src_cell, target_cell, x[y, sx, tx]*65535., 100, 100)

    elif (cmdline_arg_isset("-plot")):

        import matplotlib
        import matplotlib.pyplot

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)

        matrix = numpy.loadtxt("xtalk.matrix")
        img = ax.imshow(matrix[:8], interpolation='nearest', norm=matplotlib.colors.LogNorm())
        fig.colorbar(img)

        fig.show()
        matplotlib.pyplot.show()
        matplotlib.pyplot.close(fig)
        matplotlib.pyplot.close()

    else:
        all_entries = None

        for filename in get_clean_cmdline()[1:]:
            #filename = sys.argv[1]

            _ae  = measure_crosstalk(filename)
            all_entries = _ae if all_entries is None else numpy.append(all_entries, _ae, axis=0)

        numpy.savetxt(sys.stdout, all_entries, "%02d %02d %12.8f %5d %5d")
        numpy.savetxt("xtalk.test", all_entries)



