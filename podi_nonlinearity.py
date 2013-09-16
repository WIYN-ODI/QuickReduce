#!/usr/local/bin/python
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
import scipy.optimize

from podi_plotting import *

gain_correct_frames = False
from podi_definitions import *



import podi_plotting




colors = ['', '#900000', '#00a000', '#0000a0']


bordersize = 75
max_polyorder=0


def create_nonlinearity_data(inputfiles):

    # Open all files, and loop over all cells in all extensions:

    all_data = []

    for filename in inputfiles:

        hdulist = pyfits.open(filename)
        exptime = hdulist[0].header['EXPTIME']
        expmeas = hdulist[0].header['EXPMEAS']

        for ext in range(1, len(hdulist)):
            if (not type(hdulist[ext]) == pyfits.hdu.image.ImageHDU):
                continue

            extname = hdulist[ext].header['EXTNAME']
            data = hdulist[ext].data
            ota = int(extname[3:5])
            otax = int(extname[3])
            otay = int(extname[4])

            for cellx in range(8):
                for celly in range(8):
                
                    cell_area = cell2ota__get_target_region(cellx, celly)
                    x1, x2, y1, y2 = cell_area

                    cell_data = data[y1:y2, x1:x2]

                    cell_center = cell_data[bordersize:-bordersize,bordersize:-bordersize]

                    median_int = numpy.median(cell_center)
                    mean_int = numpy.mean(cell_center)
                    std_int = numpy.std(cell_center)
                    #print exptime, extname, ota, cellx, celly, median_int

                    thiscell = [ota, otax, otay, cellx, celly, exptime, expmeas, median_int, mean_int, std_int]
                    print thiscell

                    all_data.append(thiscell)

        hdulist.close()

    return all_data



def create_nonlinearity_fits(data, outputfits, polyorder=3, exptime_ranges=[0.1,2.5]):

    otas = set(data[:,0])
    print otas

    result_count = 0 
    result_ota = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_cellx = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_celly = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_coeffs = numpy.zeros(shape=(data.shape[0],polyorder-1))
    result_coeffuncert = numpy.zeros(shape=(data.shape[0],polyorder-1))
    
    for ota in otas:
        for cellx in range(8):
            for celly in range(8):
                stdout_write("\rFitting OTA %02d, cell %1d,%1d ..." % (ota, cellx, celly))

                this_cell = (data[:,0] == ota) & (data[:,3] == cellx) & (data[:,4] == celly)

                subset = data[this_cell]
                #print ota, cellx, celly,":\n",subset

                not_nans = numpy.isfinite(subset[:,7]) & numpy.isfinite(subset[:,9])

                if (numpy.sum(not_nans) <= 0):
                    continue

                exptime = subset[:,5][not_nans]
                medlevel = subset[:,7][not_nans]
                stdlevel = subset[:,9][not_nans]

                def fit_fct(p, x):
                    y = numpy.zeros(x.shape)
                    for i in range(p.shape[0]):
                        y += p[i] * x**(i+1)
                    return y
                def err_fct(p,x,y,err, fitrange_x, fitrange_y):
                    yfit = fit_fct(p,x)
                    in_fit_range = numpy.isfinite(x) & numpy.isfinite(y)
                    if (not fitrange_x == None):
                        in_fit_range = in_fit_range & (x >= fitrange_x[0]) & (x <= fitrange_x[1])
                    if (not fitrange_y == None):
                        in_fit_range = in_fit_range & (y >= fitrange_y[0]) & (y <= fitrange_y[1])
                    if (err == None):
                        return ((y-yfit))[in_fit_range]
                    return ((y-yfit)/err)[in_fit_range]

                pinit = numpy.zeros(polyorder)

                # fit = scipy.optimize.leastsq(err_fct, pinit,
                #                              args=(exptime, medlevel, stdlevel, exptime_ranges), 
                #                              full_output=1)

                fit = scipy.optimize.leastsq(err_fct, pinit,
                                             args=(medlevel, exptime, None, None, exptime_ranges), 
                                             full_output=1)

                pfit = fit[0]
                uncert = numpy.sqrt(numpy.diag(fit[1]))

                print ota, cellx, celly, pfit, uncert

                linear_factor = pfit[0]

                coefficients_normalized = (pfit / linear_factor)[1:]
                coefficient_errors_normalized = (uncert / linear_factor)[1:]
                
                result_ota[result_count] = ota
                result_cellx[result_count] = cellx
                result_celly[result_count] = celly
                result_coeffs[result_count] = coefficients_normalized
                result_coeffuncert[result_count] = coefficient_errors_normalized

                
                result_count += 1

    stdout_write(" done!\n")

    # Prepare all data to be written to a fits file
    print result_coeffs[:10,:]
    print result_coeffuncert[:10,:]

    if (outputfits == None):
        return

    columns = [\
        pyfits.Column(name='OTA',    format='B', array=result_ota[:], disp='ota'),
        pyfits.Column(name='CELLX',  format='B', array=result_cellx[:], disp='cell-x'),
        pyfits.Column(name='CELLY',  format='B', array=result_celly[:], disp='cell-y'),
        ]

    for i in range(polyorder-1):
        col = pyfits.Column(name="COEFF_X**%d" % (i+2), 
                            format='E',
                            array=result_coeffs[:,i],
                            disp="polynomial coeff x^%d" % (i+2)
                            )
        columns.append(col)
    for i in range(polyorder-1):
        col = pyfits.Column(name="UNCERT_COEFF_X**%d" % (i+2), 
                            format='E',
                            array=result_coeffuncert[:,i],
                            disp="uncertainty of polynomial coeff x^%d" % (i+2)
                            )
        columns.append(col)

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    primhdu = pyfits.PrimaryHDU()
    primhdu.header.update("POLYORDR", polyorder)

    hdulist = pyfits.HDUList([primhdu, tbhdu])
    hdulist.writeto(outputfits, clobber=True)
    
    return 



def load_nonlinearity_correction_table(filename, search_ota):
    
    # Load the catalog file
    hdulist = pyfits.open(filename)

    # Determine what the fittting order was
    polyorder = hdulist[0].header['POLYORDR']

    # Create an array holding all coefficients
    nonlinearity_coeffs = numpy.zeros(shape=(8,8,polyorder-1))

    # Now load the full catalog and sort the coefficients 
    # into the coefficient matrix
    ota = hdulist[1].data.field('OTA')
    cellx = hdulist[1].data.field('CELLX')
    celly = hdulist[1].data.field('CELLY')

    all_coeffs = numpy.zeros(shape=(ota.shape[0],polyorder-1))
    for order in range(polyorder-1):
        columnname = "COEFF_X**%d" % (order+2)
        all_coeffs[:,order] = hdulist[1].data.field(columnname)[:]

    in_this_ota = ota == search_ota

    cellx = cellx[in_this_ota]
    celly = celly[in_this_ota]
    all_coeffs = all_coeffs[in_this_ota]

    for i in range(cellx.shape[0]):
        nonlinearity_coeffs[cellx[i], celly[i], :] = all_coeffs[i]

    return nonlinearity_coeffs


def compute_cell_nonlinearity_correction(data, cellx, celly, all_coeffs):

    coeffs = all_coeffs[cellx, celly, :]
    #print cellx, celly, coeffs
    #print numpy.mean(data)
    
    correction = numpy.zeros(data.shape)
    for i in range(coeffs.shape[0]):
        correction += coeffs[i] * data**(i+2)

    #print numpy.mean(correction)
    return correction


if __name__ == "__main__":


    if (cmdline_arg_isset("-fit")):
        datafile = get_clean_cmdline()[1]
        data = numpy.loadtxt(datafile)

        create_nonlinearity_fits(data, "nonlinfit.fits")
        sys.exit(0)

    if (cmdline_arg_isset("-load")):
        fitsfile = get_clean_cmdline()[1]
        ota = int(get_clean_cmdline()[2])
        coeffs = load_nonlinearity_correction_table(fitsfile, ota)
        print coeffs
        sys.exit(0)

    if (cmdline_arg_isset("-correct")):
        infile = get_clean_cmdline()[1]
        hdulist = pyfits.open(infile)
        ota = int(hdulist[0].header['FPPOS'][2:4])
        catfile = get_clean_cmdline()[2]
        ota_coeffs = load_nonlinearity_correction_table(catfile, ota)
        
        for i in range(1, 65):
            stdout_write("\rcorrecting cell %d ..." % (i))
            data = hdulist[i].data
            overscan_level = numpy.median(data[:,494:])
            data -= overscan_level
            cellx = hdulist[i].header['WN_CELLX']
            celly = hdulist[i].header['WN_CELLY']
            correction = compute_cell_nonlinearity_correction(data, cellx, celly, ota_coeffs)

            #data += correction
            hdulist[i].data = data

        stdout_write(" writing ...")
        outfile = get_clean_cmdline()[3]
        hdulist.writeto(outfile, clobber=True)

        stdout_write(" done!\n\n")
        sys.exit(0)

    # Read the input directory that contains the individual OTA files
    inputfiles = get_clean_cmdline()[1:]
    all_data = create_nonlinearity_data(inputfiles)
    numpy.savetxt("nonlin.data", all_data)
