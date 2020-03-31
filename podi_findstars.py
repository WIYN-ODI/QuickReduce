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

from __future__ import print_function
import sys
import os
import astropy.io.fits as pyfits
import numpy
import scipy
#import pywcs
from astLib import astWCS

from podi_definitions import *

FWHM = 2 *math.sqrt(2*math.log(2))

def calculate_moments(data):

    # Corners are as follows:
    #
    #  3 -------- 0
    #  |          |
    #  |          |
    #  |          |
    #  2 -------- 1
    #
    #

    corner = [None]*4
    corner_size = 3
    corner[0] = data[-corner_size:, -corner_size:]
    corner[1] = data[-corner_size:,  :corner_size]
    corner[2] = data[:corner_size ,  :corner_size]
    corner[3] = data[:corner_size , -corner_size:]

    # First determine some background level, based on the mean value of the 4 corners
    corner_level = numpy.zeros(shape=(4))
    corner_variance = numpy.zeros(shape=(4))
    for i in range(4):
        corner_level[i] = numpy.mean(corner[i])
        corner_variance[0] = numpy.mean(numpy.power((corner[i] - corner_level[0]), 2))
        
    bg_level = numpy.max([0, numpy.median(corner_level)])
    bg_variance = numpy.mean(corner_variance)

    # Subtract the background
    signal = data - bg_level

    # Some stats about pixels above the detection threshold
    minimumFlux = 3 * bg_variance
    total = signal.sum()
    significant_pixel = signal > minimumFlux
    total_above_min = signal[significant_pixel].sum()
    area = numpy.sum(significant_pixel)

    # Compute the moments of the distribution, i.e. the x/y fwhms
    X, Y = numpy.indices(data.shape)
    center_x = (X*signal).sum()/total
    center_y = (Y*signal).sum()/total

    dx = X - center_x
    dy = Y - center_y
    if (total_above_min > 0):
        momentXX = numpy.sum((signal * dx * dx)[significant_pixel]) / total_above_min
        momentYY = numpy.sum((signal * dy * dy)[significant_pixel]) / total_above_min
        momentXY = numpy.sum((signal * dx * dy)[significant_pixel]) / total_above_min

        fwhm_x = numpy.sqrt(momentXX) * FWHM
        fwhm_y = numpy.sqrt(momentYY) * FWHM
        roundness = (fwhm_x - fwhm_y) / (fwhm_x + fwhm_y)
    else:
        fwhm_x, fwhm_y, roundness = 0, 0, -5
        
    #col = data[:, int(y)]
    #print col
    #width_x = numpy.sqrt(abs((numpy.arange(col.size)-y)**2*col).sum()/col.sum())
    #row = data[int(x), :]
    #width_y = numpy.sqrt(abs((numpy.arange(row.size)-x)**2*row).sum()/row.sum())

    peak = data.max()
    amplitude = peak - bg_level

    # Important: signal/noise
    s_n = amplitude / bg_variance
    
    return amplitude, center_x, center_y, fwhm_x, fwhm_y, roundness, peak, bg_level, bg_variance, s_n, area
    






def find_stars(hdu,
               binning=4,
               boxsize=24,
               dumpfile="find_stars.dump",
               verbose=True,
               detect_threshold=2,
               detect_minarea=5,
               saturation_limit=60000,
               roundness_limit=[-1,1],
               max_starcount=1e9,
               extension_id=0
               ):


    # Extract data and WCS solution from HDU
    data = hdu.data.copy().transpose()
    wcs = astWCS.WCS(hdu.header, mode="pyfits")
    
    # Prepare the data: Set all NaNs to something illegal
    data[numpy.isnan(data)] = -1e10

    # Create the binned array we will use to search for maxima
    binned = rebin_image(data, binning)
    bboxsize = boxsize / binning

    # Determine some global background level
    global_bg = numpy.median(binned[binned > 0])
    global_bg_rms = math.sqrt(global_bg) 
    print("# Found Background %d +/- %d" % (global_bg, global_bg_rms))

    # Prepare some indices so we can reconstruct the pixel position of each peak
    indices_x, indices_y = numpy.indices(binned.shape)
    indices_x_1d = indices_x.ravel()
    indices_y_1d = indices_y.ravel()

    # Create the output file that will hold the star catalog
    if (dumpfile != None):
        outcat_file = dumpfile 
        outcat = open(outcat_file, "w")

    #
    # Begin the actual work
    #

    # Sort the binned image by brightness
    binned_1d = binned.ravel()
    sorted_1d = numpy.argsort(binned_1d)

    
    nstars_found = 0
    npeaks_found = 0
    peak_value = 1e9
    
    current_peak = 0   # Holds where we are in the sorted list

    source_list = []

    min_peak_level = global_bg+detect_threshold*global_bg_rms
    #print "min peak level=",min_peak_level

    # Continue searching for peaks as long there's a chance of finding a significant one
    while(peak_value > min_peak_level and len(source_list)<max_starcount):
        current_peak += 1
        
        # Which pixel in the 1-d array are we talking about
        idx = sorted_1d[-1*current_peak]

        # Now determine the x,y positions of this pixel
        x,y = indices_x_1d[idx], indices_y_1d[idx]

        # Ignore this pixel if it's already marked as invalid
        if (binned_1d[idx] < -1e9 or binned[x,y] < -1e9):
            #current_peak += 1
            #print "skipping"
            continue

        peak_value = binned[x,y]
        if (verbose): print(current_peak, idx, binned_1d[idx])
        if (verbose): print("---> ",x, y, binned[x,y], binned_1d[idx])

        full_x, full_y = x*binning, y*binning
        if (verbose): print("maxima = ",binned_1d.max(), binned[x,y], "   @ ",x, y, "  (",full_x, full_y,")")

        bx1 = 0 if x-bboxsize < 0 else x-bboxsize
        bx2 = binned.shape[0] if x+bboxsize > binned.shape[0] else x+bboxsize
        by1 = 0 if y-bboxsize < 0 else y-bboxsize
        by2 = binned.shape[1] if y+bboxsize > binned.shape[1] else y+bboxsize
        #bx1, bx2 = max([0, x-bboxsize]), min([binned.shape[0], x+bboxsize])
        #by1, by2 = max([0, y-bboxsize]), min([binned.shape[1], y+bboxsize])
        if (verbose): print("binned box   = %4d - %4d, %4d - %4d" % (bx1, bx2, by1, by2))

        fx1 = 0 if full_x-boxsize < 0 else full_x-boxsize
        fx2 = data.shape[0] if full_x+boxsize > data.shape[0] else full_x+boxsize
        fy1 = 0 if full_y-boxsize < 0 else full_y-boxsize
        fy2 = data.shape[1] if full_y+boxsize > data.shape[1] else full_y+boxsize
        #fx1, fx2 = max([0, full_x-boxsize]), min([data.shape[0], full_x+boxsize])
        #fy1, fy2 = max([0, full_y-boxsize]), min([data.shape[1], full_y+boxsize])
        if (verbose): print("full-res box = %4d - %4d, %4d - %4d" % (fx1, fx2, fy1, fy2))

        binned[bx1:bx2, by1:by2] = -1e10

        blocking_radius = boxsize
        center_x, center_y = full_x, full_y
        npeaks_found += 1

        if (data[x,y] < saturation_limit):

            minibox = data[fx1:fx2, fy1:fy2]
            if (numpy.min(minibox) > -1e9):
                amplitude, mb_center_x, mb_center_y, fwhm_x, fwhm_y, roundness, peak, bg_level, bg_variance, s_n, area = calculate_moments(minibox)
                if (verbose): print("Found star %d at pos %02d, %02d   ===>  " % (nstars_found, full_x, full_y),) #center_x, center_y),

                if (s_n < detect_threshold):
                    if (verbose): print("too faint, s/n=",s_n,"\n")
                elif (roundness < roundness_limit[0] or roundness > roundness_limit[1]):
                    if (verbose): print("not round enough, %.3f vs %.3f" % (fwhm_x, fwhm_y),"\n")
                elif (area < detect_threshold):
                    if (verbose): print("Insufficent significant area (%d < 10)" % (area))
                else:
                    center_x, center_y = mb_center_x + fx1, mb_center_y + fy1
                    nstars_found += 1
                    #print amplitude, center_x, center_y, fwhm_x, fwhm_y, roundness, peak, bg_variance, s_n
                    if (verbose):
                        print()
                        print(" -->  pos = %8.2f %8.2f" % (center_x, center_y))
                        print(" --> fwhm = %8.2f %8.2f" % (fwhm_x, fwhm_y))
                        print(" -->  S/N = %8.2f (noise=%6.1f)" % (s_n, bg_variance))
                        print(" --> Peak = %6d (bg=%6d, total=%6d)" % (amplitude, bg_level, peak))
                        print()
                        if (dumpfile is not None):
                            print(center_y+1, center_x+1, fwhm_x, fwhm_y, s_n, file=outcat)
                    stdout_write("\rFound %d stars (%d peaks, [%d >= %d])..." % (nstars_found, npeaks_found, peak_value, min_peak_level))
                    blocking_radius = numpy.min([5 * math.sqrt(fwhm_x * fwhm_y), 3*boxsize])

                    # Add 1 pixel, as fits starts counting at 1, whereas numpy starts at 0
                    center_x += 1
                    center_y += 1
                    ra, dec = wcs.pix2wcs(center_x, center_y)
                    source_info = [ ra, dec, center_x, center_y, fwhm_x, fwhm_y, amplitude, peak, bg_level, bg_variance, s_n, area, extension_id]
                    source_list.append(source_info)
                del minibox
            else:
                if (verbose): print("Too close to a flagged region\n")
                pass
        else:
            if (verbose): print("This is something saturated!\n")
            pass

        blocking_radius = numpy.max([blocking_radius, boxsize])

        b_cenx, b_ceny = center_x / binning, center_y / binning
        b_blockrad = blocking_radius / binning
        _bx1 = 0 if b_cenx-b_blockrad < 0 else b_cenx-b_blockrad
        _bx2 = binned.shape[0] if b_cenx+b_blockrad > binned.shape[0] else b_cenx+b_blockrad
        _by1 = 0 if b_ceny-b_blockrad < 0 else b_ceny-b_blockrad
        _by2 = binned.shape[1] if b_ceny+b_blockrad > binned.shape[1] else y+bboxsize

        #if (verbose):
        #    print "BLOCKING OUT= %4d,%4d %4d,%4d" % (bx1, bx2, by1, by2)
        #    print "              %4d,%4d %4d,%4d" % (_bx1, _bx2, _by1, _by2)

        binned[_bx1:_bx2, _by1:_by2] = -1e10


    #outfile = dumpfile+".fits"
    #clobberfile(outfile)
    #prim = pyfits.PrimaryHDU(data=data, header=None)
    #out_hdulist = pyfits.HDUList([prim])
    #out_hdulist.writeto(outfile, overwrite=True)

    source_cat = numpy.array(source_list)

    stdout_write("done!\n")
    return source_cat


def findstar_allota(inputfile):
    hdulist = pyfits.open(inputfile)
    for ext in range(1, len(hdulist)):
        dumpfile = "fs_dump.%s" % (hdulist[ext].header['EXTNAME'])
        data = hdulist[ext].data
        find_stars(data, binning=4, boxsize=24, dumpfile=dumpfile, verbose=False)
    
if __name__ == "__main__":

    inputfile = sys.argv[1]

    hdulist = pyfits.open(inputfile)

    for ext in range(1, len(hdulist)):
        print("\n\n\n\n\n%s\n\n\n\n\n" % (hdulist[ext].header['EXTNAME']))

        dumpfilename = "fs_dump.%s" % (hdulist[ext].header['EXTNAME'])
        dumpfile = open(dumpfilename, "w")
        
        source_cat = find_stars(hdulist[ext], binning=4, boxsize=24, dumpfile=dumpfilename, verbose=False,
                                detect_threshold=1.5, detect_minarea=6, roundness_limit=[-0.2,+0.2])

        #ra, dec = pix2sky(source_cat[:,0], source_cat[:,1], hdulist[ext].header)

        #print source_cat[:,0:2]
        
        #wcs = pywcs.WCS(header=hdulist[ext].header)
        #ra_dec = wcs.all_pix2sky(source_cat[:,0:2], 0)
        #xxx.fits")
        #print ra_dec

        numpy.savetxt(dumpfile, source_cat[:,0:4], delimiter=" ")

        dumpfile.close()
        
    #import cProfile, pstats
    #cProfile.run('findstar_allota(inputfile)', "profiler")
    #p = pstats.Stats("profiler")
    #p.strip_dirs().sort_stats('time').print_stats()
    #p.sort_stats('time').print_stats()
