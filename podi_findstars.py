#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy

from podi_definitions import *

boxsize = 24

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
        
    #corner_level[0] = numpy.mean(data[-corner_size:, -corner_size:])
    #corner_level[1] = numpy.mean(data[-corner_size:,  :corner_size])
    #corner_level[2] = numpy.mean(data[:corner_size ,  :corner_size])
    #corner_level[3] = numpy.mean(data[:corner_size , -corner_size:])
    bg_level = numpy.max([0, numpy.median(corner_level)])

    #corner_variance = numpy.zeros(shape=(4))
    #corner_variance[0] = numpy.mean(numpy.pow((data[-corner_size:, -corner_size:] - corner_level[0]), 2))
    #corner_variance[1] = numpy.mean(numpy.pow((data[-corner_size:,  :corner_size] - corner_level[1]), 2))
    #corner_variance[2] = numpy.mean(numpy.pow((data[:corner_size ,  :corner_size] - corner_level[2]), 2))
    #corner_variance[3] = numpy.mean(numpy.pow((data[:corner_size , -corner_size:] - corner_level[3]), 2))
    bg_variance = numpy.mean(corner_variance)
    
    signal = data - bg_level
    #print data.shape, signal.shape

    minimumFlux = 3 * bg_variance
    
    total = signal.sum()
    significant_pixel = signal > minimumFlux
    total_above_min = signal[significant_pixel].sum()
    area = numpy.sum(significant_pixel)
    
    X, Y = numpy.indices(data.shape)
    #print X
    #print Y
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

    s_n = amplitude / bg_variance
    
    #print bg_level
    return amplitude, center_x, center_y, fwhm_x, fwhm_y, roundness, peak, bg_level, bg_variance, s_n, area
    
if __name__ == "__main__":

    inputfile = sys.argv[1]

    hdulist = pyfits.open(inputfile)

    binning = 4
    bboxsize = boxsize / binning
    
    for ext in range(1, len(hdulist)):
        
        data = hdulist[ext].data
        data[numpy.isnan(data)] = -1e10

        data_select = rebin_image(data, binning)
        
        # Determine some global background level
        global_bg = numpy.median(data_select[data_select > 0])
        global_bg_rms = math.sqrt(global_bg) #numpy.std(data_select[data_select > 0])
        print "#\n#\n# Found Bkacground %d +/- %d\n#\n#\n" % (global_bg, global_bg_rms)

        indices_x, indices_y = numpy.indices(data_select.shape)
        indices_x_1d = indices_x.flatten()
        indices_y_1d = indices_y.flatten()
        
        #print indices_x.shape, indices_x[0:5]
        outcat_file = "split_%02d.cat" % (ext)
        outcat = open(outcat_file, "w")

        nstars_found = 0
        peak_found = 1e9
        #nstars_found < 75 and 
        while(peak_found > global_bg+1*global_bg_rms):
            
            data_1d = data_select.ravel()
            #print data_1d.shape
            idx = numpy.argmax(data_1d)
            #print data_1d[idx]
            #print numpy.max(data_1d)

            x,y = indices_x_1d[idx], indices_y_1d[idx]
            peak_found = data_select[x,y]
            
            #print x, y, " --> ",data[x,y]

            full_x, full_y = x*binning, y*binning
            print "maxima = ",data_1d.max(), data_select[x,y], "   @ ",x, y, "  (",full_x, full_y,")"

            bx1, bx2 = max([0, x-bboxsize]), min([data_select.shape[0], x+bboxsize])
            by1, by2 = max([0, y-bboxsize]), min([data_select.shape[1], y+bboxsize])
            print "binned box   = %4d - %4d, %4d - %4d" % (bx1, bx2, by1, by2)
            
            fx1, fx2 = max([0, full_x-boxsize]), min([data.shape[0], full_x+boxsize])
            fy1, fy2 = max([0, full_y-boxsize]), min([data.shape[1], full_y+boxsize])
            print "full-res box = %4d - %4d, %4d - %4d" % (fx1, fx2, fy1, fy2)

            data_select[bx1:bx2, by1:by2] = -1e10
            
            #print lx1, lx2, ly1, ly2
            blocking_radius = boxsize
            center_x, center_y = full_x, full_y
            if (True): #data[x,y] < 60000):
                
                minibox = data[fx1:fx2, fy1:fy2]
                if (numpy.min(minibox) > -1e9):
                    amplitude, mb_center_x, mb_center_y, fwhm_x, fwhm_y, roundness, peak, bg_level, bg_variance, s_n, area = calculate_moments(minibox)
                    print "Found star %d at pos %02d, %02d   ===>  " % (nstars_found, full_x, full_y), #center_x, center_y),
                    

                    if (s_n < 2):
                        print "too faint, s/n=",s_n,"\n"
                    elif (math.fabs(roundness) > 0.3):
                        print "not round enough, %.3f vs %.3f" % (fwhm_x, fwhm_y),"\n"
                    elif (area < 10):
                        print "Insufficent significant area (%d < 10)" % (area)
                    else:
                        center_x, center_y = mb_center_x + fx1, mb_center_y + fy1
                        nstars_found += 1
                        #print amplitude, center_x, center_y, fwhm_x, fwhm_y, roundness, peak, bg_variance, s_n
                        print
                        print " -->  pos = %8.2f %8.2f" % (center_x, center_y)
                        print " --> fwhm = %8.2f %8.2f" % (fwhm_x, fwhm_y)
                        print " -->  S/N = %8.2f (noise=%6.1f)" % (s_n, bg_variance)
                        print " --> Peak = %6d (bg=%6d, total=%6d)" % (amplitude, bg_level, peak)
                        print
                        print >>outcat, center_y+1, center_x+1, fwhm_x, fwhm_y, s_n

                        blocking_radius = numpy.min([5 * math.sqrt(fwhm_x * fwhm_y), 3*boxsize])
                        #)amplitude, center_x, center_y, fwhm_x, fwhm_y, roundness, peak, bg_variance, s_n
                        # Now mask the area around this brightest pixel as invalid
                        #print numpy.array(data[x-4:x+5, y-4:y+5], dtype=numpy.int)
                        #print numpy.array(minibox, dtype=numpy.int)
                else:
                    print "Too close to a flagged region\n"
                    pass
            else:
                print "This is something saturated!\n"
                pass

            blocking_radius = numpy.max([blocking_radius, boxsize])
            
            # Now mask the area around this brightest pixel as invalid
            #print numpy.array(data[x-4:x+5, y-4:y+5], dtype=numpy.int)
            #print numpy.array(minibox, dtype=numpy.int)
            b_cenx, b_ceny = center_x / binning, center_y / binning
            b_blockrad = blocking_radius / binning
            _bx1, _bx2 = max([0, b_cenx-b_blockrad]), min([data_select.shape[0], b_cenx+b_blockrad])
            _by1, _by2 = max([0, b_ceny-b_blockrad]), min([data_select.shape[1], b_ceny+b_blockrad])
            #by1, by2 = max([0, y-bboxsize]), min([data_select.shape[1], y+bboxsize])

            print "BLOCKING OUT= %4d,%4d %4d,%4d" % (bx1, bx2, by1, by2)
            print "              %4d,%4d %4d,%4d" % (_bx1, _bx2, _by1, _by2)
            

            #data_select[bx1:bx2, by1:by2] = -1e10
            data_select[_bx1:_bx2, _by1:_by2] = -1e10

            
            
            #x, y = numpy.indices(data.shape)
            #print idx
            #print data.max(), data[x,y]

            #print
            
            #data_sortedindex = numpy.argsort(data)
            #print data_sortedindex[0:10,0:10]
        
        #outfile = "split_%02d.fits" % (ext)
        #prim = pyfits.PrimaryHDU(data=hdulist[ext].data, header=hdulist[ext].header)
        #out_hdulist = pyfits.HDUList([prim])
        #out_hdulist.writeto(outfile, clobber=True)
        
