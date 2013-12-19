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



import sys
import os
import pyfits
import numpy
import scipy
import scipy.optimize
import time
import podi_persistency
import math


def gaussian(x, pos, ampl, width):
    return ampl * exp(-(x-pos)**2/(2*width**2))

if __name__ == "__main__":

#    for i in range(len(sys.argv)):
#        print i, sys.argv[i]

    output_base = sys.argv[1]
    
    dump = open(output_base+".txt", "w")
    dump2 = open(output_base+".fit", "w")

    # Create mask from saturated pixels in reference file.
    filename = sys.argv[2]
    allinfo = []
    if (True):
        print filename,"..."
        hdulist = pyfits.open(filename)
        fppos = int(hdulist[0].header['FPPOS'][2:4])
        #ctimestamp = time.ctime(os.path.getctime(filename))
        mtimestamp = os.path.getmtime(filename)
        #print ctimestamp
        #print mtimestamp

        for ext in range(1, len(hdulist)):

            data = hdulist[ext].data

            saturated = (data >= 65535)

            overscan_level = numpy.median(data[494:, :480] * 1.0)
            data[494:, :480] -= overscan_level
            data[:,480:] -= overscan_level

            skylevel = numpy.median(data[:494, :480] * 1.0)
            data[:494, :480] -= skylevel

            hdulist[ext].data = data

            # Now extract all columns with saturated pixels and write them to file for later plotting

            rows, cols = numpy.indices(data.shape)
            row_ids, col_ids = numpy.indices((data.shape[0],1))

            mask_thisframe = numpy.zeros(shape=data.shape)
            mask_thisframe = mask_thisframe > 1
            saturated_rows = rows[saturated]
            saturated_cols = cols[saturated]

            unique_cols = set(saturated_cols)
            xy = int(hdulist[ext].header['EXTNAME'][2:4])

            for col in unique_cols:
                this_col_saturated = row_ids[saturated[:,col]]
                
                min_y = numpy.min(this_col_saturated)
                max_y = numpy.max(this_col_saturated)

                # extract the full column from the data frame
                lineprofile = data[:,col]
                rownumbers = row_ids[:,0]

                combined = numpy.empty(shape=(lineprofile.shape[0],11))
                avg = numpy.median(data[500:,col])
                n_saturated = numpy.sum(saturated[:,col])
                combined[:,0] = rownumbers[:]
                combined[:,1] = rownumbers[:]-max_y+(n_saturated/2)
                combined[:,2] = lineprofile[:]
                combined[:,3] = n_saturated
                combined[:,4] = avg
                combined[:,5] = col
                combined[:,6] = xy
                combined[:,7] = fppos


                numpy.savetxt(dump, combined)
                print >>dump, "\n\n\n\n\n"

                continue


                #s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
                # Smooth the array in 10 pixels boxes
                # Fit over range where signal > 1 * sigma_sky
                window_len = 10
                if (lineprofile.shape[0] <= window_len):
                    continue

                if (True):
                    w=numpy.ones(10,'d')
                else:
                    w=eval('numpy.'+window+'(window_len)')
                combined[:,8] = numpy.convolve(w/w.sum(),combined[:,2],mode='same')
                skysigma = math.sqrt(skylevel)
                max_range = numpy.min(combined[:,1][combined[:,8] < skysigma])

                safety = numpy.max([max_range - 30, 0.5*max_range])
                #safety = 25 #2 * (max_y - min_y) + 25
                fitdata_x = combined[safety:max_range,1]
                fitdata_y = combined[safety:max_range,2]

                pfinal = [0,0]
                fit_a, fit_b = 0,0
                fit_a_err, fit_b_err = 0,0
                if (fitdata_x.shape[0] > 10):
                    fitdata_err = numpy.sqrt(fitdata_y + skylevel)
                    fitdata_err[fitdata_err < skysigma] = skysigma

                    fitfunc = lambda p,x: p[0] * numpy.exp(-1 * x / p[1]) 
                    errfunc = lambda p,x,y,err: (y - fitfunc(p,x)) / err
                    pinit = [100, 75]
                    out = scipy.optimize.leastsq(errfunc, pinit, args=(fitdata_x, fitdata_y, fitdata_err), full_output=True)
                    pfinal = out[0]
                    covar = out[1]

                    fit_a, fit_b = pfinal[0], pfinal[1]
                    fit_a_err, fit_b_err = math.sqrt(covar[0,0]), math.sqrt(covar[1,1])

                    try:
                        print >>dump2, xy, fppos, col, fit_a, fit_b, fit_a_err, fit_b_err, max_y, safety, n_saturated, skysigma
                    except:
                        pass

                    combined[:,9] = fitfunc(pfinal, combined[:,1])
                combined[:,10] = skysigma

                thisinfo = (fppos, xy, combined, col, 
                            max_y, min_y, numpy.sum(saturated[:,col]), 
                            fit_a, fit_a_err, fit_b, fit_b_err,
                            mtimestamp, safety, max_range,
                            )
                allinfo.append(thisinfo)

                numpy.savetxt(dump, combined)
                print >>dump, "\n\n\n\n\n"
                
    sys.exit(0)

    hdulist.writeto(output_base+".fits", clobber=True)

    persfile = sys.argv[3]
    hdulist = pyfits.open(persfile)
    
    perstime = os.path.getmtime(persfile)

    dumpx = open(output_base+".pers", "w")
    for ext in range(1, len(hdulist)):
        data = numpy.array(hdulist[ext].data, dtype=numpy.float32)

        overscan_level = numpy.median(data[494:, :480])
        data[494:, :480] -= overscan_level
        data[:,480:] -= overscan_level

        skylevel = numpy.median(data[:494, :480])
        data[:494, :480] -= skylevel

        hdulist[ext].data = data
        xy = int(hdulist[ext].header['EXTNAME'][2:4])

        

        for i in range(len(allinfo)):
            #fit_a, fit_a_err, fit_b, fit_b_err
            (fppos, sat_xy, combined, col, 
             max_y, min_y, n_saturated, 
             fit_a, fit_a_err, fit_b, fit_b_err, 
             mtimestamp, safety, max_range
             ) = allinfo[i]
            if (not xy == sat_xy):
                continue

            # Now measure how much persistency we have:
            median = numpy.median(data[:min_y,col])
            mean = numpy.mean(data[:min_y,col])
            d_time = perstime - mtimestamp

            print >> dumpx, perstime, mtimestamp, d_time, median, mean, fit_a, fit_b, fit_a_err, fit_b_err, n_saturated, col


