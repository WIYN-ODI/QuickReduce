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

if __name__ == "__main__":

#    for i in range(len(sys.argv)):
#        print i, sys.argv[i]

    output_base = sys.argv[1]
    
    dump = open(output_base+".txt", "w")
    dump2 = open(output_base+".fit", "w")

    # Create mask from saturated pixels in reference file.
    for filename in sys.argv[2:]:
        print filename,"..."
        hdulist = pyfits.open(filename)
        fppos = int(hdulist[0].header['FPPOS'][2:4])

        for ext in range(1, len(hdulist)):

            data = hdulist[ext].data

            saturated = (data >= 65535)

            overscan_level = numpy.median(data[494:, :480])
            data[494:, :480] -= overscan_level
            data[:,480:] -= overscan_level

            skylevel = numpy.median(data[:494, :480])
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

                lineprofile = data[max_y:,col]
                rownumbers = row_ids[max_y:,0]
                combined = numpy.empty(shape=(lineprofile.shape[0],10))
                avg = numpy.median(data[500:,col])
                combined[:,0] = rownumbers[:]
                combined[:,1] = rownumbers[:]-max_y
                combined[:,2] = lineprofile[:]
                combined[:,3] = numpy.sum(saturated[:,col])
                combined[:,4] = avg
                combined[:,5] = col
                combined[:,6] = xy
                combined[:,7] = fppos

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
                max_range = numpy.max(combined[:,1][combined[:,8] > skysigma])

                safety = 25 #2 * (max_y - min_y) + 25
                fitdata_x = combined[safety:max_range,1]
                fitdata_y = combined[safety:max_range,2]

                if (fitdata_x.shape[0] > 10 ):
                    fitdata_err = numpy.sqrt(fitdata_y + skylevel)
                    fitdata_err[fitdata_err <= 0] = numpy.sqrt(skylevel)

                    fitfunc = lambda p,x: p[0] * numpy.exp(-1 * x / p[1]) 
                    errfunc = lambda p,x,y,err: (y - fitfunc(p,x)) / err
                    pinit = [100, 75]
                    out = scipy.optimize.leastsq(errfunc, pinit, args=(fitdata_x, fitdata_y, fitdata_err), full_output=True)
                    pfinal = out[0]
                    covar = out[1]

                    try:
                        print >>dump2, xy, fppos, col, pfinal[0], pfinal[1], math.sqrt(covar[0,0]), math.sqrt(covar[1][1]), max_y, safety
                    except:
                        pass

                    combined[:,9] = fitfunc(pfinal, combined[:,1])

                numpy.savetxt(dump, combined)
                print >>dump, "\n\n\n\n\n"

    hdulist.writeto(output_base+".fits", clobber=True)

