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

if __name__ == "__main__":

#    for i in range(len(sys.argv)):
#        print i, sys.argv[i]

    reference_file = sys.argv[1]
    extension = sys.argv[2]

    output_txt = sys.argv[3]

    # Create mask from saturated pixels in reference file.

    hdu = pyfits.open(reference_file)
    data = hdu[extension].data
    #hdu.info()

    #filetime = time.strptime(reference_file[0:15],"%Y%m%dT%H%M%S")

    overscan_level = numpy.median(data[:,500:])
    print "overscan-level=",overscan_level
                                  
    obsid = hdu[0].header['OBSID']
    ref_mjd = podi_persistency.get_mjd_from_timestamp(obsid[0:15])
    #ref_mjd = podi_persistency.get_mjd_from_timestamp(reference_file[0:15])
    #hdu[0].header['MJD-OBS']

    mask = (data >= 65500) #65535)
    saturated = (data >= 65535)
    data -= overscan_level

    bglevel = numpy.median(data[10:480,10:490])
    data -= bglevel

    hdu.close()
    del  hdu

    rows, cols = numpy.indices(data.shape)

    mask_thisframe = numpy.zeros(shape=data.shape)
    mask_thisframe = mask_thisframe > 1
    mask_time      = numpy.zeros(shape=data.shape)
    mask_time      = mask_time > 1
        #mask_time.fill(False)

    saturated_rows = rows[saturated]
    saturated_cols = cols[saturated]

    unique_cols = set(saturated_cols)

    dump1 = open(output_txt+".sat", "w")

    # New, optimized and way faster method
    row_ids, col_ids = numpy.indices((data.shape[0],1))

    satcolumn_info = []

    overscan_sky = numpy.median(data[500:,5:470])
    data[495:,0:480] -= overscan_sky

    for col in unique_cols:
        #print "working on col",col #saturated[col,:]
        this_col_saturated = row_ids[saturated[:,col]]
        #print "saturated in this col",this_col_saturated
        min_y = numpy.min(this_col_saturated)
        max_y = numpy.max(this_col_saturated)

        mask_thisframe[min_y:, col] = True
        mask_time[:max_y, col] = True

        lineprofile = data[max_y:,col]
        rownumbers = row_ids[max_y:,0]
        combined = numpy.empty(shape=(lineprofile.shape[0],6))
        avg = numpy.median(data[500:,col])
        combined[:,0] = rownumbers[:]
        combined[:,1] = rownumbers[:]-max_y
        combined[:,2] = lineprofile[:]
        combined[:,3] = numpy.sum(saturated[:,col])
        combined[:,4] = avg
        combined[:,5] = col
        
        this_sat = [col, max_y, min_y, numpy.sum(saturated[:,col]), avg]
        satcolumn_info.append(this_sat)

        numpy.savetxt(dump1, combined)
        print >>dump1, "\n\n\n"

    dump1.close()

    dump2 = open(output_txt+".sat2", "a")
    numpy.savetxt(dump2, numpy.array(satcolumn_info))

    if (len(sys.argv) <= 4):
        sys.exit(0)

    persist_data_all = []
    for file_id in range(4, len(sys.argv)):
        filename = sys.argv[file_id]
        hdulist = pyfits.open(filename)
        
        data = hdulist[extension].data
        overscan = numpy.median(data[:,500:])
        print "overscan-level=",overscan_level
        data -= overscan_level
        
        bglevel = numpy.median(data[10:480,10:490])
        data -= bglevel

        obsid = hdulist[0].header['OBSID']
        this_mjd = mjd = podi_persistency.get_mjd_from_timestamp(obsid[0:15])

        persist_data = []
        for i in range(len(satcolumn_info)):
            col, max_y, min_y, n_saturated, avg = satcolumn_info[i]

            strip = data[15:min_y-20,col]
            persistency = numpy.median(strip)

            this_pd = [this_mjd, (this_mjd - ref_mjd)*86400, persistency, col, max_y, min_y, n_saturated, avg]
            persist_data.append(this_pd)

        print >>dump2, "\n\n\n\n\n"
        #numpy.savetxt("%s___%d.persistency" % (output_txt, file_id), persist_data)
        numpy.savetxt(dump2, persist_data)

        persist_data_all.append(persist_data)

    dump2.close()

    # Now print the information into a file, sorted in blocks, 
    # one for each saturation-affected column
    dump3 = open(output_txt+".sat3", "w")
    dump4 = open(output_txt+".sat4", "w")
    timeseries = numpy.zeros(shape=(len(persist_data_all),2))
    for curcol in range(len(satcolumn_info)):
        

        for curtime in range(len(persist_data_all)):
            mjd, d_mjd, pers_level, col, max_y, min_y, n_saturated, avg = persist_data_all[curtime][curcol]
            timeseries[curtime,:] = [d_mjd, pers_level]

        # Now we also have the timeseries as numpy array, so let's fit a analytical function to it.
        # most generally, the decline should have this form:
        # p(t) = p0 + p1 * exp(-t/tau)
        # where ideally p0=0, but maybe not in reality to due background issues
        
        fitfunc = lambda p,x: p[0] * numpy.exp(-1 * x / p[1]) #+ p[2]
        errfunc = lambda p,x,y: y - fitfunc(p,x)

        pinit = [50,100,0]
        out = scipy.optimize.leastsq(errfunc, pinit, args=(timeseries[:,0], timeseries[:,1]), full_output=True)
        pfinal = out[0]
        
        fitted = numpy.zeros(shape=(timeseries.shape[0],3))
        fitted[:,0:2] = timeseries[:]
        fitted[:,2] = fitfunc(pfinal, timeseries[:,0])
        col, max_y, min_y, n_saturated, avg = satcolumn_info[curcol]
        print >>dump4, col, max_y, min_y, n_saturated, avg, pfinal[0], pfinal[1], pfinal[2]

        for curtime in range(len(persist_data_all)):
            mjd, d_mjd, pers_level, col, max_y, min_y, n_saturated, avg = persist_data_all[curtime][curcol]
            print >>dump3, col, max_y, min_y, n_saturated, avg, mjd, d_mjd, pers_level, fitted[curtime,2]
        print >>dump3, "\n\n\n\n\n\n"


    dump3.close()
    dump4.close()


    sys.exit(0)

    data[:] = 0
    data[mask] = 1

    hdulist_out = [pyfits.PrimaryHDU(data=data)]

    txtout = open(output_txt, "w")

    for file in sys.argv[5:]:
        print >>txtout, "#", file, extension

    print >>txtout, ref_mjd, 0, 65535, 65535

    # Now go through each of the files and get the median level of the formerly saturated pixels
    for file in sys.argv[5:]:
        print "Measuring file",file
        this_hdu = pyfits.open(file)
        this_data = this_hdu[extension].data
        this_mjd = mjd = podi_persistency.get_mjd_from_timestamp(file[0:15])
#this_hdu[0].header['MJD-OBS']

        saturated = this_data[mask]
        #print saturated[0]

        mean = numpy.mean(saturated)
        median = numpy.median(saturated)

        hdudump = pyfits.ImageHDU(data=saturated)
        hdudump.header.update("MJD", this_mjd)
        #hdudump.header.update("MEAN", mean)
        #hdudump.header.update("MEDIAN", median)
        
        hdulist_out.append(hdudump)

        print >>txtout, this_mjd, this_mjd-ref_mjd, mean, median

        del this_data
        this_hdu.close()
        del this_hdu

    hduout = pyfits.HDUList(hdulist_out)
    hduout.writeto(output_fits, clobber=True)
