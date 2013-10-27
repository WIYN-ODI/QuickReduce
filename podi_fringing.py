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
import pyfits
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import scipy
import scipy.stats
import scipy.optimize
import bottleneck

import Queue
import threading
import multiprocessing
import ctypes

from podi_definitions import *
import podi_imcombine
import podi_fitskybackground
import time

avg_sky_countrates = {
    "odi_i": 3.8,
    "odi_z": 4.0,
    "823_v2": 0.1,
}

number_cpus = 4

def make_fringing_template(input_filelist, outputfile, return_hdu=False, skymode='local'):

    # First loop over all filenames and make sure all files exist
    hdu_filelist = []
    for file in input_filelist:
        if (os.path.isfile(file)):
            hdu_filelist.append(pyfits.open(file))

    if (len(hdu_filelist) <= 0):
        stdout_write("No existing files found in input list, hence nothing to do!\n")
        return
    
    # Read the input parameters
    # Note that file headers are copied from the first file

    # Create the primary extension of the output file
    ref_hdulist = hdu_filelist[0]
    primhdu = pyfits.PrimaryHDU(header=ref_hdulist[0].header)

    # Add PrimaryHDU to list of OTAs that go into the output file
    out_hdulist = [primhdu]

    filtername = primhdu.header['FILTER']
    if (outputfile == "auto"):
        outputfile = "fringes__%s.fits" % (filtername)

    print "Output file=",outputfile


    #
    # Now loop over all extensions and compute the mean
    #
    for cur_ext in range(1, len(ref_hdulist)):
        
        data_blocks = []
        # Check what OTA we are dealing with
        if (not is_image_extension(ref_hdulist[cur_ext].header)):
            continue
        extname = ref_hdulist[cur_ext].header['EXTNAME']

        if (filtername in otas_for_photometry):
            useful_otas = otas_for_photometry[filtername]
            ota_id = int(extname[3:5])
            if (not ota_id in useful_otas):
                continue
        
        stdout_write("\rCombining frames for OTA %s (#% 2d/% 2d) ..." % (extname, cur_ext, len(ref_hdulist)-1))

        # Now open all the other files, look for the right extension, and copy their image data to buffer
        for file_number in range(0, len(filelist)):

            try:
                this_hdu = hdu_filelist[file_number][extname]
            except:
                continue

            # Skip all OTAs that are marked as video/guide OTAs
            cellmode = this_hdu.header['CELLMODE']
            if (cellmode.find("V") >= 0):
                continue
            
            skylevel = this_hdu.header['SKY_MEDI']
            if (skymode == 'global'):
                skylevel = hdu_filelist[file_number][0].header['SKYLEVEL']
            if ("EXPTIME" in hdu_filelist[file_number][0].header):
                exptime = hdu_filelist[file_number][0].header['EXPTIME']
                filter = hdu_filelist[file_number][0].header['FILTER']
                if (filter in avg_sky_countrates):
                    max_skylevel = 2 * avg_sky_countrates[filter] * exptime
                    if (skylevel > max_skylevel):
                        stdout_write(" (%.1f)" % (skylevel/exptime))
                        continue

            fringing = (this_hdu.data - skylevel) / skylevel
            stdout_write(" %.1f" % (skylevel/exptime))
            data_blocks.append(fringing)

            # delete the data block to free some memory, since we won't need it anymore
            del this_hdu.data

        stdout_write(" combining ...")
        #combined = podi_imcombine.imcombine_data(data_blocks, "nanmedian")
        combined = podi_imcombine.imcombine_data(data_blocks, "nanmedian.bn")

        # Create new ImageHDU
        # Insert the imcombine'd frame into the output HDU
        # Copy all headers from the reference HDU
        # stdout_write(" creating HDU ...")
        hdu = pyfits.ImageHDU(header=ref_hdulist[cur_ext].header, data=combined)

        # Append the new HDU to the list of result HDUs
        out_hdulist.append(hdu)
        stdout_write(" done!\n")

        del hdu

    # Add the assumed skylevel countrate to primary header so we can use it
    # when it comes to correcting actual data
    out_hdulist[0].header["SKYCNTRT"] = (avg_sky_countrates[filter], "average countrate of dark sky")

    return_hdu = False
    out_hdu = pyfits.HDUList(out_hdulist)
    if (not return_hdu and outputfile != None):
        stdout_write(" writing results to file %s ..." % (outputfile))
        clobberfile(outputfile)
        out_hdu.writeto(outputfile, clobber=True)
        out_hdu.close()
        del out_hdu
        del out_hdulist
        stdout_write(" done!\n")
    elif (return_hdu):
        stdout_write(" returning HDU for further processing ...\n")
        return out_hdu
    else:
        stdout_write(" couldn't write output file, no filename given!\n")

    return


def compute_fringe_scale(datahdu, fringehdu):

    extname = datahdu.header['EXTNAME']
    print "\n\n\n",extname,"\n"
    data = datahdu.data
    fringe = fringehdu.data

    # rebin data to cut down on processing time
    binning = 8
    data_binned = rebin_image(data, binning)
    fringe_binned = rebin_image(fringe, binning)

    # compute mean fringe amplitude
    mean_fringe = bottleneck.nanmean(fringe_binned)
    fringe_binned -= mean_fringe
    print "mean fringe =",mean_fringe

    skylevel = datahdu.header['SKY_MEDI']
    skysub = data_binned - skylevel

    def min_stat(scale, data, fringe):
        return bottleneck.nanmean( numpy.fabs((skysub - scale*fringe)*fringe) )
        #return bottleneck.nanmean( ((skysub - scale*fringe)*fringe)**2 )

    initial_guess = 100
    res = scipy.optimize.fmin(min_stat, initial_guess, 
                              args=(skysub, fringe_binned),
                              full_output=True)
    print res
    res_fits = [[res[0][0], res[1], res[2], res[3], res[4]]]

    return res_fits

def mpworker_fringe_scale(queue_jobs, queue_return): #datahdu, fringehdu):

    while (True):
        cmd_quit, data = queue_jobs.get()
        if (cmd_quit):
            break

        # Read all the data for the job to do
        # datahdu, fringehdu = data
        data_filename, fringe_filename, extname = data

        # Open both fits files and select the right extension
        data_hdulist = pyfits.open(data_filename)
        fringe_hdulist = pyfits.open(fringe_filename)
        datahdu = data_hdulist[extname]
        fringehdu = fringe_hdulist[extname]

        # Do the calculation
        res_fits = compute_fringe_scale(datahdu, fringehdu)

        # close the files and hopefully free up some memory
        data_hdulist.close()
        fringe_hdulist.close()

        # return the data to the main process
        return_data = (res_fits)
        queue_return.put(return_data)
        queue_jobs.task_done()

    return


def match_subtract_fringing(data_filename, fringe_filename, verbose=True, output=None):

#    if (not type(fringe_hdulist) is pyfits.HDUList):
#        # The fringe variable is a filename
#        fringe_filename = fringe_hdulist
    fringe_hdulist = pyfits.open(fringe_filename)
    data_hdulist = pyfits.open(data_filename)
        
    if (verbose): stdout_write("Creating queues\n")
    queue = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.Queue()

    if (verbose): stdout_write("Handing out work\n")
    all_results = None
    number_extensions = 0
    for ext in range(1, len(data_hdulist)):
        if (type(data_hdulist[ext]) != pyfits.hdu.image.ImageHDU):
            continue

        # Load data for each extension/OTA
        extname = data_hdulist[ext].header['EXTNAME']
        print "Queuing work on",extname

        cmd_data = (data_filename, fringe_filename, extname)
        queue.put( (False, cmd_data) )

        #queue.put( (False, (data_hdulist[extname], fringe_hdulist[extname]) ) )
        number_extensions += 1

    # Start all worker processes
    if (verbose): stdout_write("Starting workers\n")
    processes = []
    for i in range(number_cpus):
        worker_args = (queue, return_queue)
        p = multiprocessing.Process(target=mpworker_fringe_scale, args=worker_args)
        p.start()
        processes.append(p)
        time.sleep(0.01)

    # Send one termination command to each worker
    if (verbose): stdout_write("seinding quit commands\n")
    for p in processes:
        queue.put( (True, None) )

    # Collect all results from all the workers
    if (verbose): stdout_write("Collecting work\n")
    all_results = None
    for i in range(number_extensions):
        res_fits = return_queue.get()
        if (res_fits != None):
            all_results = res_fits if (all_results == None) \
                else numpy.append(all_results, numpy.array(res_fits), axis=0)

    if (verbose): stdout_write("computing scaling\n")
    print all_results[:, 0:3]

    quality_sorted = numpy.sort(all_results[:,1])
    use_for_median = all_results[:,1] < quality_sorted[5]
        
    median_scaling = numpy.median(all_results[:,0][use_for_median])
    std_scaling = numpy.std(all_results[:,0][use_for_median])

    if (verbose): stdout_write("doing reduction\n")

    for ext in range(1, len(data_hdulist)):
        if (type(data_hdulist[ext]) != pyfits.hdu.image.ImageHDU):
            continue

        extname = data_hdulist[ext].header['EXTNAME']
        data_hdulist[ext].data -= (fringe_hdulist[extname].data * median_scaling)

        data_hdulist[ext].header["FRNG_SCL"] = (median_scaling, "fringe scaling")
        data_hdulist[ext].header["FRNG_STD"] = (std_scaling, "fringe scaling std.dev.")

    if (output != None):
        data_hdulist.writeto(output, clobber=True)
        data_hdulist.close()
        return median_scaling, std_scaling, None

    print data_hdulist
    return median_scaling, std_scaling, data_hdulist




def get_fringe_scaling(data, fringe, region_file):
    """
    This routine implements the technique for determining the optimal
    fringe scaling outlined in Snodgrass & Carry 2013, ESO Messenger 152, 14.

    In short, it determines the mean value in a number of regions selected
    visually to represent dark- and bright spots in the fringe map. The 
    difference between bright and dark represents the fringe amplitude. The 
    same measurements are taken for the same regions in the data frame,
    informing about the fringe amplitude in the data frame. The ratio between
    the two amplitudes represents the required fringe scaling factor.

    Input variables are:
    - the data frame as 2-d numpy array
    - the fringe map as 2-d numpy array
    - the filename of a ds9 region file defining the regions

    Output:
    A vector of measurements, column 6 of which is the scaling factor for
    each region.
    """
    
    if (not os.path.isfile(region_file)):
        return None

    # Read the region file
    regfile = open(region_file, "r")
    entries = regfile.readlines()
    regfile.close()

    data = numpy.array(data, dtype=numpy.float32)
    fringe = numpy.array(fringe, dtype=numpy.float32)

    all_vectors = []
    for line in entries:
        #print line
        if (line.find("vector(") < 0):
            continue
        start = line.find("vector(") + 7
        end = line.find(")", start)
        #print line[start:end]," --> ",
        items = line[start:end].split(",")
        #print items

        origx = int(float(items[0]))
        origy = int(float(items[1]))
        length = float(items[2])
        angle = float(items[3]) 

        dx = length * math.cos(math.radians(angle))
        dy = length * math.sin(math.radians(angle))

        darkx, darky = origx, origy
        lightx, lighty = int(origx+dx), int(origy+dy)

        #print darkx, darky, lightx, lighty

        boxwidth = 10
        data_dark = bottleneck.nanmedian(data[darky-boxwidth:darky+boxwidth,darkx-boxwidth:darkx+boxwidth])
        data_light = bottleneck.nanmedian(data[lighty-boxwidth:lighty+boxwidth,lightx-boxwidth:lightx+boxwidth])

        fringe_dark = bottleneck.nanmedian(fringe[darky-boxwidth:darky+boxwidth,darkx-boxwidth:darkx+boxwidth])
        fringe_light = bottleneck.nanmedian(fringe[lighty-boxwidth:lighty+boxwidth,lightx-boxwidth:lightx+boxwidth])

        data_diff = data_light - data_dark
        fringe_diff = fringe_light - fringe_dark
        scaling = data_diff / fringe_diff
        #print data_light, data_dark, fringe_light, fringe_dark, data_diff, fringe_diff, scaling

        one_vector = [data_light, data_dark, fringe_light, fringe_dark, data_diff, fringe_diff, scaling]
        all_vectors.append(one_vector)

    vecs = numpy.array(all_vectors)

    return vecs


if __name__ == "__main__":

    if (cmdline_arg_isset("-singles")):
        for filename in get_clean_cmdline()[1:]:
            outputfile = filename[:-5]+".fringe.fits"
            if (cmdline_arg_isset("-noclobber") and os.path.isfile(outputfile)):
                stdout_write("\n%s already exists, skipping!\n\n" % (outputfile))
                continue

            stdout_write("Converting %s to fringemask ...\n" % (filename))
            hdulist = pyfits.open(filename)

            out_hdu = [pyfits.PrimaryHDU(header=hdulist[0].header)]

            for ext in range(len(hdulist)):
                if (not is_image_extension(hdulist[ext].header)):
                    continue
                
                # Skip all OTAs that are marked as video/guide OTAs
                cellmode = hdulist[ext].header['CELLMODE']
                if (cellmode.find("V") >= 0):
                    continue

                skylevel = hdulist[ext].header['SKY_MEDI']
                
                fringing = (hdulist[ext].data - skylevel) / skylevel
                stdout_write("   %s = %.1f\n" % (hdulist[ext].header['EXTNAME'], skylevel))

                out_hdu.append(pyfits.ImageHDU(header=hdulist[ext].header,
                                               data=fringing))

            stdout_write("writing (%s)..." % (outputfile))
            out_hdulist = pyfits.HDUList(out_hdu)
            out_hdulist.writeto(outputfile, clobber=True)
            stdout_write(" done!\n\n")

        sys.exit(0)
                               
    
    if (cmdline_arg_isset("-make_template")):
        outputfile = get_clean_cmdline()[1]
        filelist = get_clean_cmdline()[2:]
        operation = cmdline_arg_set_or_default("-op", "mean")
        operation = cmdline_arg_set_or_default("-op", "nanmedian.bn")
        make_fringing_template(filelist, outputfile, operation)

        sys.exit(0)

    if (cmdline_arg_isset("-samplesky")):
        import podi_fitskybackground
        boxsize = int(cmdline_arg_set_or_default("-boxsize", "15"))
        count = int(cmdline_arg_set_or_default("-count", "12500"))

        print "Using %d boxes with +/- %d pixel width" % (count, boxsize)

        filename = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]

        print filename,"-->",outputfile
        hdulist = pyfits.open(filename)
        for i in range(3,len(hdulist)):
            if (True): #is_image_extension(hdulist[i].header)):
                extname = hdulist[i].header['EXTNAME']
                stdout_write("\rWorking on %s ..." % (extname))
                regions = numpy.array(podi_fitskybackground.sample_background(hdulist[i].data, None, None, min_found=count, boxwidth=boxsize))
                
                median = numpy.median(regions[:,4])
                std = numpy.std(regions[:,4])


                histogram, binedges = numpy.histogram(regions[:,4], bins=500, range=(median-20*std,median+20*std))
                nice_histogram = numpy.empty(shape=(histogram.shape[0], 3))
                nice_histogram[:,0] = binedges[:-1]
                nice_histogram[:,1] = binedges[1:]
                nice_histogram[:,2] = histogram[:]
                
                dumpfile = "%s.%s" % (outputfile, extname)
                df = open(dumpfile, "w")
                #numpy.savetxt(df, regions)
                #print >>df, "\n\n\n\n\n\n\n"
                #numpy.savetxt(df, nice_histogram)
                #print >>df, "\n\n\n\n\n\n\n"

                validpixels = hdulist[i].data
                histogram, binedges = numpy.histogram(validpixels, bins=1000, range=(median-10*std,median+10*std))
                nice_histogram = numpy.empty(shape=(histogram.shape[0], 3))
                nice_histogram[:,0] = binedges[:-1]
                nice_histogram[:,1] = binedges[1:]
                nice_histogram[:,2] = histogram[:]

                numpy.savetxt(df, nice_histogram)

                df.close()
            

    if (cmdline_arg_isset("-correct")):

        inputframe = get_clean_cmdline()[1]
        template = get_clean_cmdline()[2]
        scaling = int(get_clean_cmdline()[3])
        output = get_clean_cmdline()[4]

        inputhdu = pyfits.open(inputframe)
        templatehdu = pyfits.open(template)

        for i in range(len(inputhdu)):
            if (not is_image_extension(inputhdu[i].header)):
                continue

            extname = inputhdu[i].header['EXTNAME']
            stdout_write("\rWorking on %s ... " % (extname))

            try:
                fringemap = templatehdu[extname].data * scaling
                inputhdu[i].data -= fringemap
            except:
                print "Problem with extension",extname
                pass
                continue

        stdout_write(" writing ...")

        inputhdu.writeto(output, clobber=True)
        
        stdout_write(" done!\n")

    if (cmdline_arg_isset("-optimize")):
        dataframe = get_clean_cmdline()[1]
        fringemap = get_clean_cmdline()[2]
        
        datahdu = pyfits.open(dataframe)
        fringehdu = pyfits.open(fringemap)

        for ext in range(1,14): #len(datahdu)):

            # Load data for each extension/OTA
            extname = datahdu[ext].header['EXTNAME']
            print "\n\n\n",extname,"\n"
            data = datahdu[extname].data
            fringe = fringehdu[extname].data

            # rebin data to cut down on processing time
            binning = 16
            data_binned = rebin_image(data, binning)
            fringe_binned = rebin_image(fringe, binning)
            
            output_hdulist = [pyfits.PrimaryHDU(data = fringe_binned)]
            output_hdulist[0].header["EXTNAME"] = "MAP"

            valid_both = numpy.isfinite(data_binned) & numpy.isfinite(fringe_binned)
            data_binned = data_binned[valid_both]
            fringe_binned = fringe_binned[valid_both]

            skylevel = datahdu[ext].header['SKY_MEDI']

            data_binned -= skylevel

            percent_limit = 15

            # Now select top 10 and bottom 10% of all intensities in the fringe map
            top10 = scipy.stats.scoreatpercentile(fringe_binned.ravel(), 100-percent_limit)
            bottom10 = scipy.stats.scoreatpercentile(fringe_binned.ravel(), percent_limit)
            print top10, bottom10
            #sys.exit(0)

            select_top10 = fringe_binned >= top10
            select_bottom10 = fringe_binned <= bottom10
            #fringe_top10 = numpy.median(fringe_binned[select_top10])
            #fringe_bottom10 = numpy.median(fringe_binned[select_bottom10])

#            x = fringe.copy()
#            output_hdulist.append(pyfits.ImageHDU(data=x))
#            x = fringe.copy()
#            x[x < top10] = numpy.NaN
#            output_hdulist.append(pyfits.ImageHDU(data=x))
#            x = fringe.copy()
#            x[x > bottom10] = numpy.NaN
#            output_hdulist.append(pyfits.ImageHDU(data=x))

#            x = data.copy()
#            output_hdulist.append(pyfits.ImageHDU(data=x))
#            x = data.copy()
#            x[fringe < top10] = numpy.NaN
#            output_hdulist.append(pyfits.ImageHDU(data=x))
#            x = data.copy()
#            x[fringe > bottom10] = numpy.NaN
#            output_hdulist.append(pyfits.ImageHDU(data=x))

            
            def clip3sigma(data):
                med, sigma = 0, 1e9
                valid = numpy.isfinite(data)
                for rep in range(3):
                    med = numpy.median(data[valid])
                    sigma = 0.5 * (scipy.stats.scoreatpercentile(data[valid], 84) - 
                                   scipy.stats.scoreatpercentile(data[valid], 16))
                    valid = (data < med + 3*sigma) & (data > med - 3*sigma)
                return numpy.median(data[valid])

            fringe_top10 = clip3sigma(fringe_binned[select_top10])
            fringe_bottom10 = clip3sigma(fringe_binned[select_bottom10])

            # Now get intensities in the original frame for the same pixels
            #sky_top10 = numpy.median(data_binned[select_top10])
            #sky_bottom10 = numpy.median(data_binned[select_bottom10])
            sky_top10 = clip3sigma(data_binned[select_top10])
            sky_bottom10 = clip3sigma(data_binned[select_bottom10])

            #sys.exit(0)

            print top10, bottom10
            print fringe_top10, fringe_bottom10, fringe_top10-fringe_bottom10
            print sky_top10, sky_bottom10, sky_top10-sky_bottom10
            fringe_scaling = (sky_top10-sky_bottom10)/(fringe_top10-fringe_bottom10)
            print "scaling = ", fringe_scaling

            def save_histogram(target, data, nbins=100):
                valid = numpy.isfinite(data)
                med = numpy.median(data[valid])
                ls, us = scipy.stats.scoreatpercentile(data[valid], 16), scipy.stats.scoreatpercentile(data[valid], 84)
                sigma = 0.5 * (us - ls)
                min = med - 3*sigma
                max = med + 3*sigma
                histogram, binedges = numpy.histogram(data[valid], bins=nbins, range=(min,max))
                nice_histogram = numpy.empty(shape=(histogram.shape[0], 3))
                nice_histogram[:,0] = binedges[:-1]
                nice_histogram[:,1] = binedges[1:]
                nice_histogram[:,2] = histogram[:]
                numpy.savetxt(target, nice_histogram)

            dump = open("xxx.dump"+extname, "w")
            save_histogram(dump, fringe_binned[select_top10])
            print >>dump, "\n\n\n\n\n"
            save_histogram(dump, fringe_binned[select_bottom10])
            print >>dump, "\n\n\n\n\n"
            save_histogram(dump, data_binned[select_top10])
            print >>dump, "\n\n\n\n\n"
            save_histogram(dump, data_binned[select_bottom10])
            dump.close()

#            output_hdulist.append(pyfits.ImageHDU(data=data))
#            #output_hdulist.append(pyfits.ImageHDU(data=(data - fringe_scaling * fringe)))
#            output_hdulist.append(pyfits.ImageHDU(data=(data - 0.5*fringe_scaling * fringe)))

#            outhdulist = pyfits.HDUList(output_hdulist)
#            outhdulist.writeto("output.fits", clobber=True)

            datahdu[ext].data -= skylevel #fringe_scaling * fringe

        datahdu.writeto("corrected.fits", clobber=True)

        sys.exit(0)



    if (cmdline_arg_isset("-rmfringe")):
        dataframe = get_clean_cmdline()[1]
        fringemap = get_clean_cmdline()[2]
        
        datahdu = pyfits.open(dataframe)
        fringehdu = pyfits.open(fringemap)

        all_results = None
        for ext in range(1,14): #len(datahdu)):

            # Load data for each extension/OTA
            extname = datahdu[ext].header['EXTNAME']
            print "\n\n\n",extname,"\n"
            data = datahdu[extname].data
            fringe = fringehdu[extname].data

            # rebin data to cut down on processing time
            binning = 8
            data_binned = rebin_image(data, binning)
            fringe_binned = rebin_image(fringe, binning)

            # compute mean fringe amplitude
            mean_fringe = bottleneck.nanmean(fringe_binned)
            fringe_binned -= mean_fringe
            print "mean fringe =",mean_fringe

            skylevel = datahdu[ext].header['SKY_MEDI']
            skysub = data_binned - skylevel

            def min_stat(scale, data, fringe):
                return bottleneck.nanmean( numpy.fabs((skysub - scale*fringe)*fringe) )
                #return bottleneck.nanmean( ((skysub - scale*fringe)*fringe)**2 )

            initial_guess = 100
            res = scipy.optimize.fmin(min_stat, initial_guess, 
                                      args=(skysub, fringe_binned),
                                      full_output=True)
            print res
            res_fits = [[res[0][0], res[1], res[2], res[3], res[4]]]
            all_results = numpy.array(res_fits) if all_results == None else numpy.append(all_results, res_fits, axis=0)

#            xopt, fopt, iter, funcalls, warnflag, allvecs = res
#            print xopt
 
#            for s in range(0,2000,25):
#                print extname, s, min_stat(s, skysub, fringe_binned)

            #print "success=",res.success
            #print "res.X=",res.x
            #print "msg=",res.message
            #print "\n\n"

            

        print all_results

        # Now that we got all best-fit factors, keep only the 
        # best 5 fits
        quality_sorted = numpy.sort(all_results[:,1])
        use_for_median = all_results[:,1] < quality_sorted[5]
        
        median_scaling = numpy.median(all_results[:,0][use_for_median])
        std_scaling = numpy.std(all_results[:,0][use_for_median])

        print median_scaling, std_scaling


        #datahdu.writeto("corrected.fits", clobber=True)

        sys.exit(0)




    if (cmdline_arg_isset("-matchsubtract")):
        dataframe = get_clean_cmdline()[1]
        fringemap = get_clean_cmdline()[2]
        
        #datahdu = pyfits.open(dataframe)
        #fringehdu = pyfits.open(fringemap)
        #match_subtract_fringing(datahdu, fringehdu)

        datahdu = match_subtract_fringing(dataframe, fringemap, output="matchsubtract.fits")
        
        #datahdu.writeto("matchsubtract.out.fits", clobber=True)


    if (cmdline_arg_isset("-esomethod")):
        fringe_frame = get_clean_cmdline()[1]
        data_frame = get_clean_cmdline()[2]
        output_filename = get_clean_cmdline()[3]
        
        data_hdulist = pyfits.open(data_frame)
        fringe_hdulist = pyfits.open(fringe_frame)

        filter_name = data_hdulist[0].header['FILTER']
        print "\nThis is filter",filter_name,"\n"

        all_vecs = None
        for ext in range(1, len(data_hdulist)):
            if (type(data_hdulist[ext]) != pyfits.hdu.image.ImageHDU):
                continue
            if (not data_hdulist[ext].header['CELLMODE'].find("V") < 0):
                # This is marked as a guide CCD and most likely useless
                continue

            extname = data_hdulist[ext].header['EXTNAME']

            data = data_hdulist[extname].data
            fringe = fringe_hdulist[extname].data

            region_file = "fringe__%s__%s.reg" % (filter_name, extname[0:5])
            vecs = get_fringe_scaling(data, fringe, region_file) 

            if (not vecs == None):
                all_vecs = vecs if all_vecs == None else numpy.append(all_vecs, vecs, axis=0)

        scaling_factors = all_vecs[:,6]
        valid = scaling_factors > 0
                
        median_scaling = numpy.median(scaling_factors[valid])
        print "median scaling=",median_scaling

        for rep in range(3):
            lsig = scipy.stats.scoreatpercentile(scaling_factors[valid], 16)
            hsig = scipy.stats.scoreatpercentile(scaling_factors[valid], 84)
            median = numpy.median(scaling_factors[valid])
            sigma = 0.5 * (hsig - lsig)
            #print median, sigma
            valid = (scaling_factors > median-3*sigma) & (scaling_factors < median+3*sigma)

        #all_vecs[:,6][valid == False] *= -1.
        #numpy.savetxt("all_vecs.dat", all_vecs)

        final_scaling = numpy.median(scaling_factors[valid])
        print "final scaling (OTA %s): %.3f" % (extname, final_scaling)

        # Now do the correction
        for ext in range(1, len(data_hdulist)):
            if (type(data_hdulist[ext]) != pyfits.hdu.image.ImageHDU):
                continue
            data_hdulist[extname].data -= (fringe_hdulist[extname].data * final_scaling)

        data_hdulist.writeto(output_filename, clobber=True)
