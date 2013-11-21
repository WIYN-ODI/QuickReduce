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


import Queue
import threading
import multiprocessing
import ctypes

fix_cpu_count = False
number_cpus = 2
max_cpu_count = -1

try:
    number_cpus = multiprocessing.cpu_count()
    print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
    if (number_cpus > max_cpu_count and max_cpu_count > 1):
        number_cpus = max_cpu_count
        print "... but using only %d of them!" % (number_cpus)
except:
    pass

#number_cpus = 1

from podi_definitions import *
import bottleneck
verbose = cmdline_arg_isset("-verbose")

def parallel_compute(queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist, operation):
    #queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist = worker_args

    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len_filelist))
    result_buffer = shmem_as_ndarray(shmem_results).reshape((size_x, size_y))

    while (True):
        cmd_quit, line = queue.get()
        if (cmd_quit):
            queue.task_done()
            break

        if (operation == "median"):
            result_buffer[line,:] = numpy.median(buffer[line,:,:], axis=1)

        elif (operation == "medsigclip"):
            # Do not use (yet), is slow as hell 
            # (maskedarrays are pure python, not C as all the rest)

            #print buffer[line,:,:].shape
            _sigma_plus  = numpy.ones(shape=(buffer.shape[1],buffer.shape[2])) * 1e9
            _sigma_minus = numpy.ones(shape=(buffer.shape[1],buffer.shape[2])) * 1e9
            _median = numpy.median(buffer[line,:,:], axis=1)

            nrep = 3
            valid_pixels = numpy.ma.MaskedArray(buffer[line,:,:])

            for rep in range(nrep):

                _median_2d = _median.reshape(_median.shape[0],1).repeat(buffer.shape[2], axis=1)
                _min = _median_2d - 3 * _sigma_minus
                _max = _median_2d + 3 * _sigma_plus

                #valid_pixels = numpy.ma.masked_inside(buffer[line,:,:], _min, _max)
                valid = (buffer[line,:,:] > _min) & (buffer[line,:,:] < _max)

                valid_pixels = numpy.ma.array(buffer[line,:,:], mask=valid)
                #valid_pixels = numpy.ma.MaskedArray(buffer[line,:,:], valid)

                #print _min.shape, valid.shape, valid_pixels.shape

                #if (numpy.sum(valid, axis=1).any() <= 0):
                #    break

                #_median = numpy.median(buffer[line,:,:][valid], axis=1)
                _median = numpy.median(valid_pixels, axis=1)
                if (rep < nrep-1):
                    #_sigma_plus = scipy.stats.scoreatpercentile(buffer[line,:,:][valid], 84) - _median
                    #_sigma_minus = _median - scipy.stats.scoreatpercentile(buffer[line,:,:][valid], 16) 
                    _sigma_plus = scipy.stats.scoreatpercentile(valid_pixels, 84) - _median
                    _sigma_minus = _median - scipy.stats.scoreatpercentile(valid_pixels, 16) 

            result_buffer[line,:] = _median

        elif (operation == "medclip"):
            intermediate = numpy.sort(buffer[line,:,:], axis=1)
            result_buffer[line,:] = numpy.median(intermediate[:,1:-2], axis=1)

        elif (operation == "min"):
            result_buffer[line,:] = numpy.min(buffer[line,:,:], axis=1)

        elif (operation == "max"):
            result_buffer[line,:] = numpy.max(buffer[line,:,:], axis=1)

        elif (operation == "nanmean"):
            result_buffer[line,:] = scipy.stats.nanmean(buffer[line,:,:], axis=1)

        elif (operation == "nanmedian"):
            result_buffer[line,:] = scipy.stats.nanmedian(buffer[line,:,:], axis=1)

        elif (operation == "nanmedian.bn"):
            x = numpy.array(buffer[line,:,:], dtype=numpy.float32)
            result_buffer[line,:] = bottleneck.nanmedian(x, axis=1)
            x = None
            del x
        elif (operation == "nanmean.bn"):
            x = numpy.array(buffer[line,:,:], dtype=numpy.float32)
            result_buffer[line,:] = bottleneck.nanmean(x, axis=1)
            x = None
            del x
        else:
            result_buffer[line,:] = numpy.mean(buffer[line,:,:], axis=1)             
            

        queue.task_done()

    buffer = None
    shmem_buffer = None
    del shmem_buffer
    del buffer
    sys.exit(0)

    return



def imcombine_data(datas, operation="nanmean"):

    # Allocate enough shared memory to load a single OTA from all files. The shared part is
    # important to make communication between the main and the slave processes possible.
    size_x, size_y = datas[0].shape[0], datas[0].shape[1]
    total_pixels = size_x*size_y*len(datas)
    # print "total pixel count",total_pixels
    shmem_buffer = multiprocessing.RawArray(ctypes.c_float, total_pixels) #size_x*size_y*len(datas))

    # Extract the shared memory buffer as numpy array to make things easier
    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len(datas)))

    # Set the full buffer to NaN
    buffer[:,:,:] = numpy.NaN

    # Now open all the other files, look for the right extension, and copy their image data to buffer
    for data_id in range(len(datas)):
        # stdout_write("copying %d" % (data_id))
        buffer[:,:,data_id] = datas[data_id][:,:]

    sizes = (size_x, size_y, len(datas))
    combined = imcombine_sharedmem_data(shmem_buffer, operation, sizes)

    del shmem_buffer
    return combined

def imcombine_sharedmem_data(shmem_buffer, operation, sizes):

    size_x, size_y, n_frames = sizes
    shmem_results = multiprocessing.RawArray(ctypes.c_float, size_x*size_y)

    #
    # Set up the parallel processing environment
    #
    queue = multiprocessing.JoinableQueue()

    #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
    processes = []
    for i in range(number_cpus):
        worker_args = (queue, shmem_buffer, shmem_results,
                       size_x, size_y, n_frames, operation)
        p = multiprocessing.Process(target=parallel_compute, args=worker_args)
        p.start()
        processes.append(p)

    # Now compute median/average/sum/etc
    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, n_frames))
    for line in range(buffer.shape[0]):
        #print "Adding line",line,"to queue"
        queue.put((False,line))

    # Tell all workers to shut down when no more data is left to work on
    for i in range(number_cpus):
        queue.put((True,None))

    # Once all command are sent out to the workers, join them to speed things up
    try:
        queue.join()
    except KeyboardInterrupt:
        for p in processes:
            p.terminate()
        sys.exit(-1)

    for p in processes:
        p.terminate()

    results = numpy.copy(shmem_as_ndarray(shmem_results).reshape((size_x, size_y)))

    del shmem_results
    del queue
    del buffer

    return results



def imcombine_subprocess(extension, filelist, shape, operation, queue, verbose):

    #
    # Allocate enough shared momory to hold all frames
    #
    size_x, size_y = shape[0], shape[1]
    shmem_buffer = multiprocessing.RawArray(ctypes.c_float, size_x*size_y*len(filelist))

    # Extract the shared memory buffer as numpy array to make things easier
    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len(filelist)))

    # Set the full buffer to NaN
    buffer[:,:,:] = numpy.NaN

    # Now open all files, look for the right extension, and copy their image data to buffer
    for file_number in range(len(filelist)):
        filename = filelist[file_number]
        hdulist = pyfits.open(filename)
        for i in range(1, len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue
            fppos = hdulist[i].header['EXTNAME']
            if (fppos == extension):
                buffer[:,:,file_number] = hdulist[i].data[:,:]
                break
        hdulist[i].data = None
        hdulist.close()
        del hdulist
        if (verbose): stdout_write("\n   Added file %s ..." % (filename))

    # stdout_write("\n   Starting imcombine for real ...")
    combined = imcombine_sharedmem_data(shmem_buffer, operation=operation, sizes=(size_x, size_y, len(filelist)))

    # put the imcombine'd data into the queue to return them to the main process
    queue.put(combined)

    # and kill this process, returning all its memory
    sys.exit(0)



def imcombine(input_filelist, outputfile, operation, return_hdu=False):
    # First loop over all filenames and make sure all files exist
    filelist = []
    for file in input_filelist:
        if (os.path.isfile(file)):
            filelist.append(file)

    if (len(filelist) <= 0):
        stdout_write("No existing files found in input list, hence nothing to do!\n")
        return
    elif (len(filelist) == 1):
        # stdout_write("Only 1 file to combine, save the hassle and copy the file!\n")
        hdulist = pyfits.open(filelist[0])
        if (return_hdu):
            return hdulist
        hdulist.writeto(outputfile)
        return
    
    # Read the input parameters
    # Note that file headers are copied from the first file
    reference_filename = filelist[0]
    ref_hdulist = pyfits.open(reference_filename)

    # Create the primary extension of the output file
    # Copy all headers from the reference HDU
    primhdu = pyfits.PrimaryHDU(header=ref_hdulist[0].header)

    # Add PrimaryHDU to list of OTAs that go into the output file
    out_hdulist = [primhdu]

    #
    # Now loop over all extensions and compute the mean
    #
    for cur_ext in range(1, len(ref_hdulist)):

        data_blocks = []
        # Check what OTA we are dealing with
        if (not is_image_extension(ref_hdulist[cur_ext])):
            continue
        ref_fppos = ref_hdulist[cur_ext].header['EXTNAME']

        stdout_write("\rCombining frames for OTA %s (#% 2d/% 2d) ..." % (ref_fppos, cur_ext, len(ref_hdulist)-1))

        #
        # Add some weird-looking construct to move the memory allocation and actual 
        # imcombine into a separate process. This has to do with if/how/when python releases 
        # memory (or not), causing a massive short-term memory leak.
        # With the additional process, the memory is onwed by the other process and the memory
        # is freed once we destroy this helper process.
        #
        return_queue = multiprocessing.JoinableQueue()
        worker_args=(ref_fppos, filelist, ref_hdulist[cur_ext].data.shape, operation, return_queue, verbose)
        p = multiprocessing.Process(target=imcombine_subprocess, args=worker_args)
        p.start()
        combined = return_queue.get()
        p.terminate()
        del p

        if (verbose): stdout_write(" done, creating fits extension ...")

        # Create new ImageHDU, insert the imcombined's data and copy the 
        # header from the reference frame
        hdu = pyfits.ImageHDU(data=combined, header=ref_hdulist[cur_ext].header)

        # Append the new HDU to the list of result HDUs
        out_hdulist.append(hdu)

        if (verbose): stdout_write(" done\n")

    #
    # At this point, add a fits table listing all filenames that went into this combination
    #
    filenames_only = []
    for file in filelist:
        dirname, filename = os.path.split(file)
        filenames_only.append(filename)

    out_hdulist[0].header["NCOMBINE"] = (len(filenames_only), "number of combined files")

    columns = [
        pyfits.Column(name='filenumber', format="I4", array=range(1,len(filenames_only)+1)),
        pyfits.Column(name='filename', format="A100", array=filenames_only),
        ]
    coldefs = pyfits.ColDefs(columns)
    tablehdu = pyfits.new_table(coldefs, tbtype="BinTableHDU")
    tablehdu.header["EXTNAME"] = "FILELIST"
    out_hdulist.append(tablehdu)

    #
    # All work done now, prepare to return the data or write it to disk
    #
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

if __name__ == "__main__":

    outputfile = get_clean_cmdline()[1]

    filelist = get_clean_cmdline()[2:]

    operation = cmdline_arg_set_or_default("-op", "mean")

    imcombine(filelist, outputfile, operation)
