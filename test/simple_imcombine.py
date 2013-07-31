#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy
import scipy.stats


import Queue
import threading
import multiprocessing
import ctypes
import bottleneck

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

   

def parallel_compute(queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist, operation):
    #queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist = worker_args

    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len_filelist))
    result_buffer = shmem_as_ndarray(shmem_results).reshape((size_x, size_y))

    while (True):
        cmd_quit, line = queue.get()
        if (cmd_quit):
            queue.task_done()
            return

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
            #print "nanmedian"
            result_buffer[line,:] = scipy.stats.nanmedian(buffer[line,:,:], axis=1)

        elif (operation == "nanmean.bn"):
            x = numpy.array(buffer[line,:,:], dtype=numpy.float32)
            result_buffer[line,:] = bottleneck.nanmean(x, axis=1)

        elif (operation == "nanmedian.bn"):
            #print "nanmedian"
            x = numpy.array(buffer[line,:,:], dtype=numpy.float32)
            result_buffer[line,:] = bottleneck.nanmedian(x, axis=1)
            #result_buffer[line,:] = scipy.stats.nanmedian(buffer[line,:,:], axis=1)

        elif (operation == "nansum.bn"):
            x = numpy.array(buffer[line,:,:], dtype=numpy.float32)
            result_buffer[line,:] = bottleneck.nansum(x, axis=1)

        else:
            result_buffer[line,:] = numpy.mean(buffer[line,:,:], axis=1)             
            

        queue.task_done()



def imcombine(input_filelist, outputfile, operation):
    queue = multiprocessing.JoinableQueue()

    # First loop over all filenames and make sure all files exist
    filelist = []
    for file in input_filelist:
        if (os.path.isfile(file)):
            filelist.append(file)

    if (len(filelist) <= 0):
        stdout_write("No existing files found in input list, hence nothing to do!\n")
        return
    elif (len(filelist) == 1):
        stdout_write("Only 1 file to combine, save the hassle and copy the file!\n")
        hdulist = pyfits.open(filelist[0])
        hdulist.writeto(outputfile)
        return
    
    # Read the input parameters
    # Note that file headers are copied from the first file
    reference_filename = filelist[0]
    ref_hdulist = pyfits.open(reference_filename)

    # Create the primary extension of the output file
    primhdu = pyfits.PrimaryHDU()

    # Copy all headers from the reference HDU
    cards = ref_hdulist[0].header.ascardlist()
    for c in cards:
        primhdu.header.update(c.key, c.value, c.comment)

    # Add PrimaryHDU to list of OTAs that go into the output file
    out_hdulist = [primhdu]

    #
    # Now loop over all extensions and compute the mean
    #
    for cur_ext in range(1, len(ref_hdulist)):
        # Check what OTA we are dealing with

        stdout_write("\rCombining frames for extension %2d / %2d ..." % (cur_ext+1, len(ref_hdulist)))

        # Allocate enough shared memory to load a single OTA from all files. The shared part is
        # important to make communication between the main and the slave processes possible.
        size_x, size_y = ref_hdulist[cur_ext].data.shape[0], ref_hdulist[cur_ext].data.shape[1]
        shmem_buffer = multiprocessing.RawArray(ctypes.c_float, size_x*size_y*len(filelist))
        shmem_results = multiprocessing.RawArray(ctypes.c_float, size_x*size_y)
        
        # Extract the shared memory buffer as numpy array to make things easier
        buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len(filelist)))

        # Set the full buffer to NaN
        buffer[:,:,:] = numpy.NaN
        
        # Copy the reference data
        buffer[:,:,0] = ref_hdulist[cur_ext].data[:,:]
        del ref_hdulist[cur_ext].data

        # Now open all the other files, look for the right extension, and copy their image data to buffer
        for file_number in range(1, len(filelist)):
            filename = filelist[file_number]
            print "Loading frame",filename
            hdulist = pyfits.open(filename)

            buffer[:,:,file_number] = hdulist[cur_ext].data[:,:]
            hdulist.close()
            del hdulist

        #
        # Set up the parallel processing environment
        #
        #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
        processes = []
        for i in range(number_cpus):
            worker_args = (queue, shmem_buffer, shmem_results,
                           size_x, size_y, len(filelist), operation)
            p = multiprocessing.Process(target=parallel_compute, args=worker_args)
            p.start()
            processes.append(p)

        # Now compute median/average/sum/etc
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

        # Create new ImageHDU
        hdu = pyfits.ImageHDU()

        # Insert the imcombine'd frame into the output HDU
        hdu.data = numpy.copy(shmem_as_ndarray(shmem_results).reshape((size_x, size_y)))

        # Copy all headers from the reference HDU
        cards = ref_hdulist[cur_ext].header.ascardlist()
        for c in cards:
            hdu.header.update(c.key, c.value, c.comment)

        # Append the new HDU to the list of result HDUs
        out_hdulist.append(hdu)

        del hdu
        del shmem_buffer
        del shmem_results

    stdout_write(" writing results to file %s ..." % (outputfile))
    out_hdu = pyfits.HDUList(out_hdulist)
    clobberfile(outputfile)
    out_hdu.writeto(outputfile, clobber=True)
    out_hdu.close()
    del out_hdu
    del out_hdulist    
    stdout_write(" done!\n")

if __name__ == "__main__":

    outputfile = sys.argv[1]

    filelist = sys.argv[3:]
    #get_clean_cmdline()[2:]

    operation = sys.argv[2] #cmdline_arg_set_or_default("-op", "mean")
    print "operation:", operation
    print filelist

    imcombine(filelist, outputfile, operation)
