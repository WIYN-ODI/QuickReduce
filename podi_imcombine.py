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


from podi_definitions import *

   

def parallel_compute(queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist):
    #queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist = worker_args

    buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len_filelist))
    result_buffer = shmem_as_ndarray(shmem_results).reshape((size_x, size_y))
    
    while (True):
        cmd_quit, line = queue.get()
        if (cmd_quit):
            queue.task_done()
            return
        else:
            result_buffer[line,:] = numpy.mean(buffer[line,:,:], axis=1)
            queue.task_done()



def imcombine(filelist, outputfile):
    queue = multiprocessing.JoinableQueue()

    # Read the input parameters
    # Note that file headers are copied from the first file
    reference_filename = filelist[0]
    ref_hdulist = pyfits.open(reference_filename)
    filter = ref_hdulist[1].header['FILTER']

    # Create the primary extension of the output file
    primhdu = pyfits.PrimaryHDU()
    out_hdulist = [primhdu]

    #
    # Now loop over all extensions and compute the mean
    #
    for cur_ext in range(1, len(ref_hdulist)):
        # Check what OTA we are dealing with
        ref_fppos = ref_hdulist[cur_ext].header['FPPOS']

        stdout_write("\rWorking on OTA %s (#% 2d/% 2d) ..." % (ref_fppos, cur_ext, len(ref_hdulist)-1))

        # Allocate enough shared memory to load a single OTA from all files. The ahred part is
        # important to make communication between the main and the salve processes possible.
        size_x, size_y = ref_hdulist[cur_ext].data.shape[0], ref_hdulist[cur_ext].data.shape[1]
        shmem_buffer = multiprocessing.RawArray(ctypes.c_float, size_x*size_y*len(filelist))
        shmem_results = multiprocessing.RawArray(ctypes.c_float, size_x*size_y)
        
        # Extratc the shared memory buffer as numpy array to make things easier
        buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len(filelist)))

        # Set the full buffer to NaN
        buffer[:,:,:] = numpy.NaN
        
        # Copy the reference data
        buffer[:,:,0] = ref_hdulist[cur_ext].data[:,:]
        del ref_hdulist[cur_ext].data

        # Now open all the other files, look for the right extension, and copy their image data to buffer
        for file_number in range(1, len(filelist)):
            filename = filelist[file_number]
            hdulist = pyfits.open(filename)
            for i in range(1, len(hdulist)):
                fppos = hdulist[i].header['FPPOS']
                if (fppos == ref_fppos):
                    buffer[:,:,file_number] = hdulist[i].data[:,:]
                    break
            hdulist.close()
            del hdulist

        #
        # Set up the parallel processing environment
        #
        #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
        processes = []
        for i in range(number_cpus):
            worker_args = (queue, shmem_buffer, shmem_results,
                           size_x, size_y, len(filelist))
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

    filelist = sys.argv[2:]

    imcombine(filelist, outputfile)
