#! /usr/bin/env python3
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

This module contains all functionality for combining a number of input frames
into one output frame by performing a user-sdefined function on the range of
input values of each pixel.

To improve performance, the actual image-combination step is executed in
parallel, split by line. Furthermore, each extension is computed in a new
process, enabling better memory management by the operating system. This enables
memory-limited operation with very little memory loss/leakage between
extensions.

Example:
--------
    For unbinned pODI (13 OTAs) frames, memory demand is roughly 
    1 GB + 65MB x numer of frames. Assuming 6 GB of available memory, this would
    practically (avoiding memory swapping) to 75 input frames.


Standalone option:
------------------

    Combine a set of input images and write output to disk

    ``podi_imcombine.py (-op=xxx) output.fits file1.fits file2.fits ...``


"""


import sys
import os
import astropy.io.fits as pyfits
import numpy
#numpy.seterr(divide='ignore', invalid='ignore')
import scipy
import scipy.stats

import podi_logging
import logging
import time

from podi_definitions import *
from podi_commandline import *
import podi_associations

import queue
import threading
import multiprocessing
import ctypes
import collections

from sharedmemory import SharedMemory

fix_cpu_count = False
number_cpus = 2
max_cpu_count = -1

try:
    number_cpus = multiprocessing.cpu_count()
    # print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
    if (number_cpus > max_cpu_count and max_cpu_count > 1):
        number_cpus = max_cpu_count
        # print "... but using only %d of them!" % (number_cpus)
except:
    pass

#number_cpus = 1

from podi_definitions import *
from podi_commandline import *
import bottleneck
verbose = cmdline_arg_isset("-verbose")

import podi_cython

def weighted_mean(_line):
    max_weight = 50
    
    # print _line.shape
    
    median_2d = bottleneck.nanmedian(_line, axis=1).reshape(_line.shape[0],1).repeat(_line.shape[1], axis=1)
    std = bottleneck.nanstd(_line, axis=1)
    std_2d = std.reshape(_line.shape[0],1).repeat(_line.shape[1], axis=1)
    
    weight_2d = numpy.fabs(std_2d / (_line - median_2d))
#    weight_2d[weight_2d > max_weight] = max_weight
    weight_2d[numpy.isinf(weight_2d)] = max_weight
    
    for i in range(3):
        avg = bottleneck.nansum(_line*weight_2d, axis=1)/bottleneck.nansum(weight_2d, axis=1)
        avg_2d = avg.reshape(_line.shape[0],1).repeat(_line.shape[1], axis=1)
        
        std = numpy.sqrt(bottleneck.nansum(((_line - avg_2d)**2 * weight_2d), axis=1)/bottleneck.nansum(weight_2d, axis=1))
        std_2d = std.reshape(_line.shape[0],1).repeat(_line.shape[1], axis=1)
        
        weight_2d = numpy.fabs(std_2d / (_line - avg_2d))
        #weight_2d[weight_2d > max_weight] = max_weight
        weight_2d[numpy.isinf(weight_2d)] = max_weight
    
    return bottleneck.nansum(_line*weight_2d, axis=1)/bottleneck.nansum(weight_2d, axis=1)

def parallel_compute(queue, return_queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist, operation):
    #queue, shmem_buffer, shmem_results, size_x, size_y, len_filelist = worker_args

    # buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len_filelist))
    buffer = shmem_buffer.to_ndarray()
    # result_buffer = shmem_as_ndarray(shmem_results).reshape((size_x, size_y))
    result_buffer = shmem_results.to_ndarray()

    logger = logging.getLogger("ParallelImcombine")
    logger.debug("Operation: %s, #samples/pixel: %d" % (operation, len_filelist))


    while (True):
        line = queue.get()
        if (line is None):
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

        elif (operation == "sigclipx"):
            stdout_write(".")
            rep_count = 2

            _line = buffer[line,:,:].astype(numpy.float32)
            # print _line.shape

            mask = numpy.isfinite(_line)
            #print "line.shape=",_line.shape
            # numpy.savetxt("line_block_%d.dat" % (line), _line)

            def sigclip_pixel(pixelvalue):
                mask = numpy.isfinite(pixelvalue)
                old_mask = mask
                rep = 0
                while (rep < rep_count and numpy.sum(mask) > 3):
                    old_mask = mask

                    mss = scipy.stats.scoreatpercentile(pixelvalue[mask], [16,50,84])
                    
                    lower = mss[1] - 3 * (mss[1] - mss[0]) # median - 3*sigma
                    upper = mss[1] + 3 * (mss[2] - mss[1]) # median + 3*sigma

                    mask = (pixelvalue > lower) & (pixelvalue < upper)

                    rep += 1
                    if (rep == rep_count or numpy.sum(mask) < 3):
                        mask = old_mask

                return numpy.mean(pixelvalue[mask])

            result_buffer[line,:] = [sigclip_pixel(_line[x,:]) for x in range(_line.shape[0])]


        elif (operation == "sigmaclipmean"):
            _line = buffer[line,:,:].astype(numpy.float64)
            output = numpy.zeros(shape=(_line.shape[0]))
            podi_cython.sigma_clip_mean(_line, output)
            result_buffer[line,:] = output

        elif (operation == "sigmaclipmedian"):
            _line = buffer[line,:,:].astype(numpy.float64)
            output = numpy.zeros(shape=(_line.shape[0]))
            podi_cython.sigma_clip_median(_line, output)
            result_buffer[line,:] = output

        elif (operation == "weightedmean"):
            _line = buffer[line,:,:].astype(numpy.float32)
            result_buffer[line,:] = weighted_mean(_line)

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
            
        return_queue.put(line)
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
    shmem_buffer = SharedMemory(ctypes.c_float, (size_x, size_y, len(datas)))
    # multiprocessing.RawArray(ctypes.c_float, total_pixels) #size_x*size_y*len(datas))

    # Extract the shared memory buffer as numpy array to make things easier
    buffer = shmem_buffer.to_ndarray()
    # shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, len(datas)))

    # Set the full buffer to NaN
    buffer[:,:,:] = numpy.nan

    # Now open all the other files, look for the right extension, and copy their image data to buffer
    for data_id in range(len(datas)):
        # stdout_write("copying %d" % (data_id))
        buffer[:,:,data_id] = datas[data_id][:,:]

    sizes = (size_x, size_y, len(datas))
    combined = imcombine_sharedmem_data(shmem_buffer, operation, sizes)
    shmem_buffer.free()

    del shmem_buffer
    return combined

def imcombine_sharedmem_data(shmem_buffer, operation, sizes):

    size_x, size_y, n_frames = sizes
    shmem_results = SharedMemory(ctypes.c_float, (size_x,size_y))
    # multiprocessing.RawArray(ctypes.c_float, size_x*size_y)

    logger = logging.getLogger("CombineMgr")

    #
    # Set up the parallel processing environment
    #
    queue = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.Queue()

    # Now compute median/average/sum/etc
    # buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, n_frames))
    buffer = shmem_buffer.to_ndarray()
    for line in range(buffer.shape[0]):
        #print "Adding line",line,"to queue"
        queue.put(line)
        
    lines_done = numpy.zeros((buffer.shape[0]), dtype=numpy.bool)
    lines_read = 0

    #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
    processes = []
    for i in range(number_cpus):
        worker_args = (queue, return_queue,
                       shmem_buffer, shmem_results,
                       size_x, size_y, n_frames, operation)
        p = multiprocessing.Process(target=parallel_compute, args=worker_args)
        p.start()
        processes.append(p)

    while (lines_read < buffer.shape[0] and numpy.sum(lines_done) < buffer.shape[0]):
        try:
            line = return_queue.get(timeout=5)
            lines_read += 1
            try:
                lines_done[line] = True
            except:
                pass
        except Queue.Empty:
            logger.error("Encountered timeout while combinging data")
            # something bad has happened to one of the workers
            # find one of the lines that has not been processed yet
            missing_lines = (numpy.arange(buffer.shape[0]))[~lines_done]
            logger.info("Re-queuing %d lines for processing" % (missing_lines.shape[0]))
            for line in missing_lines:
                queue.put(line)
        except:
            podi_logging.log_exception()

       
        
    # Tell all workers to shut down when no more data is left to work on
    logger.debug("telling all workers to shut down!")
    for i in range(number_cpus):
        logger.debug("telling all worker %d to shut down!" % (i))
        queue.put((None))

    # Once all command are sent out to the workers, join them to speed things up
    logger.debug("Terminating workers!")
    for p in processes:
        p.terminate()
        p.join(timeout=1)

    results = numpy.copy(shmem_results.to_ndarray()) #).reshape((size_x, size_y)))
    shmem_results.free()

    del shmem_results
    del queue
    del buffer

    return results



def imcombine_subprocess(extension, filelist, shape, operation, queue, verbose,
                         subtract=None, scale=None):

    logger = logging.getLogger("ImCombine")

    #
    # Allocate enough shared momory to hold all frames
    #
    size_x, size_y, n_frames = shape[0], shape[1], shape[2]
    shmem_buffer = SharedMemory(ctypes.c_float, (size_x, size_y, n_frames))
    # multiprocessing.RawArray(ctypes.c_float, size_x*size_y*n_frames) #len(filelist))

    # Extract the shared memory buffer as numpy array to make things easier
    # buffer = shmem_as_ndarray(shmem_buffer).reshape((size_x, size_y, n_frames))
    buffer = shmem_buffer.to_ndarray()

    # Set the full buffer to NaN
    buffer[:,:,:] = numpy.nan

    # Now open all files, look for the right extension, and copy their image data to buffer
    cur_frame = 0
    for file_number in range(len(filelist)):
        filename = filelist[file_number]
        hdulist = pyfits.open(filename)
        for i, ext in enumerate(hdulist):
            if (not is_image_extension(ext)):
                continue
            fppos = ext.name #header['EXTNAME']

            if (not fppos == extension):
                continue

            # Get data for the right extension in this frame
            framedata = ext.data[:,:]
            
            # optionally, apply the scaling and subtraction correction
            if (subtract is not None):
                try:
                    framedata -= float(subtract)
                except ValueError:
                    if (subtract in hdulist[0].header):
                        # print "subtracting",hdulist[0].header[subtract]
                        framedata -= hdulist[0].header[subtract]                
            if (scale is not None):
                try:
                    framedata *= float(scale)
                except ValueError:
                    if (scale in hdulist[0].header):
                        framedata *= hdulist[0].header[scale]         

            # store the (corrected) image data for parallel processing
            buffer[:,:,cur_frame] = framedata
            cur_frame += 1
            break

        ext.data = None
        hdulist.close()
        del hdulist
        if (verbose): stdout_write("\n   Added file %s ..." % (filename))

    if (n_frames > 1):
        # stdout_write("\n   Starting imcombine for real ...")
        combined = imcombine_sharedmem_data(shmem_buffer, operation=operation, sizes=(size_x, size_y, n_frames))
    else:
        logger.debug("Only a single frame contributes to this OTA, skipping combine and copying input to output")
        combined = numpy.array(buffer[:,:,0])

    shmem_buffer.free()

    # put the imcombine'd data into the queue to return them to the main process
    queue.put(combined)

    # and kill this process, returning all its memory
    sys.exit(0)



def imcombine(input_filelist, outputfile, operation, return_hdu=False,
              subtract=None, scale=None, gather_all_otas=True):

    logger = logging.getLogger("ImCombine")

    # First loop over all filenames and make sure all files exist
    filelist = []
    for file in input_filelist:
        if (os.path.isfile(file)):
            filelist.append(file)

    if (len(filelist) <= 0):
        logger.error("No existing files found in input list, hence nothing to do!\n")
        return None

    logger.debug("Stacking the following files:\n -- %s" % ("\n -- ".join(filelist)))

    # elif (len(filelist) == 1):
    #     # stdout_write("Only 1 file to combine, save the hassle and copy the file!\n")
    #     hdulist = pyfits.open(filelist[0])
    #     if (return_hdu):
    #         return hdulist
    #     hdulist.writeto(outputfile, overwrite=True)
    #     return
    
    # Read the input parameters
    # Note that file headers are copied from the first file
    reference_filename = filelist[0]
    ref_hdulist = pyfits.open(reference_filename)

    # Create the primary extension of the output file
    # Copy all headers from the reference HDU
    primhdu = pyfits.PrimaryHDU(header=ref_hdulist[0].header)

    # Add PrimaryHDU to list of OTAs that go into the output file
    out_hdulist = [primhdu]

    ref_hdulist.close()
    del ref_hdulist

    #
    # Compile a list of OTAs that will be in the output frame
    #
    otas_found = []
    ota_sizes = {}
    ref_header = {}
    if (gather_all_otas):
        logger.info("Checking all files to compile comprehensive list of available OTAs")
    master_associations = None
    for fn in filelist:
        hdulist = pyfits.open(fn)
        for ext in hdulist:
            if (is_image_extension(ext)):
                try:
                    otas_found.append(ext.name)
                    if (ext.name not in ota_sizes):
                        ota_sizes[ext.name] = ext.data.shape
                        ref_header[ext.name] = ext.header
                except:
                    podi_logging.log_exception()

        #
        # Collect all individual association tables to create 
        # a master association table combining all input data
        #
        assoc_table = podi_associations.read_associations(hdulist)
        if (assoc_table is None):
            logger.info("No association data available")
            assoc_table = {'input_simple': [fn]}
        if (master_associations is None):
            master_associations = assoc_table
        else:
            master_associations = podi_associations.collect_reduction_files_used(
                master_associations, assoc_table)

        hdulist.close()
        del hdulist
        if (not gather_all_otas):
            break

    otas_to_combine = set(otas_found)
    ota_counter = collections.Counter(otas_found)
    #print otas_to_combine
    #print ota_counter


    #
    # Now loop over all extensions and compute the mean
    #
    # for cur_ext in range(0, len(ref_hdulist)):
    for cur_ext, extname in enumerate(otas_to_combine):


        data_blocks = []
        # Check what OTA we are dealing with
        # if (not is_image_extension(ref_hdulist[cur_ext])):
        #     continue
        # ref_fppos = ref_hdulist[cur_ext].name #header['EXTNAME']

        #stdout_write("\rCombining frames for OTA %s (#% 2d/% 2d) ..." % (ref_fppos, cur_ext+1, len(ref_hdulist)))
        #logger.debug("Combining frames for OTA %s (#% 2d/% 2d) ..." % (ref_fppos, cur_ext+1, len(ref_hdulist)))
        logger.info("Combining frames for OTA %s (#% 2d/% 2d) ..." % (extname, cur_ext+1, len(otas_to_combine)))

        #
        # Add some weird-looking construct to move the memory allocation and actual 
        # imcombine into a separate process. This has to do with if/how/when python releases 
        # memory (or not), causing a massive short-term memory leak.
        # With the additional process, the memory is onwed by the other process and the memory
        # is freed once we destroy this helper process.
        #
        return_queue = multiprocessing.JoinableQueue()
        #worker_args=(ref_fppos, filelist, ref_hdulist[cur_ext].data.shape, operation, return_queue, verbose)

        kw_args = {
            'extension': extname,
            'filelist':  filelist, 
            'shape':     (ota_sizes[extname][0], ota_sizes[extname][1], ota_counter[extname]), #ref_hdulist[cur_ext].data.shape,
            'operation': operation, 
            'queue':     return_queue, 
            'verbose':   verbose,
            'subtract':  subtract,
            'scale':     scale,
        }
        # p = multiprocessing.Process(target=imcombine_subprocess, args=worker_args)
        p = multiprocessing.Process(target=imcombine_subprocess, kwargs=kw_args)
        p.start()
        combined = return_queue.get()
        p.terminate()
        del p

        if (verbose): stdout_write(" done, creating fits extension ...")
        logger.debug("done with computing, creating fits extension ...")
        # Create new ImageHDU, insert the imcombined's data and copy the 
        # header from the reference frame
        hdu = pyfits.ImageHDU(data=combined, header=ref_header[extname])

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

    logger.debug("Adding association data to output file")

    out_hdulist[0].header["NCOMBINE"] = (len(filenames_only), "number of combined files")

    columns = [
        pyfits.Column(name='filenumber', format="I4", array=range(1,len(filenames_only)+1)),
        pyfits.Column(name='filename', format="A100", array=filenames_only),
        ]
    coldefs = pyfits.ColDefs(columns)
    tablehdu = pyfits.BinTableHDU.from_columns(coldefs)
    tablehdu.header["EXTNAME"] = "FILELIST"
    out_hdulist.append(tablehdu)

    assoc_tbhdu = podi_associations.create_association_table(
        master_associations)
    out_hdulist.append(assoc_tbhdu)

    #
    # All work done now, prepare to return the data or write it to disk
    #
    out_hdu = pyfits.HDUList(out_hdulist)
    if (not return_hdu and outputfile != None):
        logger.debug(" writing results to file %s ..." % (outputfile))
        clobberfile(outputfile)
        try:
            out_hdu.writeto(outputfile, overwrite=True, checksum=True)
        except TypeError:
            # this most likely is this error:
            # TypeError: object of type 'NoneType' has no len()
            # related to the checksum calculation
            clobberfile(outputfile)
            out_hdu.writeto(outputfile, overwrite=True)
        except pyfits.VerifyError:
            logger.warning("Encountered FITS verification error, writing anyway")
            try:
                out_hdu.writeto(outputfile, overwrite=True, checksum=False, output_verify='ignore')
            except:
                raise
        except:
            raise
        out_hdu.close()
        del out_hdu
        del out_hdulist
        stdout_write(" done!\n")
    elif (return_hdu):
        logger.debug(" returning HDU for further processing ...")
        return out_hdu
    else:
        logger.debug(" couldn't write output file, no filename given!")

    return None

if __name__ == "__main__":

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)

    outputfile = get_clean_cmdline()[1]

    filelist = handle_filelists(get_clean_cmdline()[2:])

    operation = cmdline_arg_set_or_default("-op", "mean")

    subtract = cmdline_arg_set_or_default("-subtract", None)
    scale = cmdline_arg_set_or_default("-scale", None)


    print("Combining into %s:\n -- %s" % (outputfile, "\n -- ".join(filelist)))

    imcombine(filelist, outputfile, operation, subtract=subtract, scale=scale)

    podi_logging.shutdown_logging(options)
