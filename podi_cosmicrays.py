#!/usr/bin/env python

"""

Wrapper around the LACosmics python implementation by Malte Tewes. Also provides
a stand-alone funtionality, based on Malte's demo program.

"""

import ext_cosmics as cosmics
import pyfits
import numpy
from podi_definitions import *
from podi_commandline import *
import itertools
import podi_sitesetup as sitesetup
import time
import multiprocessing
import podi_cython
import podi_collectcells
import podi_logging
import logging


def remove_cosmics(data, n_iterations=4, method="cy", 
                   gain=1.5, readnoise=9,
                   sigclip = 5.0, sigfrac = 0.3, objlim = 5.0,
                   verbose=False,
                   saturation_limit=65000,
                   binning=1,
):


    corrected = data.copy()
    mask = numpy.zeros_like(data)

    logger = logging.getLogger("CosmicRayX(%s)" % (method))
    logger.debug("Parameters: gain=%.3f readnoise=%.3f sigclip=%.2f, sigfrac=%.2f objlim=%.2f sat-limit=%d" % (
        gain, readnoise, sigclip, sigfrac, objlim, saturation_limit))

    # Break up each entire OTA into its cells to see if this is faster
    for cx, cy in itertools.product(range(8),repeat=2):
        if (verbose): stdout_write("\rworking on cell %d,%d" % (cx, cy))

        x1, x2,y1, y2 = cell2ota__get_target_region(cx, cy, binning=binning)
        cell_data = data[y1:y2, x1:x2]

        if (method == "py"):
            saturation_limit = saturation_limit if saturation_limit > 0 else 1e10
            c = cosmics.cosmicsimage(cell_data, 
                                     gain=gain, readnoise=readnoise, 
                                     sigclip=sigclip, sigfrac=sigfrac, objlim=objlim,
                                     verbose=False)
            c.run(maxiter=n_iterations)

            # # Re-insert the data into the full frame
            corrected[y1:y2, x1:x2] = c.cleanarray
            mask[y1:y2, x1:x2][c.mask] = 1
            logger.debug("Found %d cosmics in cell %d,%d" % (numpy.sum(c.mask), cx, cy))

        elif (method == "cy"):
            # Run cosmic ray removal
            crj = podi_cython.lacosmics(cell_data.astype(numpy.float64), 
                                        gain=gain, readnoise=readnoise, 
                                        niter=int(n_iterations),
                                        sigclip=sigclip, sigfrac=sigfrac, objlim=objlim,
                                        saturation_limit=saturation_limit,
                                        verbose=False)
            cell_cleaned, cell_mask, cell_saturated = crj

            # Re-insert data into full OTA
            corrected[y1:y2, x1:x2] = cell_cleaned
            mask[y1:y2, x1:x2] = cell_mask
            
            logger.debug("Found %d cosmics in cell %d,%d" % (numpy.sum(cell_mask>0), cx, cy))
        else:
            corrected[y1:y2, x1:x2] = cell_data
            logger.warning("Unknown Crj-X method: %s" % (method))

    return corrected, mask


def removecosmics_mp(input_queue, return_queue):

    while (True):
        task = input_queue.get()
        if (task == None):
            input_queue.task_done()
            break

        data, niter, id, method, verbose, setup = task
        gain, readnoise, sigclip, sigfrac, objlim, saturation_limit, binning = setup

        print "Starting cosmic removal on one extension"

        fixed, mask = remove_cosmics(data, niter, method, 
                                     gain=gain, readnoise=readnoise,
                                     sigclip=sigclip, sigfrac=sigfrac, objlim=objlim,
                                     saturation_limit=saturation_limit,
                                     binning=binning,
                                     verbose=verbose)

        # Send return data
        ret = (fixed, mask, id)
        return_queue.put(ret)

        input_queue.task_done()
        continue

    return


if __name__ == "__main__":

    if (len(sys.argv) <= 1):
        print __doc__
        print cosmics.__doc__

    else:

        options = podi_collectcells.read_options_from_commandline()
        podi_logging.setup_logging(options)

        inputfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(inputfile)
        
        #import yappi
        #yappi.start()

        method = cmdline_arg_set_or_default("-method", "cy")
        niter = int(cmdline_arg_set_or_default("-niter", 4))
        verbose = cmdline_arg_isset("-verbose")

        sigclip = float(cmdline_arg_set_or_default("-sigclip", 5.0))
        sigfrac = float(cmdline_arg_set_or_default("-sigfrac", 0.3))
        objlim = float(cmdline_arg_set_or_default("-objlim", 5.0))
        saturation_limit = float(cmdline_arg_set_or_default("-saturate", 60000))

        job_queue = multiprocessing.JoinableQueue()
        return_queue = multiprocessing.JoinableQueue()

        n_jobs = 0
        for i in range(len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue

            gain = hdulist[i].header['GAIN'] if 'GAIN' in hdulist[i].header else 1.5
            readnoise = hdulist[i].header['RDNOISEE'] if 'RDNOISEE' in hdulist[i].header else 8.5
            binning = hdulist[0].header['BINNING'] if 'BINNING' in hdulist[0].header else 1

            setup = gain, readnoise, sigclip, sigfrac, objlim, saturation_limit, binning
            task = (hdulist[i].data, niter, i, method, verbose, setup)
            job_queue.put(task)
            n_jobs += 1


        # Now spawn of the workers
        processes = []
        worker_args = (job_queue, return_queue)
        for i in range(sitesetup.number_cpus):
            p = multiprocessing.Process(target=removecosmics_mp, args=worker_args)
            p.start()
            processes.append(p)
            job_queue.put(None)
            time.sleep(0.01)
            
            # ignore the mask for now

        # Put results back into file
        for i in range(n_jobs):
            ret = return_queue.get()
            fixed, mask, id = ret
            hdulist[id].data = fixed
            return_queue.task_done()

        # Join processes
        for p in processes:
            p.join()

        outputfile = get_clean_cmdline()[2]
        hdulist.writeto(outputfile, clobber=True)

        #yappi.get_thread_stats().print_all()
        #yappi.get_func_stats().debug_print() #print_all()

        podi_logging.shutdown_logging(options)

