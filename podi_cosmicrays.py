#!/usr/bin/env python

"""

Wrapper around the LACosmics python implementation by Malte Tewes. Also provides
a stand-alone funtionality, based on Malte's demo program.

"""

import ext_cosmics as cosmics
import pyfits
import numpy
from podi_definitions import *
import itertools
import podi_sitesetup as sitesetup
import time
import multiprocessing
import podi_cython



def remove_cosmics(hdu, n_iterations=4, method="cy"):

    gain = hdu.header['GAIN']
    readnoise = 8 #hdu.header['RON']

    corrected = hdu.data.copy()
    mask = numpy.zeros_like(hdu.data)

    # Break up each entire OTA into its cells to see if this is faster
    for cx, cy in itertools.product(range(8),repeat=2):
        # print "working on cell",cx, cy

        x1, x2,y1, y2 = cell2ota__get_target_region(cx, cy, binning=1)
        cell_data = hdu.data[y1:y2, x1:x2]

        if (method == "py"):
            c = cosmics.cosmicsimage(cell_data, 
                                     gain=gain, readnoise=readnoise, 
                                     sigclip = 5.0, sigfrac = 0.3, objlim = 5.0,
                                     verbose=False)
            c.run(maxiter=n_iterations)

            # # Re-insert the data into the full frame
            corrected[y1:y2, x1:x2] = c.cleanarray
            mask[y1:y2, x1:x2][c.mask] = 1
        elif (method == "cy"):

            # Run cosmic ray removal
            crj = podi_cython.lacosmics(cell_data.astype(numpy.float64), 
                                        gain=gain, readnoise=readnoise, 
                                        niter=int(n_iterations),
                                        sigclip = 5.0, sigfrac=0.3, objlim=5.0,
                                        verbose=False)
            cell_cleaned, cell_mask, cell_saturated = crj

            # Re-insert data into full OTA
            corrected[y1:y2, x1:x2] = cell_cleaned
            mask[y1:y2, x1:x2] = cell_mask
            
        else:
            corrected[y1:y2, x1:x2] = cell_data
            print "Unknown method"

    return corrected, mask

    # # Prepare all structures
    # c = ext_lacosmics.cosmicsimage(hdu.data, 
    #                                gain=gain, readnoise=readnoise, 
    #                                sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)

    # # Run the cosmic removal
    # c.run(maxiter=4)

    # return c.cleanarray, c.mask


def removecosmics_mp(input_queue, return_queue):

    while (True):
        task = input_queue.get()
        if (task == None):
            input_queue.task_done()
            break

        hdu, niter, id, method = task
        # print hdu
        # print niter
        # print id
        # print method

        if ('OTA' in hdu.header):
            print "Removing cosmics from OTA",hdu.header['OTA'],"using method",method
        else:
            print "removing cosmics"

        fixed, mask = remove_cosmics(hdu, niter, method)

        ota = hdu.header['OTA']
        pyfits.HDUList([pyfits.PrimaryHDU(data=fixed, header=hdulist[i].header)]).writeto("cleaned_OTA%d.%s.fits" % (ota, method), clobber=True)
        pyfits.HDUList([pyfits.PrimaryHDU(data=mask, header=hdulist[i].header)]).writeto("mask_OTA%d.%s.fits" % (ota, method), clobber=True)

        # Send return data
        ret = (fixed, mask, id)
        return_queue.put(ret)

        input_queue.task_done()
        continue

    return


if __name__ == "__main__":

    if (cmdline_arg_isset("-standalone")):
        inputfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(inputfile)

        #import yappi
        #yappi.start()

        method = cmdline_arg_set_or_default("-method", "cy")
        niter = int(cmdline_arg_set_or_default("-niter", 4))

        job_queue = multiprocessing.JoinableQueue()
        return_queue = multiprocessing.JoinableQueue()

        n_jobs = 0
        for i in range(len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue

            task = (hdulist[i], niter, i, method)
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

    else:
        print __doc__
        print cosmics.__doc__
