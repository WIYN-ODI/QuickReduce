#!/usr/bin/env python3

"""

Wrapper around the LACosmics python implementation by Malte Tewes. Also provides
a stand-alone funtionality, based on Malte's demo program.

"""

import ext_cosmics as cosmics
import astropy.io.fits as pyfits
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
        if (task is None):
            input_queue.task_done()
            break

        data, niter, id, method, verbose, setup = task
        gain, readnoise, sigclip, sigfrac, objlim, saturation_limit, binning = setup

        print("Starting cosmic removal on one extension")

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


def csv_to_list(x):
    out = []
    items = x.split(",")
    for item in items:
        out.append(float(item))
    return out

fwhm_to_sigma = 2 * math.sqrt(2 * math.log(2))

if __name__ == "__main__":

    if (len(sys.argv) <= 1):
        print(__doc__)
        print(cosmics.__doc__)

    elif (cmdline_arg_isset("-simple")):

        options = podi_collectcells.read_options_from_commandline()
        podi_logging.setup_logging(options)

        inputfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(inputfile)

        bglevel = float(cmdline_arg_set_or_default('-bg', 0.0))

        for i in range(len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue

            hdulist[i].data += bglevel
            crj = podi_cython.lacosmics(hdulist[i].data.astype(numpy.float64), 
                                        gain=3.3, readnoise=20, 
                                        niter=4,
                                        sigclip=5.0, sigfrac=0.5, objlim=5.0,
                                        saturation_limit=25000,
                                        verbose=False)
            cell_cleaned, cell_mask, cell_saturated = crj

            hdulist[i].data = cell_cleaned - bglevel


        outputfile = get_clean_cmdline()[2]
        hdulist.writeto(outputfile, overwrite=True)
        podi_logging.shutdown_logging(options)

    elif (cmdline_arg_isset("-check")):

        bglevels = csv_to_list(cmdline_arg_set_or_default("-bg", '0'))
        peaks = csv_to_list(cmdline_arg_set_or_default("-peak", '100'))
        cosmicfluxes = csv_to_list(cmdline_arg_set_or_default("-crflux", '1000'))
        fwhms = csv_to_list(cmdline_arg_set_or_default("-fwhm", '0.4,0.6'))
        nrandom = int(cmdline_arg_set_or_default("-n", 100))
        imagesize = int(cmdline_arg_set_or_default("-imgsize", 101))
        if (imagesize % 2 == 0):
            imagesize += 1
        centering = float(cmdline_arg_set_or_default("-center", 0.2))
        niter = int(cmdline_arg_set_or_default("-niter", 1))

        sigclip = float(cmdline_arg_set_or_default("-sigclip", 5.0))
        sigfrac = float(cmdline_arg_set_or_default("-sigfrac", 0.3))
        objlim  = float(cmdline_arg_set_or_default("-objlim", 5.0))
                                                        
        readnoise = 8
        gain = 1.3

        primhdu = pyfits.PrimaryHDU()
        hdulist_in = [primhdu]
        hdulist_out = [primhdu]

        # Allocate some memory to hold all the individual postage stamps
        n_thumbs = int(math.ceil(math.sqrt(nrandom)))
        tile_in = numpy.empty((n_thumbs*imagesize, n_thumbs*imagesize))
        tile_out = numpy.empty((n_thumbs*imagesize, n_thumbs*imagesize))

#        centering = 0.5

        for fwhm in fwhms:
            fwhm_pixels = fwhm / 0.11
            #imagesize = int(20 * fwhm_pixels)+1
            center = 0.5*(imagesize-1)

            for peak in peaks:
            
                img = numpy.zeros((imagesize, imagesize))
                print(img.shape)

                # Now convolve the source with the gaussian PSF
                gauss_sigma = fwhm_pixels / fwhm_to_sigma
                _x, _y = numpy.indices(img.shape)
                _x -= center
                _y -= center
                _r = numpy.hypot(_x,_y)
                psf = numpy.exp(-_r**2/(2*gauss_sigma**2)) * peak

                # img[center,center] = peak

                # try:
                #     psf = scipy.ndimage.filters.gaussian_filter(
                #         input=img,
                #         sigma = gauss_sigma,
                #         order = 0,
                #         mode = 'constant',
                #         cval = 0,
                #         truncate = 10
                #     )
                # except:
                #     psf = scipy.ndimage.filters.gaussian_filter(
                #         input=img,
                #         sigma = gauss_sigma,
                #         order = 0,
                #         mode = 'constant',
                #         cval = 0,
                #     )
                
                for bg in bglevels:

                    #
                    # Create a fake image with the source at the center
                    #
                    with_bg = psf + bg

                    for cosmicflux in cosmicfluxes:

                        # Reset the output buffer
                        tile_in[:,:] = numpy.nan
                        tile_out[:,:] = numpy.nan

                        print(bg, peak, cosmicflux, fwhm)
                        
                        #
                        # Add poisson noise to the image
                        #
                        noise_level = numpy.sqrt(with_bg * gain + readnoise**2)
                        print(noise_level.shape)

                        # Now pick positions for the cosmics
                        # cosmics will occupy a fraction of the central area of 
                        # each postage stamp
                        cr_pos = numpy.array(numpy.random.random((nrandom,2))
                                             * centering * imagesize
                                             + 0.5*(1.-centering)*imagesize,
                                             dtype=int)

                        for n in range(nrandom):

                            # Compute some noise and add it to the PSF with background
                            noise = numpy.random.randn(noise_level.shape[0], noise_level.shape[1]) * noise_level
                            noisy_image = with_bg + noise

                            # Add the cosmic
                            noisy_image[cr_pos[n,0], cr_pos[n,1]] += cosmicflux

                            # Save image for later
                            tile_x = n % n_thumbs
                            tile_y = int(math.floor(n/n_thumbs))

                            tile_in[tile_x*imagesize:(tile_x+1)*imagesize, 
                                    tile_y*imagesize:(tile_y+1)*imagesize] = noisy_image

                            flux_in = numpy.sum(noisy_image)

                            #
                            # Run CRJ detection/rejection
                            #
                            if (not cmdline_arg_isset("-nx")):
                                crj = podi_cython.lacosmics(noisy_image.astype(numpy.float64), 
                                                            gain=gain, 
                                                            readnoise=readnoise, 
                                                            niter=niter, #int(n_iterations),
                                                            sigclip=sigclip, 
                                                            sigfrac=sigfrac, 
                                                            objlim=objlim,
                                                            saturation_limit=65000, #saturation_limit,
                                                            verbose=cmdline_arg_isset("-verbose"))
                                fixed, mask, cell_saturated = crj
                            else:

                                for i in range(niter):
                                    crj = podi_cython.lacosmics(noisy_image.astype(numpy.float64), 
                                                                gain=gain, 
                                                                readnoise=readnoise, 
                                                                niter=1, #niter, #int(n_iterations),
                                                                sigclip=sigclip, 
                                                                sigfrac=sigfrac, 
                                                                objlim=objlim,
                                                                saturation_limit=65000, #saturation_limit,
                                                                verbose=False)
                                    fixed, mask, cell_saturated = crj
                                    noisy_image[:,:] = fixed[:,:]

                            # fixed, mask = remove_cosmics(noisy_image, n_iterations=1, method='cy', 
                            #                              gain=gain, 
                            #                              readnoise=readnoise,
                            #                              # These are the important parameters
                            #                              sigclip = 5.0, 
                            #                              sigfrac = 0.3, 
                            #                              objlim = 5.0,
                            #                              # ...
                            #                              saturation_limit=65000.,
                            #                              binning=1,
                            #                              verbose=True)



                            tile_out[tile_x*imagesize:(tile_x+1)*imagesize, 
                                    tile_y*imagesize:(tile_y+1)*imagesize] = fixed

                            # Compute how much flux we removed vs. 
                            # how much flux is in the cosmic
                            flux_out = numpy.sum(fixed)
                            print("%7.1f %7.1f --> %7.1f == %6.3f" % (
                                flux_in, flux_out, flux_in-flux_out, (flux_in-flux_out)/cosmicflux
                            ))

                        #all_images = numpy.append(all_images, noisy_image, axis=0)
                        imghdu = pyfits.ImageHDU(data=tile_in.copy())
                        imghdu.header['OBJECT'] = "Input: bg %d, flux %.1f, fwhm=%.2f" % (bg, peak, fwhm)
                        hdulist_in.append(imghdu)

                        imghdu_fixed = pyfits.ImageHDU(data=tile_out.copy())
                        imghdu.header['OBJECT'] = "CRJ fixed: bg %d, flux %.1f, fwhm=%.2f" % (bg, peak, fwhm)
                        hdulist_out.append(imghdu_fixed)


                        #
                        # Check if the cosmic was found 
                        #

                    # for cosmicflux
                # for bg
            # for peak

        # next fwhm

        hdulist_in = pyfits.HDUList(hdulist_in)
        hdulist_in.writeto("crj_debug_in.fits", overwrite=True)

        hdulist_out = pyfits.HDUList(hdulist_out)
        hdulist_out.writeto("crj_debug_out.fits", overwrite=True)
            

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
        hdulist.writeto(outputfile, overwrite=True)

        #yappi.get_thread_stats().print_all()
        #yappi.get_func_stats().debug_print() #print_all()

        podi_logging.shutdown_logging(options)

