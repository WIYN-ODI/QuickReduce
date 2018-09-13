#!/usr/bin/env python3


import os
import sys
import pyfits
import numpy
import logging

from astLib import astWCS
import scipy.spatial
import itertools

# QR imports
from podi_definitions import is_image_extension
import podi_fitskybackground
import podi_photflat

debug = False

def scale_subtract_background_model(in_filename, ic_file, out_filename,
                                    per_ota=True,
                                    remove_gradient=None,
                                    reuse_samples=False,
                                    twod_model=True,
                                    logger=None):

    if (logger is None):
        logger = logging.getLogger("ScaleSubtractBG")

    if (type(in_filename) == str):
        imghdu = pyfits.open(in_filename)
    else:
        imghdu = in_filename

    if (type(ic_file) is str):
        ic_hdu = pyfits.open(ic_file)
    else:
        ic_hdu = ic_file

    # now go over each OTA, and find sky-measurement boxes
    # ext_skylevels = imghdu['SKYLEVEL']

    final_grid_sky = {}
    pixel_sampling = 500
    _y, _x = numpy.indices((9,9))
    grid_x = _x * pixel_sampling
    grid_y = _y * pixel_sampling
    # if (debug): print(grid_x)
    logger.info("Finding optimal scaling for background model")
    n_sky_samples = 500
    smoothing_length = 8.  # 8 arcmin
    all_ratios = None
    all_sky_samples = None
    all_model_samples = None
    ota_ratios = {}
    for ext in imghdu:
        if (not is_image_extension(ext)):
            continue
        # also check if we have this extension in the IC frame
        try:
            ic_ext = ic_hdu[ext.name]
        except:
            continue

        logger.debug("Working on %s" % (ext.name))
        wcs = astWCS.WCS(ext.header, mode='pyfits')

        sky_samples = numpy.array(
            podi_fitskybackground.sample_background(
                data=ext.data, wcs=wcs,
                starcat=None,
                min_found=n_sky_samples, boxwidth=30,
                min_box_spacing=3))
        logger.debug(sky_samples.shape)

        # now repeat measurements in the IC data
        ic_samples = numpy.array(
            podi_fitskybackground.sample_background(
                data=ic_ext.data,
                box_center=sky_samples[:, 2:4],
                wcs=None, starcat=None,
            )
        )
        logger.debug(ic_samples.shape)

        ratio = sky_samples[:,4] / ic_samples[:, 4]
        logger.debug("%s %s %s" % (ratio.shape, numpy.nanmedian(ratio), numpy.nanstd(ratio)))

        all_ratios = ratio if all_ratios is None else numpy.append(all_ratios, ratio)
        if (all_sky_samples is None):
            all_sky_samples = sky_samples
            all_model_samples = ic_samples
        else:
            all_sky_samples = numpy.append(all_sky_samples, sky_samples, axis=0)
            all_model_samples = numpy.append(all_model_samples, ic_samples, axis=0)

        if (per_ota):
            ota_ratios[ext.name] = ratio

        # also calculate the corners for the final grid
        final_grid_sky[ext.name] = numpy.array(wcs.pix2wcs(grid_x.ravel(), grid_y.ravel()))\
            .reshape((grid_x.shape[0], grid_x.shape[1], 2))
        #print(final_grid_sky[ext.name])
        #numpy.savetxt("gridraw", final_grid_sky[ext.name].reshape((-1,2)))
        #break

    if (debug):
        numpy.savetxt("ratios.dmp", all_ratios)
        numpy.savetxt("bgsamples.sky", all_sky_samples)

        all_model_samples[:, 0:2] = all_sky_samples[:, 0:2]
        numpy.savetxt("bgsamples.model", all_model_samples)

        numpy.savetxt("all_sky.dmp", all_sky_samples)
        numpy.savetxt("all_models.dmp", all_model_samples)

    all_model_samples[:, 0:2] = all_sky_samples[:, 0:2]

    #
    # now correct data
    #
    logger.info("Subtracting background model")
    bgsub_samples = None
    per_ota=False
    if (per_ota):
        for extname in ota_ratios:
            scaling = numpy.nanmedian(ota_ratios[extname])
            imghdu[extname].data -= ic_hdu[extname].data * scaling
            logger.debug("Scaling model by %f for %s" % (scaling, extname))
    else:
        global_ratio = numpy.nanmedian(all_ratios)
        logger.info("Using global sky scaling = %f" % (global_ratio))
        for ext in imghdu:
            if (not is_image_extension(ext) or not ext.name in ic_hdu):
                continue
            ext.data -= ic_hdu[ext.name].data * global_ratio
            logger.debug("Scaling model by %f for %s" % (global_ratio, ext.name))

    all_sky_samples[:,4] -= global_ratio * all_model_samples[:,4]
    if (debug):
        numpy.savetxt("all_sky_sub.dmp", all_sky_samples)
        print(all_sky_samples[:10, 4])

    if (twod_model):
        logger.info("Running 2-d model")
        # print(final_grid_sky)
        projected_sky = all_sky_samples.copy()
        median_declination = numpy.median(projected_sky[:, 1])
        median_ra = numpy.median(projected_sky[:,0])
        cos_declination = numpy.cos(numpy.radians(median_declination))
        # projected_sky[:,0] *= cos_declination
        projected_sky[:,0] = (projected_sky[:,0] - median_ra) * numpy.cos(numpy.radians(projected_sky[:,1])) + median_ra
        if (debug): numpy.savetxt("all_sky_proj.dmp", projected_sky)


        coord_tree = scipy.spatial.cKDTree(projected_sky[:, 0:2])
        max_samples_return = 1500
        gridded_sky = []
        skyframe = [pyfits.PrimaryHDU(header=imghdu[0].header)]

        for extname in final_grid_sky:
            radec_grid = final_grid_sky[extname]
            # numpy.savetxt("grid_%s" % (extname), radec_grid.reshape((-1,2)))
            logger.info("Calculating background model for %s" % (extname))

            # equally deproject all grid boxes
            radec_grid = final_grid_sky[extname]
            grid_ra = radec_grid[:,:,0]
            grid_dec = radec_grid[:,:,1]
            grid_ra = (grid_ra - median_ra) * numpy.cos(numpy.radians(grid_dec)) + median_ra
            # radec_grid[:,0] = (radec_grid[:,0] - median_ra) * numpy.cos(numpy.radians(radec_grid[:,1])) + median_ra
            # print(extname, radec_grid[:3,:3,:])
            # numpy.savetxt("grid_%s" % (extname), numpy.array([
            #     grid_ra.flatten(), grid_dec.flatten()]).T)

            local_sky_array = numpy.zeros(grid_x.shape)

            for x,y in itertools.product(range(radec_grid.shape[0]), range(radec_grid.shape[1])):
                gra, gdec = grid_ra[y,x], grid_dec[y,x] #radec_grid[y,x]
                # print(gra,gdec)

                d, i = coord_tree.query([gra,gdec],
                    distance_upper_bound=smoothing_length/60.,
                    k=max_samples_return,
                    p=2)
                valid_indices = numpy.isfinite(d)
                idx = i[valid_indices]

                local_sky_samples = projected_sky[idx]
                local_sky = numpy.nanmedian(local_sky_samples[:, 4])
                # TODO: outlier rejection

                gridded_sky.append([gra, gdec, local_sky, local_sky_samples.shape[0]])
                local_sky_array[x,y] = local_sky

            # numpy.savetxt("grid_%s" % (extname), radec_grid.reshape((-1,2)))

            logger.info("Expanding %s to full-res" % (extname))
            fullframe = podi_photflat.expand_to_fullres(local_sky_array, blocksize=pixel_sampling,
                                            mag2flux=False)
            skyframe.append(pyfits.ImageHDU(header=imghdu[extname].header,
                                            data=fullframe))

            imghdu[extname].data -= fullframe
            # print("\n"*3)

        gridded_sky = numpy.array(gridded_sky)
        if (debug): numpy.savetxt("sky.grid", gridded_sky)

        #logger.info("Saving skyframe")
        sky_hdu = pyfits.HDUList(skyframe)
        #sky_hdu.writeto("skyframe.fits", clobber=True)

    #
    # save output
    #
    if (out_filename is not None):
        logger.info("Saving sky-model output to %s" % (out_filename))
        imghdu.writeto(out_filename, clobber=True)

    if (twod_model):
        return imghdu, sky_hdu

    return imghdu

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("IC-back")

    in_filename = sys.argv[1]
    ic_file = sys.argv[2]
    out_filename = sys.argv[3]

    try:
        per_ota = (sys.argv[4] == "ota")
    except:
        per_ota = False


    scale_subtract_background_model(
        in_filename=in_filename,
        ic_file=ic_file,
        out_filename=out_filename,
        per_ota=per_ota)

