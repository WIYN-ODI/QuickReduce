#!/usr/bin/env python3


import os
import sys
import pyfits
import numpy
import logging

from astLib import astWCS

# QR imports
from podi_definitions import is_image_extension
import podi_fitskybackground



def scale_subtract_background_model(in_filename, ic_file, out_filename,
                                    per_ota=True,
                                    remove_gradient=None,
                                    reuse_samples=False,
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

    logger.info("Finding optimal scaling for background model")
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
                min_found=200, boxwidth=30,
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

    # numpy.savetxt("ratios.dmp", all_ratios)


    #
    # now correct data
    #
    logger.info("Subtracting background model")
    bgsub_samples = None
    if (per_ota):
        for extname in ota_ratios:
            scaling = numpy.nanmedian(ota_ratios[extname])
            imghdu[extname].data -= ic_hdu[extname].data * scaling
            logger.debug("Scaling model by %f for %s" % (scaling, extname))
    else:
        global_ratio = numpy.nanmedian(all_ratios)
        for ext in imghdu:
            if (not is_image_extension(ext) or not ext.name in ic_hdu):
                continue
            ext.data -= ic_hdu[ext.name].data * global_ratio
            logger.debug("Scaling model by %f for %s" % (global_ratio, ext.name))

            numpy.savetxt("all_sky.dmp", all_sky_samples)
            numpy.savetxt("all_models.dmp", all_model_samples)
            all_sky_samples[:,4] -= global_ratio * all_model_samples[:,4]
            numpy.savetxt("all_sky_sub.dmp", all_sky_samples)

    #
    # save output
    #
    if (out_filename is not None):
        logger.info("Saving output")
        imghdu.writeto(out_filename, clobber=True)

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

