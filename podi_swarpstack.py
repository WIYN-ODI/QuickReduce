#!/usr/bin/env python3
#
# Copyright 2012-2013 Ralf Kotulla & WIYN Observatory
#                     University of Wisconsin - Milwaukee & Madison
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestions on how to improve the code or its 
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

Function
--------------------------------------------------------------------------------

**podi_swarpstack.py** handles most of the most common stacking problems. It 
supports a number of use-cases, from simle stacking a bunch of files, adding
more frames to an additional stack, as well as matching a new frame to a 
reference stack.

The general procedure is to create a single rectified image for each of the
input frames first. In a second step, swarp combines all these rectified images
into the final step. While this is overkill for stellar fields, this approach
ensures better sky background subtraction across OTAs in the case of large
galaxies or objects. 

Before each run, swarpstack determines the final size of the resulting stacked
frame. Each of the input frames is then rectified, distortion corrected and
interpolated onto this final grid, eliminating double-interpolations.



How to run podi_swarpstack
--------------------------------------------------------------------------------

Simple stacking of a bunch of frames

    podi_swarpstack.py output_stack.fits file1.fits file2.fits


Adding frames to an already existing stack

    podi_swarpstack.py -add output_stack.fits file3.fits

Matching a new frame to a reference frame

    podi_swarpstack.py -reference=reference.fits output_stack.fits file*.fits



Additional command line options
--------------------------------------------------------------------------------

* **-bgsub**

  Activate background subtraction. Without this frame no background subtraction
  is performed. While this is the safest options for some cases it leads to 
  rather ratty looking stacks if the sky-level was at all variable during the 
  exposure (e.g. with moon, twilight, scattered light)


* **-reusesingles**

  By default, the single recified images are re-created when swarpstack is being
  run. With this option swarpstack checks if a recitified image already exists
  and uses the existing file instead of re-creating a new one. This is
  particularly helpful and in fact the default if adding new frames to an
  already existing stack.

  
* **-dimension=X**

  If given without the -bgsub parameter this option has no effect. If -bgsub is
  enabled, the background filtering options are adjusted automatically to avoid
  over-subtracting objects of the specified dimenson or smaller, hence avoiding
  subtracting the galacy continuum as background fluctuation. The value X
  specifies the galaxy/source dimension in arc-minutes.

* **-pixelscale=X**

  Chooses the specified pixel scale for the output stack, with X given in 
  arcsec/pixel.

* **-nonsidereal=a,b,c**

  Apply a non-sidereal motion correction to all input frames before stacking.
  a is the motion in dRA*cos(dec) in arcseconds/hour, b is the motion dDec, also 
  in arcsec/hour, and c is the MJD reference, either a valid FITS frame or the
  MJD as floating point number

"""


from __future__ import print_function
import os
import sys
import astropy.io.fits as pyfits
import subprocess
import math
import shutil

from podi_commandline import *
from podi_definitions import *
#from podi_collectcells import *
import podi_sitesetup as sitesetup
import podi_illumcorr
import multiprocessing
from podi_fitskybackground import sample_background_using_ds9_regions
import podi_associations
import podi_photflat

import podi_logging
import logging
import socket
import tempfile
import shutil
import warnings
import time
import itertools
import bottleneck
import ephem
import ic_background
from podi_focalplanelayout import FocalPlaneLayout

# try:
#     # sys.path.append("/work/quickreduce/merge_master_py3/test/")
#     # print(sys.path)
#     import ephemerides
#     import podi_ephemerides
#     print("Added the test durectory for ephems")
# except:
#     pass





def mp_prepareinput(input_queue, output_queue, swarp_params, options, apf_data=None):

    while (True):

        cmd = input_queue.get()
        if (cmd is None):
            input_queue.task_done()
            break

        (input_file, fileid) = cmd
        logger = logging.getLogger("MP-Prep( %s )" % (os.path.basename(input_file)))
        
        ret = {
            "master_reduction_files": {},
            "corrected_file": None,
            'exptime': 0,
            'mjd_obs_start': 0,
            'mjd_obs_end': 0,
            'nonsidereal-dradec': None,
        }

        try:
            hdulist = pyfits.open(input_file)
        except IOError:
            logger.error("Can't open file %s" % (input_file))
            output_queue.put(None)
            input_queue.task_done()
            continue

        mjd_obs_start = hdulist[0].header['MJD-OBS']
        exptime = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header else \
                  hdulist[0].header['EXPTIME']
        mjd_obs_end = hdulist[0].header['MJD-OBS'] + (exptime/86400.)

        # Save these values for the return queue
        ret['exptime'] = exptime
        ret['mjd_obs_start'] = mjd_obs_start
        ret['mjd_obs_end'] = mjd_obs_end

        corrected_filename = None

        # Keep track of what files are being used for this stack
        master_reduction_files_used = podi_associations.collect_reduction_files_used(
            {}, {"calibrated": input_file} ) #ret['corrected_file']})

        if (swarp_params['preswarp_only']):
            logger.info("No prepping, working towards pre-swarp only")
            ret['corrected_file'] = input_file
            ret['gain'] = 1.0
            ret['skylevel'] = 0.0
            ret['weight'] = 1.
            ret["master_reduction_files"] = master_reduction_files_used
            hdulist.close()
            output_queue.put(ret)
            input_queue.task_done()
            continue

        #
        # Compute on how to scale the flux values
        #
        fluxscale_value = numpy.NaN
        magzero = hdulist[0].header['PHOTZP_X'] if 'PHOTZP_X' in hdulist[0].header else -99.
        if (magzero > 0 and not swarp_params['no-fluxscale']):
            fluxscale_value = math.pow(10, 0.4*(swarp_params['target_magzero']-magzero))
        else:
            exptime = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header else (
                hdulist[0].header['EXPTIME'] if 'EXPTIME' in hdulist[0].header else 1.0)
            fluxscale_value = 1./exptime

        # Assemble the temporary filename for the corrected frame
        suffix = None
        # Now construct the output filename
        if (suffix is None and swarp_params['no-fluxscale'] and not numpy.isnan(fluxscale_value)):
            suffix = "exptimenorm"
        if (suffix is None and not numpy.isnan(swarp_params['target_magzero']) and not numpy.isnan(fluxscale_value)):
            suffix = "fluxnorm"
        if (suffix is None and swarp_params['use_nonsidereal']):
            suffix = "nonsidereal"
        if (suffix is None and swarp_params['use_ephemerides']):
            suffix = "ephemerides"
        if (suffix is None and options['skip_otas'] != []):
            suffix = "otaselect"
        if (suffix is None and options['illumcorr']):
            suffix = "illumcorr"
        if (suffix is None and options['bpm_dir'] is not None):
            suffix = "bpmfixed"
        if (suffix is None and not swarp_params['subtract_back'] == 'swarp'):
            suffix = "skysub"
        if (suffix is None and options['photflat'] is not None):
            suffix = "photflat"

        if (suffix is not None):
            corrected_filename = "%(single_dir)s/%(obsid)s.%(suffix)s.%(fileid)d.fits" % {
                "single_dir": swarp_params['unique_singledir'],
                "obsid": hdulist[0].header['OBSID'],
                "suffix": suffix,
                "fileid": fileid,
            }

        logger.info("Applying preparations ...")

        #
        # Check if we need to apply any corrections
        #
        if (corrected_filename is None or
            (os.path.isfile(corrected_filename) and swarp_params['reuse_singles'])
        ):
            # Either we don't need to apply any corrections or we can re-use an 
            # older file with these corrections already applied
            logger.info("No correction needs to be applied!")
            ret['corrected_file'] = input_file
        else:
            gain = hdulist[0].header['GAIN']

            if (swarp_params['use_nonsidereal']):
                logger.debug("Applying the non-sidereal motion correction")

                # First apply the non-sidereal correction to all input frames
                # Keep frames that are already modified
                from podi_collectcells import apply_nonsidereal_correction
                # print options['nonsidereal']
                nonsidereal_dradec = apply_nonsidereal_correction(hdulist, options, logger)
                ret["nonsidereal-dradec"] = nonsidereal_dradec

                try:
                    if (os.path.isfile(options['nonsidereal']['ref'])):
                        master_reduction_files_used = podi_associations.collect_reduction_files_used(
                            master_reduction_files_used, 
                            {"nonsidereal-reference": options['nonsidereal']['ref']})
                except:
                    pass

            if (swarp_params['use_ephemerides']):
                # get MJD of current frame
                mjd_thisframe = hdulist[0].header['MJD-OBS'] + 0.5*hdulist[0].header['EXPTIME']/2./86400
                mjd_ref = swarp_params['ephemerides']['ref-mjd']

                # print "\n"*10
                # print "ref-mjd",mjd_ref
                # print "this mjd",mjd_thisframe
                # print "min mjd", numpy.min(swarp_params['ephemerides']['data'][:,0])
                # print "max mjd", numpy.max(swarp_params['ephemerides']['data'][:,0])
                # print "\n"*10

                # now compute the Ra/Dec of the target in both the reference 
                # frame and in this frame
                ephem_data = swarp_params['ephemerides']['data']
                logger.debug("EPHEM_DATA: %s" % (ephem_data))
                ra_from_mjd = scipy.interpolate.interp1d( ephem_data[:,0], ephem_data[:,1], kind='linear' )
                dec_from_mjd = scipy.interpolate.interp1d( ephem_data[:,0], ephem_data[:,2], kind='linear' )

                ra_ref = ra_from_mjd(mjd_ref)
                ra_this = ra_from_mjd(mjd_thisframe)

                dec_ref = dec_from_mjd(mjd_ref)
                dec_this = dec_from_mjd(mjd_thisframe)

                # print "\n"*5, "ra//dec = ", ra_ref, dec_ref, "\n"*5
                # The Ra/Dec correction is how much the object has moved 
                # (as derived from the ephemerides) between the reference MJD 
                # and the timestamp of this frame
                d_ra = ra_ref - ra_this
                d_dec = dec_ref - dec_this
                # print d_ra, d_dec
                d_days = mjd_thisframe - mjd_ref
                logger.debug("Applying ephemerid correction (%+.6f deg, %+.6f deg, dT=%.3f days)" % (
                    d_ra, d_dec, d_days))

                # Now apply these corrections to all extensions with an
                # apparently valid WCS system
                orig_ra, orig_dec = None, None
                for ext in hdulist:
                    if ('CRVAL1' in ext.header and
                        'CRVAL2' in ext.header):
                        # print ext.header['CRVAL1'], ext.header['CRVAL2'], (ext.header['EXTNAME'] if 'EXTNAME' in ext.header else "??")
                        orig_ra = ext.header['CRVAL1'] if orig_ra is None else orig_ra
                        orig_dec = ext.header['CRVAL2'] if orig_dec is None else orig_dec
                        ext.header['CRVAL1'] += d_ra
                        ext.header['CRVAL2'] += d_dec
                        # print ext.header['CRVAL1'], ext.header['CRVAL2'], (ext.header['EXTNAME'] if 'EXTNAME' in ext.header else "??")

                logger.debug("Pre-correction Ra/Dec was: %12.7f  %+12.7f" % (orig_ra, orig_dec))
                logger.debug("Post-corrected Ra/Dec is: %12.7f %+12.7f" % (orig_ra + d_ra, orig_dec + d_dec))
                
                ret["nonsidereal-dradec"] = numpy.array([numpy.NaN, -d_dec, -d_ra])


            if (options['skip_otas'] != []):
                logger.debug("Skipping some OTAs")
                ota_list = []
                for ext in hdulist:
                    ota = -1
                    try:
                        ota = int(ext.header['EXTNAME'][3:5])
                    except:
                        pass
                    if (ota in options['skip_otas']):
                        logger.debug("skipping ota %s as requested" % (ext.header['EXTNAME']))
                        continue
                    ota_list.append(ext)

                # Save the modified OTA list for later
                hdulist = pyfits.HDUList(ota_list)

            #
            # Delete the dead OTA 2,1 for data taken past Feb 1st, 2019 (MJD ~ 58515)
            # TODO: Add command line option to override this function
            #
            if (hdulist[0].header['MJD-OBS'] > 58515):
                # this data has a dead OTA 2,1
                ota_list = []
                for ext in hdulist:
                    if (ext.name == "OTA21.SCI"):
                        continue
                    ota_list.append(ext)
                hdulist = pyfits.HDUList(ota_list)

            if (options['bpm_dir'] is not None):
                logger.debug("Applying bad-pixel masks")
                for ext in range(len(hdulist)):
                    if (not is_image_extension(hdulist[ext])):
                        continue

                    fppos = None
                    if ('FPPOS' in hdulist[ext].header):
                        fppos = hdulist[ext].header['FPPOS']
                    if (fppos is not None):
                        region_file = "%s/bpm_%s.reg" % (options['bpm_dir'], fppos)
                        if (os.path.isfile(region_file)):
                            mask_broken_regions(hdulist[ext].data, region_file)
                            master_reduction_files_used = podi_associations.collect_reduction_files_used(
                                master_reduction_files_used, {"bpm": region_file})

        
            # Loop over all extensions and only select those that are not marked as guide chips
            if (True): #options['skip_otas'] != []):
                logger.debug("Sorting out guide-OTAs")
                ota_list = []
                for ext in hdulist:
                    if ('CELLMODE' in ext.header and
                        ext.header['CELLMODE'].find("V") >= 0):
                        logger.debug("skipping ota %s as requested" % (ext.header['EXTNAME']))
                        continue
                    if (swarp_params['skip_guide_otas']):
                        if (is_guide_ota(hdulist[0], ext)):
                            logger.info("OTA %s is likely guide-OTA, skipping" % (ext.name))
                            continue

                    ota_list.append(ext)

                # Save the modified OTA list for later
                hdulist = pyfits.HDUList(ota_list)

            # Mask out saturated pixels
            if (swarp_params['mask_saturated'] is not None):
                for ext in hdulist:
                    if (not is_image_extension(ext)):
                        continue
                    ext.data[ext.data > swarp_params['mask_saturated']] = numpy.NaN

            if (options['illumcorr_dir'] is not None):
                illum_file = podi_illumcorr.get_illumination_filename(
                    options['illumcorr_dir'], hdulist[0].header['FILTER'], hdulist[0].header['BINNING'])
                logger.debug("Applying illumination correction (%s)" % (illum_file))
                master_reduction_files_used = podi_associations.collect_reduction_files_used(
                    master_reduction_files_used, {"illumination": illum_file})
                podi_illumcorr.apply_illumination_correction(hdulist, illum_file)

            if (swarp_params['wipe_cells'] is not None):
                binning = hdulist[0].header['BINNING']
                # XXX ADD SUPPORT FOR SOFTWARE BINNING AS WELL HERE !!!
                logger.info("Wiping out specified cells")
                for ext in hdulist:
                    if (not is_image_extension(ext)):
                        continue
                    wipecells(ext, swarp_params['wipe_cells'], binning=binning)

            # clear out the bottom few rows in each OTA cell to account for the fat-zero problem
            if (swarp_params['cheese_grate'] is not None):
                if (swarp_params['cheese_grate'] == 'auto'):
                    # TODO: change the equation for the auto scaling
                    _skylevel = hdulist[0].header['SKYLEVEL']
                    n_lines = int((120 - _skylevel) / 10 * 1.5)
                    if (n_lines <= 0):
                        n_lines = 0
                    logger.info("Apply de-cheese-grating [auto-mode]: %d lines", n_lines)
                else:
                    n_lines = swarp_params['cheese_grate']
                    logger.info("Apply de-cheese-grating [manual mode]: %d lines", n_lines)
                fpl = FocalPlaneLayout(hdulist)
                for ext in hdulist:
                        if (not is_image_extension(ext)):
                            continue
                        if (fpl.get_detector_generation(ext.header) >= 7):
                            # this is a lot-7 detector, nothing to do
                            logger.info("Skipping de-cheese-grate for lot-7 OTA %s" % (ext.name))
                            continue

                        for _x in range(8):
                            for _y in range(8):
                                x1,x2,y1,y2 = cell2ota__get_target_region(_x, _y, trimcell=0)
                                ext.data[y1:y1+n_lines+1, x1:x2+1] = numpy.NaN

            skylevel = 0.
            if (not swarp_params['subtract_back'] == False and
                not swarp_params['subtract_back'] in ['swarp', "_REGIONS_"] and
                not swarp_params['subtract_back'].startswith("IC,")):
                skylevel = numpy.NaN
                if (swarp_params['subtract_back'] in hdulist[0].header):
                    skylevel = hdulist[0].header[swarp_params['subtract_back']]
                else:
                    try:
                        skylevel = float(swarp_params['subtract_back'])
                    except ValueError:
                        logger.warning("Could not determine sky-level (%s), skipping sky-subtraction" % (
                            swarp_params['subtract_back']))
                    except:
                        raise
                if (not numpy.isnan(skylevel)):
                    for ext in hdulist:
                        if (not is_image_extension(ext)):
                            continue
                        ext.data -= skylevel
                        logger.debug("Subtracting skylevel (%f) from extension %s" % (skylevel, ext.name))
            elif (swarp_params['subtract_back'].startswith("IC,")):
                logger.info("Using scale & subtract model for BG subtraction")
                ic_file = swarp_params['subtract_back'][3:]
                if (not os.path.isfile(ic_file)):
                    logger.critical("Unable to open IC file for sky-subtraction: %s" % (ic_file))
                    skylevel = 0
                else:
                    # Now we do have an illumination correction template for sky-subtraction
                    hdulist, skyframe = ic_background.scale_subtract_background_model(
                        in_filename=hdulist,
                        ic_file=ic_file,
                        out_filename=None, per_ota=True,
                        twod_model=True,
                        logger=logger,
                    )
                    skylevel = hdulist[0].header['SKYLEVEL']
                    skyframe.writeto(corrected_filename[:-5]+".skybg.fits", overwrite=True)
            elif (swarp_params['subtract_back'] == 'swarp'):
                skylevel = hdulist[0].header['SKYLEVEL']
            elif (swarp_params['subtract_back'] == "_REGIONS_"):
                #
                # Get measurements for each of the sky-regions in the current OTA
                #
                all_sky_samples = None
                for ext in hdulist:
                    if (not is_image_extension(ext)):
                        continue
                    sky_samples = sample_background_using_ds9_regions(ext, swarp_params['sky_regions'])
                    if (type(sky_samples) == type(None)):
                        continue
                    #print sky_samples.shape
                    #print all_sky_samples.shape if all_sky_samples != None else "XXX"
                    all_sky_samples = sky_samples if type(all_sky_samples) == type(None) else \
                                      numpy.append(all_sky_samples, sky_samples, axis=0)
                try:
                    numpy.savetxt(os.path.basename(input_file)+".sky", all_sky_samples)
                    logger.info("Found %d sky-samples in user-defined regions" % (all_sky_samples.shape[0]))
                except:
                    pass
                #
                # Now we have all sky-samples across the entire focal plane
                #
                if (all_sky_samples is not None and
                    all_sky_samples.size > 0):
                    cleaned = three_sigma_clip(all_sky_samples[:,2])
                    skylevel = numpy.median(cleaned)
                    skynoise = numpy.std(cleaned)
                    logger.info("custom: %f vs %f (+/- %f)" % (skylevel, hdulist[0].header['SKYLEVEL'], skynoise))
                else:
                    skylevel = hdulist[0].header['SKYLEVEL']
                    logger.warning("Not enough sky-samples from ds9 regions, using default skylevel (%f)" % (
                        skylevel))
                for ext in hdulist:
                    if (not is_image_extension(ext)):
                        continue
                    ext.data -= skylevel
                    logger.debug("Subtracting ds9/user skylevel (%f) from extension %s" % (skylevel, ext.name))
            else:
                skylevel = hdulist[0].header['SKYLEVEL']

            if (options['photflat'] is not None):
                photflat_filename = options['photflat']
                if (os.path.isfile(photflat_filename)):
                    # apply correction
                    photflat_hdu = pyfits.open(photflat_filename)
                    master_reduction_files_used = podi_associations.collect_reduction_files_used(
                        master_reduction_files_used, {"photflat": photflat_filename})
                    for ext in hdulist:
                        try:
                            photflat_ext = photflat_hdu[ext.name]
                            ext.data /= photflat_ext.data
                            logger.debug("Successfully applied photometric flatfield for %s" % (ext.name))
                        except (KeyError, ValueError, TypeError):
                            continue

            if (swarp_params['autophotflat'] and apf_data is not None):
                # determine which frames to use for the auto-photflat
                this_mjd = hdulist[0].header['MJD-OBS']
                use_apf_data = apf_data.copy()
                use_apf_mjds = numpy.array([float(f) for f in use_apf_data[:,0]])
                logger.info("MJD for %s is %f" % (input_file, this_mjd))
                this_index = numpy.arange(use_apf_data.shape[0])[(use_apf_data[:,1] == input_file)]
                # print this_index
                # print "timeframe:", swarp_params['apf_timeframe']
                # print "nframes:", swarp_params['apf_nframes']

                if (swarp_params['apf_timeframe'] is not None):
                    _mjd = use_apf_mjds
                    # print _mjd
                    select = numpy.fabs(_mjd-this_mjd) <= swarp_params['apf_timeframe']
                    # print select
                    logger.info("Found %d frames within %.2f days of mjd=%.6f" % (
                        numpy.sum(select), swarp_params['apf_timeframe'], this_mjd
                    ))
                elif (swarp_params['apf_nframes'] is not None):
                    idx = numpy.arange(use_apf_data.shape[0])
                    select = numpy.fabs(idx-this_index) <= swarp_params['apf_nframes']
                else:
                    select = numpy.ones((use_apf_data.shape[0]), dtype=numpy.bool)
                if (swarp_params['apf_noauto']):
                    select[this_index] = False
                files_for_apf = [str(f) for f in use_apf_data[select, 1]]
                logger.info("Using these files to compute phot.flat for %s:\n%s" % (input_file, "\n".join(list(files_for_apf))))

                # Now that we have a list of files, create the photometric
                # flatfield
                custom_photflat, apf_extras = podi_photflat.create_photometric_flatfield(
                    filelist = files_for_apf,
                    smoothing=swarp_params['apf_resolution'],
                    strict_ota=False,
                    return_interpolator=True,
                    parallel=True,
                    n_processes=1,
                )
                # save the photflat
                output_bn = os.path.splitext(outputfile)[0]
                custom_photflat.writeto(
                    "custom_apf_for_%s_in_%s.fits" % (hdulist[0].header['OBSID'], output_bn),
                    overwrite=True
                )


                        # XXX
            if (options['illumcorr_dir'] is not None and
                swarp_params['un_illumcorr']):
                logger.info("Un-doing illumination correction (%s)" % (illum_file))
                master_reduction_files_used = podi_associations.collect_reduction_files_used(
                    master_reduction_files_used, {"un_illumination": illum_file})
                podi_illumcorr.apply_illumination_correction(hdulist, illum_file, invert=False)

            if (not numpy.isnan(fluxscale_value)):
                logger.debug("Applying flux-scaling (%.10e)" % (fluxscale_value))
                for ext in hdulist:
                    if (not is_image_extension(ext)):
                        continue
                    ext.data *= fluxscale_value
                    logger.debug("Applying flux-scaling (%.10e) to extension %s" % (fluxscale_value, ext.name))

                # Apply fluxscaling to GAIN and SKYLEVLE as well
                gain /= fluxscale_value
                skylevel *= fluxscale_value

            # Check if the corrected file already exists - if not create it
            #if (not os.path.isfile(corrected_filename)):
            logger.debug("Writing correctly prepared file--> %s" % (corrected_filename))

            clobberfile(corrected_filename)
            hdulist.writeto(corrected_filename, overwrite=True)

            # Now change the filename of the input list to reflect 
            # the corrected file
            ret['corrected_file'] = corrected_filename
            ret['gain'] = gain
            ret['skylevel'] = skylevel

        #
        # Now also create a relative weight map for this frame
        # scaling factor is the exposure time of each frame
        #
        weight_hdulist = [hdulist[0]] # copy the primary header
        for ext in hdulist[1:]:
            if (not is_image_extension(ext)):
                continue

            weight_data = numpy.ones(ext.data.shape, dtype=numpy.float32) #* 100.
            ret['weight'] = 1.
            if (not numpy.isnan(fluxscale_value)):
                weight_data /= fluxscale_value
                ret['weight'] = 1./fluxscale_value
            weight_data[numpy.isnan(ext.data)] = 0.

            weight_img = pyfits.ImageHDU(header=ext.header, data=weight_data) 
            weight_hdulist.append(weight_img)


        # convert extension list to proper HDUList ...
        weight_hdulist = pyfits.HDUList(weight_hdulist)

        # ... and write extension to file
        weight_filename = ret['corrected_file'][:-5]+".weight.fits"
        clobberfile(weight_filename)
        weight_hdulist.writeto(weight_filename, overwrite=True)
        logger.info("Wrote input weight map to %s" % (weight_filename))

        # Finally, close the input file
        hdulist.close()

        #
        # Now we have the filename of the file to be used for the swarp-input
        #
        ret["master_reduction_files"] = master_reduction_files_used

        logger.debug("Sending return value to master process")
        output_queue.put(ret)
        input_queue.task_done()

        # print input_file, "\n", ret["master_reduction_files"]

    # end of routine



def prepare_input(inputlist, swarp_params, options,
                  apf_data=None):

    logger = logging.getLogger("PrepFiles")

    #
    # initialize queues for commands and return-values
    #
    in_queue = multiprocessing.JoinableQueue()
    out_queue = multiprocessing.Queue()

    #
    # fill queue with files to be processed
    #
    n_jobs = 0
    existing_inputlist = []
    for fn in inputlist: #i in range(len(inputlist)):
        if (not os.path.isfile(fn)):
            continue
        try:
            hdulist = pyfits.open(fn)
        except IOError:
            logger.error("Can't open file %s" % (fn))
            continue
        
        hdulist.close()
        existing_inputlist.append(fn)

    if (len(existing_inputlist) <= 0):
        logger.error("No valid files found")
        return

    wcs_inputlist = []
    photcal_inputlist = []
    for idx, fn in enumerate(existing_inputlist):
        hdulist = pyfits.open(fn)
        #
        # Perform some checks to only include valid frames in the stack
        #
        # Frame needs to have valid WCS solution
        if ('WCSCAL' in hdulist[0].header and
            not hdulist[0].header['WCSCAL'] and 
            not swarp_params['ignore_quality_checks']):
            logger.info("Excluding frame (%s) due to faulty WCS calibration" % (fn))
            #good_inputlist[idx] = None
            continue

        wcs_inputlist.append(fn)

        # and proper photometric calibration
        if ('MAGZERO' in hdulist[0].header and
            hdulist[0].header['MAGZERO'] <= 0 and
            not swarp_params['no-fluxscale'] and
            not swarp_params['ignore_quality_checks']):
            logger.info("Excluding frame (%s) due to missing photometric calibration" % (fn))
            #good_inputlist[idx] = None
            continue

        photcal_inputlist.append(fn)

    if (len(photcal_inputlist) > 0):
        logger.debug("Restricting input file list to files with valid photoemtric calibration")
        inputlist = photcal_inputlist
    elif (len(wcs_inputlist) > 0):
        logger.warning("No files with photometric calibration found, reverting to list of WCS-calibrated files")
        inputlist = wcs_inputlist
    else:
        logger.warning("No files with WCS and/or photometry found, reverting to unfiltered inputlist")
        inputlist = existing_inputlist

    for idx, fn in enumerate(inputlist):
        in_queue.put((fn, idx+1))
        n_jobs += 1

    logger.info("Queued %d jobs ..." % (n_jobs))

    #
    # Start worker processes
    #
    worker_args = (in_queue, out_queue, swarp_params, options, apf_data)
    processes = []
    for i in range(sitesetup.number_cpus):
        p = multiprocessing.Process(target=mp_prepareinput, args=worker_args)
        p.start()
        processes.append(p)

        # also add a quit-command for each process
        in_queue.put(None)
        
    #
    # wait until all work is done
    #
    in_queue.join()

    #
    # return the list of corrected files.
    #

    # Keep track of what files are being used for this stack
    master_reduction_files_used = {}
    corrected_file_list = []
    nonsidereal_offsets = []

    stack_start_time = 1e9
    stack_end_time = -1e9
    stack_total_exptime = 0
    stack_framecount = 0

    gain_list = numpy.zeros((n_jobs))
    skylevel_list = numpy.zeros((n_jobs))
    weight_list = numpy.zeros((n_jobs))

    for i in range(n_jobs):
        ret = out_queue.get()
        logger.debug("Received results from job %d" % (i+1))

        gain_list[i] = ret['gain']
        skylevel_list[i] = ret['skylevel']
        weight_list[i] = ret['weight']

        # Also set some global stack-related parameters that we will add to the 
        # final stack at the end
        # mjd_obs_start = hdulist[0].header['MJD-OBS']
        # exptime = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header else \
        #           hdulist[0].header['EXPTIME']
        # mjd_obs_end = hdulist[0].header['MJD-OBS'] + (exptime/86400.)

        # logger.debug("Exposure time: %f, MJD=%f" % (exptime, mjd_obs_start))

        stack_total_exptime += ret['exptime']
        stack_framecount += 1

        stack_start_time = numpy.min([stack_start_time, ret['mjd_obs_start']])
        stack_end_time = numpy.max([stack_end_time, ret['mjd_obs_end']])

        master_reduction_files_used = podi_associations.collect_reduction_files_used(master_reduction_files_used, 
                                                                   ret['master_reduction_files'])

        corrected_file_list.append(ret['corrected_file'])
        nonsidereal_offsets.append(ret['nonsidereal-dradec'])

    photom_list = (gain_list, skylevel_list, weight_list)

    #
    # By now all frames have all corrections applied,
    # so we can go ahead and stack them as usual
    #
    
    # Make sure to join/terminate all processes
    for p in processes:
        p.join()
        
    logger.info("All files prepared!")

    return (corrected_file_list, 
            stack_total_exptime, 
            stack_framecount, 
            stack_start_time, 
            stack_end_time, 
            master_reduction_files_used,
            nonsidereal_offsets,
            photom_list)



def cleanup_singles(unique_singledir, logger):

    try:
        shutil.rmtree(unique_singledir)
    except:
        logger.error("There was a problem with recursively deleting the temp directory")
        podi_logging.log_exception()
        pass

    return



def mp_swarp_single(sgl_queue, dum):

    while(True):
        cmd = sgl_queue.get()
        if (cmd is None):
            sgl_queue.task_done()
            break

        swarp_cmd, prepared_file, single_file, swarp_params, nonsidereal_dradec, create_mask = cmd
        logger = logging.getLogger("MPSwarpSgl(%s)" % (os.path.basename(single_file)))

        hdulist = pyfits.open(prepared_file)
        obsid = hdulist[0].header['OBSID']
        logger.info("Starting work on focal-plane frame ...")

        single_created_ok = False
        try:
            logger.debug(" ".join(swarp_cmd.split()))
            ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
            (swarp_stdout, swarp_stderr) = ret.communicate()
            if (sitesetup.log_shell_output):
                logger.debug("\nCommand:\n%s\n--> Returncode: %d\n---\nStd.Out:\n%s\n---\nStd.Err:\n%s\n---", 
                             swarp_cmd, ret.returncode, swarp_stdout, swarp_stderr)
            # logger.debug("swarp stdout:\n"+swarp_stdout)
            # if (len(swarp_stderr) > 0 and ret.returncode != 0):
            #     logger.warning("swarp stderr:\n"+swarp_stderr)
            # elif (sitesetup.log_shell_output):
            #     logger.debug("swarp stderr:\n"+swarp_stderr)

            # Add some basic headers from input file to the single file
            # this is important for the differencing etc.
            logger.info("adding header to single file (%s)" % (single_file))
            hdu_single = pyfits.open(single_file,  mode='update')
            for hdrkey in [
                    'TARGRA', 'TARGDEC',
                    'FILTER', 'FILTERID', 'FILTDSCR', 
                    'OBSID', 'OBJECT', 
                    'EXPTIME',
                    'DATE-OBS', 'TIME-OBS', 'MJD-OBS']:
                if (hdrkey in hdulist[0].header):
                    key, val, com = hdulist[0].header.cards[hdrkey]
                    hdu_single[0].header[key] = (val, com)
            # hdu_single.writeto(single_file, clobber=True)
            hdu_single.flush()
            hdu_single.close()

            single_created_ok = True
            #print "\n".join(swarp_stderr)
            # single_prepared_files.append(single_file)
        except OSError as e:
            podi_logging.log_exception()
            print >>sys.stderr, "Execution failed:", e

        #
        # Apply the mask if requested
        #
        if (single_created_ok and create_mask and swarp_params['mask-list']):
            # not swarp_params['mask'] is None and

            #
            # Based on the MJD of this frame, determine the best mask-frame
            #
            try:
                frame_mjd = hdulist[0].header['MJD-OBS']
            except:
                frame_mjd = 0
            diff_mjd = numpy.fabs(swarp_params['mask-mjds'] - frame_mjd)
            # Find the frame with the smallest delta_mjd
            idx = numpy.argmin(diff_mjd)
            mask_file = swarp_params['mask-list'][idx]
            logger.info("Using mask %s for frame %s" % (mask_file, prepared_file))
 
            #
            # Apply the non-sidereal correction to the mask
            # (only in non-sidereal mode)
            #
            this_mask = single_file[:-5]+".maskraw.fits"
            logger.info("This mask: %s" % (this_mask))
            mask_hdu = pyfits.open(mask_file) #swarp_params['mask'])

            logger.info("Applying non-sid corr: %s" % (str(nonsidereal_dradec)))
            if (nonsidereal_dradec is not None):
                
                d_radec = numpy.array(nonsidereal_dradec)

                # Correct the declination
                mask_hdu[0].header['CRVAL2'] -= d_radec[1]
                # correct RA, compensating for cos(declination)

                if (d_radec.shape[0] == 3):
                    mask_hdu[0].header['CRVAL1'] -= d_radec[2]
                else:
                    cos_dec = math.cos(math.radians(mask_hdu[0].header['CRVAL2']))
                    mask_hdu[0].header['CRVAL1'] -= d_radec[0] / cos_dec

            clobberfile(this_mask)
            mask_hdu.writeto(this_mask, overwrite=True)
            mask_hdu.close()
            logger.info("wrote raw mask with fudged WCS: %s" % (this_mask))

            #
            # Swarp the mask to the identical pixelgrid as the single frame
            #
            hdu_single = pyfits.open(single_file,  mode='update')
            out_crval1 = hdu_single[0].header['CRVAL1']
            out_crval2 = hdu_single[0].header['CRVAL2']
            out_naxis1 = hdu_single[0].header['NAXIS1']
            out_naxis2 = hdu_single[0].header['NAXIS2']
            hdu_single.close()

            mask_aligned = single_file[:-5]+".mask.fits"
            swarp_mask = """
                %(swarp)s -c %(swarp_default)s 
                -IMAGEOUT_NAME %(mask_aligned)s
                -WEIGHTOUT_NAME /dev/null
                -CENTER_TYPE MANUAL
                -CENTER %(center_ra)f,%(center_dec)f
                -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d
                -RESAMPLE_DIR %(resample_dir)s
                -RESAMPLING_TYPE BILINEAR
                -PIXEL_SCALE %(pixelscale)f \
                -PIXELSCALE_TYPE %(pixelscale_type)s \
                -COMBINE Y \
                -COMBINE_TYPE AVERAGE \
                -SUBTRACT_BACK N \
                %(mask_raw)s
            """ % {
                'swarp': sitesetup.swarp_exec,
                'swarp_default': "%s/config/swarp.default" % (sitesetup.exec_dir),
                'mask_raw': this_mask,
                'mask_aligned': mask_aligned,
                'center_ra': out_crval1,
                'center_dec': out_crval2,
                'imgsizex': out_naxis1,
                'imgsizey': out_naxis2,
                'resample_dir': swarp_params['unique_singledir'],
                'pixelscale': swarp_params['pixelscale'],
                'pixelscale_type': "MANUAL",
                }

            # print "\n"*3," ".join(swarp_mask.split()),"\n"*3
            logger.info("Matching global mask to frame ...")
            try:
                start_time = time.time()
                logger.debug("Matching global mask to frame:\n"+" ".join(swarp_mask.split()))
                ret = subprocess.Popen(swarp_mask.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
                (swarp_stdout, swarp_stderr) = ret.communicate()
                if (sitesetup.log_shell_output):
                    logger.debug("\nCommand:\n%s\n--> Returncode: %d\n---\nStd.Out:\n%s\n---\nStd.Err:\n%s\n---", 
                                 swarp_cmd, ret.returncode, swarp_stdout, swarp_stderr)
                # if (len(swarp_stderr) > 0 and ret.returncode != 0):
                #     logger.warning("swarp stderr:\n"+swarp_stderr)
                # else:
                #     logger.debug("swarp stderr:\n"+swarp_stderr)
                end_time = time.time()
                logger.debug("Creating mask for X finished successfully after %.2d seconds" % (end_time-start_time))
            except:
                pass

            logger.info("Done with aligned mask: %s" % (mask_aligned))

            #
            # Multiply the weight mask of this single frame with the mask, 
            # thus eliminating all sources we want to mask out
            #
            logger.info("Applying mask to weightmap (%s)" % (single_file))
            mask_hdu = pyfits.open(mask_aligned)
            weightmap_file = single_file[:-5]+".weight.fits"
            weightmap_hdu = pyfits.open(weightmap_file, mode='update')

            weightmap_hdu[0].data[mask_hdu[0].data > 0] = 0.

            weightmap_hdu.flush()
            weightmap_hdu.close()

        hdulist.close()

        if (swarp_params['aggressive_clean']):
            clobberfile(prepared_file)
            prepared_file_weight = prepared_file[:-5]+".weight.fits"
            clobberfile(prepared_file_weight)

            # XXX ALSO DELETE INTERMEDIATE MASK FILES 

        sgl_queue.task_done()


def create_mask(fitsfile, swarp_params):

    logger.info("Creating mask (from %s)" % (fitsfile))

    unique_singledir = swarp_params['unique_singledir']
    #
    # Extract only the first extension
    # 
    hdulist = pyfits.open(fitsfile)

    # Get the mid-stack MJD date to find the closest mask in time when 
    # multiple masks are given
    mjd = hdulist[0].header['MJD-OBS']
    if ('MJD-MID' in hdulist[0].header): mjd = hdulist[0].header['MJD-MID']

    _, fitsbase = os.path.split(os.path.abspath(fitsfile))
    image_only_fits = "%s/%s.primaryonly.fits" % (unique_singledir, fitsbase)
    pyfits.HDUList([hdulist[0]]).writeto(image_only_fits, overwrite=True)

    #
    # We use source-extractor to create the mask
    #        
    sex_default = "%s/config/swarpstack_mask.conf" % (sitesetup.exec_dir)
    params_default = "%s/config/swarpstack_mask.params" % (sitesetup.exec_dir)
    segmentation_file = "%s/%s.segmentation.fits" % (unique_singledir, fitsbase)

    sex_cmd = """
    %(sex)s 
    -c %(config)s 
    -PARAMETERS_NAME %(params)s
    -DETECT_THRESH %(nsigma)f 
    -DETECT_MINAREA %(minarea)f
    -CHECKIMAGE_TYPE SEGMENTATION
    -CHECKIMAGE_NAME %(segmentation_file)s
    %(inputfits)s
    """ % {
        'sex': sitesetup.sextractor,
        'config': sex_default,
        'params': params_default,
        'nsigma': swarp_params['mask-nsigma'],
        'minarea': swarp_params['mask-npix'],
        'segmentation_file': segmentation_file,
        'inputfits': image_only_fits,
        #swarp_params['mask-fits'],
        }

    # print "\n"*5,sex_cmd,"\n"*5
    # print " ".join(sex_cmd.split())

    #
    # Run SourceExtractor to compute the segmentation map
    #
    try:
        start_time = time.time()
        logger.debug("Computing segmentation map:\n%s" % (" ".join(sex_cmd.split())))
        
        ret = subprocess.Popen(sex_cmd.split(), 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE) #, shell=True)
        (sex_stdout, sex_stderr) = ret.communicate()
        if (sitesetup.log_shell_output):
            logger.debug("\nCommand:\n%s\n--> Returncode: %d\n---\nStd.Out:\n%s\n---\nStd.Err:\n%s\n---", 
                         sex_cmd, ret.returncode, sex_stdout, sex_stderr)

        # logger.debug("sex stdout:\n"+sex_stdout)
        # if (len(sex_stderr) > 0 and ret.returncode != 0):
        #     logger.warning("sex stderr:\n"+sex_stderr)
        # else:
        #     logger.debug("sex stderr:\n"+sex_stderr)
        end_time = time.time()
        logger.info("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
    except OSError as e:
        podi_logging.log_exception()
        print >>sys.stderr, "Execution failed:", e

    #
    # Now convert the segmentation mask into a object mask
    # 
    seghdu = pyfits.open(segmentation_file)
    weight = numpy.zeros(shape=seghdu[0].data.shape, dtype=numpy.float32)
    weight[seghdu[0].data > 0] = 1.

    maskhdu = pyfits.HDUList([pyfits.PrimaryHDU(data=weight, header=seghdu[0].header.copy())])
    maskfile = "%s/%s.mask.fits" % (unique_singledir, fitsbase)
    clobberfile(maskfile)
    maskhdu.writeto(maskfile, overwrite=True)

    seghdu.close()
    maskhdu.close()
    del seghdu
    del maskhdu

    return maskfile, mjd


def swarpstack(outputfile, 
               inputlist, 
               swarp_params, 
               options, 
               keep_intermediates=False, 
               unique_dir=None,
               ):

    logger = logging.getLogger("SwarpStack - Prepare")

    # Figure out the config path
    swarp_default = "%s/config/swarp.default" % (sitesetup.exec_dir)
    logger.debug("Using swarp-default in %s" % (swarp_default))

    if (len(inputlist) <= 0):
        logger.error("No (valid) input files specified!")
        return

    master_reduction_files_used = {}
    
    ############################################################################
    #
    # Construct a unique name to hold intermediate files
    #
    ############################################################################
    if (unique_dir is not None and os.path.isdir(unique_dir)):
        unique_singledir = unique_dir
        keep_intermediates = True
    else:
        process_id = os.getpid()
        hostname = socket.gethostname()
        # Create a temporary directory
        unique_singledir = tempfile.mkdtemp(dir=sitesetup.swarp_singledir,
                                            prefix="%s-%05d----" % (hostname, process_id))
    # Save the directory name in swarp_params
    logger.info("Storing intermediate files in %s ..." % (unique_singledir))
    swarp_params['unique_singledir'] = unique_singledir


    ############################################################################
    #
    # Prepare the mask file(s) if requested
    #
    ############################################################################
    if (swarp_params['mask-fits'] is not None and
        swarp_params['mask'] is None):

        mask_list = []
        mask_mjds = []
        for mf in swarp_params['mask-fits']:
            maskfile, mjd = create_mask(mf, swarp_params)
            mask_list.append(maskfile)
            mask_mjds.append(mjd)

        swarp_params['mask-list'] = mask_list
        swarp_params['mask-mjds'] = numpy.array(mask_mjds)
            
        # swarp_params['mask'] = maskfile

    logger.info("Removing some OTAs from input: %s" % (str(options['skip_otas'])))

    logger.info("Using de-cheese-grate parameter: %s" % (str(swarp_params['cheese_grate'])))

    ############################################################################
    # Prepare the meta-data we need for the auto-phot.-flatfielding
    ############################################################################
    apf_data = None
    if (swarp_params['autophotflat']):
        # User asked for it, read all required data, including photometric
        # calibration data and time-stamps
        print("gathering all data for auto-phot.flat")
        apf_data = []
        for fn in inputlist:
            if (not os.path.isfile(fn)):
                continue
            print(fn)
            hdulist = pyfits.open(fn)
            mjd_obs = hdulist[0].header['MJD-OBS']
            apf_data.append([mjd_obs, fn])
        apf_data = numpy.array(apf_data)

        mjd_sort = numpy.argsort(apf_data[:,0])
        apf_data = apf_data[mjd_sort]

        print(apf_data)
        logger.debug("done with data gathering")
        pass

    ############################################################################
    #
    # Generate the illumination correction file if this is requested
    #
    ############################################################################
    if (swarp_params['illumcorr'] == "autogenerate" and
        swarp_params['illumcorrfiles'] != None):
        
        # Make sure all files exist
        ic_filelist = []
        for fn in swarp_params['illumcorrfiles'].split(","):
            if (os.path.isfile(fn)):
                ic_filelist.append(fn)
        
        logger.debug("Using the following files to create an illumination correction:\n -- %s" % (
            "\n -- ".join(ic_filelist)))

        ic_filename = "%s.illumcorr.fits" % (outputfile[:-5])
        podi_illumcorr.prepare_illumination_correction(
            filelist=ic_filelist,
            outfile=ic_filename,
            tmpdir=unique_singledir,
            redo=True)

        # Change the internal parameter to use the newly generated
        # illumination correction file during stacking
        params['illumcorr'] = ic_filename
        options['illumcorr_dir'] = ic_filename

    ###########################################################################
    #
    # Handle the optional user-requested zeropoint determination based on the 
    # -best or median options
    #
    ###########################################################################
    logger.info("Checking flux-scaling factors")
    all_zp = numpy.empty((len(inputlist)))
    all_zp[:] = numpy.NaN
    all_saturation = numpy.zeros_like(all_zp)
    for idx, fn in enumerate(inputlist):
        if (os.path.isfile(fn)):
            _hdu = pyfits.open(fn)
            magzero = _hdu[0].header['PHOTZP_X'] if 'PHOTZP_X' in _hdu[0].header else numpy.NaN
            logger.debug("ZP (%s) = %.4f" % (fn, magzero))
            all_zp[idx] = magzero

            # read saturation level from file, or default to 58K
            all_saturation[idx] = \
                _hdu[0].header['SATURATE'] if 'SATURATE' in _hdu[0].header \
                    else 58000.

    logger.debug(str(all_zp))
    all_zp = all_zp[numpy.isfinite(all_zp)]
    all_saturation = all_saturation[numpy.isfinite(all_zp)]

    if (swarp_params['target_magzero'] in ['best', 'median']):
        # Load all frames, and get a list of all available magzero headers
        if (all_zp.shape[0] == 0):
            final_magzero = 25.
            logger.warning("Didn't find any valid zeropoints")
        elif (swarp_params['target_magzero'] == "best"):
            final_magzero = numpy.max(all_zp)
            logger.debug("Selecting ''best'' phot. ZP of %.3f" % (final_magzero))
        elif (swarp_params['target_magzero'] == "median"):
            final_magzero = numpy.median(all_zp)
            logger.debug("Selecting ''median'' phot. ZP of %.3f" % (final_magzero))
        else:
            logger.error("Problem determing which phot. ZP to use")
            final_magzero = 25.0
        swarp_params['target_magzero'] = final_magzero
    elif (type(swarp_params['target_magzero']) == str and
          os.path.isfile(swarp_params['target_magzero'])):
        _hdu = pyfits.open(swarp_params['target_magzero'])
        magzero = _hdu[0].header['PHOTZP_X'] if 'PHOTZP_X' in _hdu[0].header else 25.0
        logger.debug("Matching ZP to %s ==> ZP = %.4f" % (swarp_params['target_magzero'], magzero))
        swarp_params['target_magzero'] = magzero
    else:
        try:
            swarp_params['target_magzero'] = float(swarp_params['target_magzero'])
        except:
            logger.error("Invalid parameter or file not found: %s" % (str(swarp_params['target_magzero'])))
            swarp_params['target_magzero'] = 25.0
    logger.info("Scaling all frames to common phot. ZP of %.4f" % (swarp_params['target_magzero']))

    # re-calculate the flux-scaled saturation levels
    flux_scaling = numpy.power(10, 0.4*(swarp_params['target_magzero']-all_zp))
    all_saturation *= flux_scaling
    logger.info("Final saturation levels: %s" % (str(all_saturation)))

    ############################################################################
    #
    # Prepare all QR'ed input files, applying additional corrections where needed
    #
    ############################################################################

    stack_start_time = 1e9
    stack_end_time = -1e9
    stack_total_exptime = 0
    stack_framecount = 0

    # print "input=",inputlist
    # print "output=",outputfile

    modified_files, stack_total_exptime, stack_framecount, \
        stack_start_time, stack_end_time, master_reduction_files_used, \
        nonsidereal_offsets, photom_lists = \
        prepare_input(inputlist, swarp_params, options,
                      apf_data)

    # print modified_files
    inputlist = modified_files

    # print "\n\n".join(inputlist)

    #print photom_lists
    gain_list, skylevel_list, weight_list = photom_lists


    add_only = swarp_params['add'] and os.path.isfile(outputfile)
    if (add_only):
        logger.info("Activating ADD mode")

    if (outputfile.endswith(".fits")):
        outputfile = outputfile[:-5]

    header_only_file = "%s/preswarp.fits" % (unique_singledir)
    logger.debug("Using header-only-file: %s" % (header_only_file))

    # Make sure the reference file is a valid file
    if (swarp_params['reference_file'] is not None):
        if (os.path.isfile(swarp_params['reference_file'])):
            logger.info("Using %s as reference file" % (swarp_params['reference_file']))
        else:
            logger.error("Could not find specified reference file (%s)" % (swarp_params['reference_file']))
            swarp_params['reference_file'] = None

    logging.debug("Using modified input list: %s" % (str(inputlist)))

    #
    # Open the first frame in the list, and preserve its header.
    # This header will be (partially) merged into the final output
    #
    hdulist = pyfits.open(inputlist[0])
    reference_header = hdulist[0].header.copy()
    hdulist.close()

    ############################################################################
    #
    # Figure out the pixel-grid and sky-coverage of the final stack 
    #
    ############################################################################
    logger.info("reference_file = %s" % (str(swarp_params['reference_file'])))
    if (swarp_params['cutout_list'] is not None):

        if (swarp_params['pixelscale'] <= 0):
            swarp_params['pixelscale'] = 0.11

        out_crval1 = swarp_params['cutout_list']['ra']
        out_crval2 = swarp_params['cutout_list']['dec']

        #
        # compute cutout dimensions in pixels
        #
        imgsize_arcsec = swarp_params['cutout_list']['size']
        imgsize_pixels = int(math.ceil(imgsize_arcsec / swarp_params['pixelscale']))
        out_naxis1 = imgsize_pixels
        out_naxis2 = imgsize_pixels
        logger.info("Setting custom image cutout region: ra=%f dec=%f imgsize=%dpx (%f'')" % (
            out_crval1, out_crval2, imgsize_pixels, imgsize_arcsec))

    elif (add_only or swarp_params['reference_file'] is not None):
        #
        # This is the simpler add-only mode
        #

        logger.info("Reading sky-coverage from reference frame")

        if (swarp_params['reference_file'] is not None):
            output_info = pyfits.open(swarp_params['reference_file'])
        else:
            # Open the existing output header and get data from there
            output_info = pyfits.open(outputfile+".fits")

        logger.info("Stack information...")
        logger.info("   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header))
        logger.info("       Output center: %(CRVAL1)10.6f / %(CRVAL2)10.6f" % (output_info[0].header))
        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']

        output_info.close()

        if (swarp_params['pixelscale'] <= 0):
            swarp_params['pixelscale'] = math.fabs(output_info[0].header['CD1_1']) * 3600.
            logger.info("Computing pixelscale from data: %.4f arcsec/pixel" % (swarp_params['pixelscale']))
    else:
        #
        # This is the regular start-from-scratch mode
        #

        # Set some Swarp options
        swarp_opts = """ \
               -IMAGEOUT_NAME %(imageout)s \
               -WEIGHTOUT_NAME %(weightout)s \
               -COMBINE_TYPE %(combine_type)s \
              """ % {
                  'imageout': header_only_file,
                  'weightout': "%s/preswarp.weight.fits" % (unique_singledir),
                  'combine_type': 'AVERAGE',
              }

        if (swarp_params['pixelscale'] > 0):
            swarp_opts += " -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.4f " % (swarp_params['pixelscale'])

            #
            # Factor in the potential change in pixelscale in the skylevel and 
            # gain computations
            #
            pixelscale_raw = 0.11
            pixelscale_binning = (swarp_params['pixelscale'] / pixelscale_raw)**2
            # do not multiply gain with binning, since flux/pixel is summed, not averages
            skylevel_list *= pixelscale_binning
            logger.info("Adjusting gain and skylevel for change in pixelscale (%.2f''/px --> bin=%.2f)" % (
                swarp_params['pixelscale'], pixelscale_binning))
        else:
            logger.info("Pixelscale remains unchanged, no need to adjust gain or skylevel")


        if (swarp_params['no-fluxscale']):
            swarp_opts += " -FSCALE_KEYWORD none "

        swarp_opts += " -SUBTRACT_BACK %s " % ("Y" if swarp_params['subtract_back']=='swarp' else "N")

        logger.debug("SWARP options for pre-stack:\n"+" ".join(swarp_opts.split()))


        # 
        # First create only the output header so we can pass some information 
        # to the user
        #
        swarp_cmd = "%(swarp)s %(opts)s -HEADER_ONLY Y %(files)s" % {
            'swarp': sitesetup.swarp_exec,
            'opts': swarp_opts,
            'files': " ".join(inputlist),
        }
        logger.debug("swarp_cmd=\n"+swarp_cmd)


        try:
            logger.info("Computing sky-coverage ...")
            logger.debug("Computing preswarp:\n%s" % (" ".join(swarp_cmd.split())))
            ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE) #, shell=True)
            # if retcode < 0:
            #     print >>sys.stderr, "Child was terminated by signal", -retcode
            # else:
            #     print >>sys.stderr, "Child returned", retcode
            #print retcode.stdout.readlines()
            #print retcode.stderr.readlines()
            (swarp_stdout, swarp_stderr) = ret.communicate()
            if (sitesetup.log_shell_output):
                logger.debug("\nCommand:\n%s\n--> Returncode: %d\n---\nStd.Out:\n%s\n---\nStd.Err:\n%s\n---", 
                             swarp_cmd, ret.returncode, swarp_stdout, swarp_stderr)

            # logger.debug("swarp stdout:\n"+swarp_stdout)
            # if (len(swarp_stderr) > 0 and ret.returncode != 0):
            #     logger.warning("swarp stderr:\n"+swarp_stderr)
            # else:
            #     logger.debug("swarp stderr:\n"+swarp_stderr)
        except OSError as e:
            podi_logging.log_exception()
            print("Execution failed:", e, file=sys.stderr)

        #
        # some information about the resulting stack is in the output-file
        #

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                output_info = pyfits.open(header_only_file)
        except IOError:
            podi_logging.log_exception()
            logger.error("Couldn't open the pre-swarp file, aborting")
            return
            
        logger.info("Stack information...")
        logger.info("   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header))

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']

        output_info.close()

        outimage_npixels = (out_naxis1 * out_naxis2) 
        if (outimage_npixels > 1e9):
            outimage_npixels /= 1e9
            if (swarp_params['huge_frame_allowed']):
                logger.warning("The output file is going to be huge (%.3f GigaPixels), continuing" % (outimage_npixels))
            else:
                logger.error("The output file exceeds the maximum allowed image size (%.3f GigaPixels)" % (outimage_npixels))
                cleanup_singles(swarp_params['unique_singledir'], logger)
                return None
        
        if (swarp_params['pixelscale'] <= 0):
            swarp_params['pixelscale'] = math.fabs(output_info[0].header['CD1_1']) * 3600.
            #pixelscale = (output_info[0].header['CD1_1'] * output_info[0].header['CD2_2'] \
            #             - output_info[0].header['CD1_2'] * output_info[0].header['CD2_1']) * 3600.
            logger.info("Computing pixelscale from data: %.4f arcsec/pixel" % (swarp_params['pixelscale']))
    
    if (swarp_params['preswarp_only']):
        # For only the pre-swarp frame, we are done here
        outputfile = outputfile+".fits"
        clobberfile(outputfile)
        #hdulist = pyfits.open(header_only_file)
        # add association table
        #hdulist.writeto(outputfile)
        shutil.copy(header_only_file, outputfile)
        return

    #############################################################################
    #
    # Prepare the individual frames, rectified and re-projected 
    # to the final grid
    #
    #############################################################################
    logger = logging.getLogger("SwarpStack - Singles")
    single_prepared_files = []

    # Prepare the worker queue
    sgl_queue = multiprocessing.JoinableQueue()

    # print inputlist
    # sys.exit(0)

    logger.info("Preparing focal-plane frames for %d frames..." % (len(inputlist)))
    for i in range(len(inputlist)):
        prepared_file = inputlist[i]
        # for prepared_file in inputlist:
        nonsidereal_dradec = nonsidereal_offsets[i]
        # print prepared_file, nonsidereal_dradec

        hdulist = pyfits.open(prepared_file)
        obsid = hdulist[0].header['OBSID']
        hdulist.close()

        fluxscale_kw = 'XXXXXXXX'
        magzero = hdulist[0].header['PHOTZP_X'] if 'PHOTZP_X' in hdulist[0].header else -99.
        # print magzero, swarp_params['no-fluxscale']
        # if (magzero > 0 and not swarp_params['no-fluxscale']):
        #     fluxscale_value = math.pow(10, 0.4*(swarp_params['target_magzero']-magzero))
        # else:
        #     exptime = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header else (
        #         hdulist[0].header['EXPTIME'] if 'EXPTIME' in hdulist[0].header else 1.0)
        #     fluxscale_value = 1./exptime

        # assemble all swarp options for that run
        dic = {'singledir': unique_singledir, #sitesetup.swarp_singledir,
               'obsid': obsid,
               'pixelscale': swarp_params['pixelscale'],
               'pixelscale_type': "MANUAL" if swarp_params['pixelscale'] > 0 else "MEDIAN",
               'center_ra': out_crval1,
               'center_dec': out_crval2,
               'imgsizex': out_naxis1,
               'imgsizey': out_naxis2,
               'resample_dir': unique_singledir,
               'inputfile': prepared_file,
               'swarp_default': swarp_default,
               'fluxscale_kw': fluxscale_kw, #'none' if swarp_params['no-fluxscale'] else 'FLXSCALE'
               'fluxscale_value': 1.0, #fluxscale_value,
               'fileid': i+1,
               'delete_tmpfiles': "N" if swarp_params['keep_resamp_files'] else "Y",
           }

        single_file = "%(singledir)s/%(obsid)s.%(fileid)d.fits" % dic

        swarp_opts = """
                 -c %(swarp_default)s 
                 -IMAGEOUT_NAME %(singledir)s/%(obsid)s.%(fileid)d.fits 
                 -WEIGHTOUT_NAME %(singledir)s/%(obsid)s.%(fileid)d.weight.fits 
                 -PIXEL_SCALE %(pixelscale)f 
                 -PIXELSCALE_TYPE %(pixelscale_type)s 
                 -COMBINE Y 
                 -COMBINE_TYPE WEIGHTED
                 -CENTER_TYPE MANUAL 
                 -CENTER %(center_ra)f,%(center_dec)f 
                 -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d 
                 -RESAMPLE_DIR %(resample_dir)s 
                 -SUBTRACT_BACK N 
                 -FSCALE_KEYWORD %(fluxscale_kw)s 
                 -FSCALE_DEFAULT %(fluxscale_value).10e 
                 -WEIGHT_TYPE MAP_WEIGHT
                 -WEIGHT_SUFFIX .weight.fits 
                 -RESCALE_WEIGHTS N
                 -DELETE_TMPFILES %(delete_tmpfiles)s \
                 %(inputfile)s 
                 """ % dic

#                 -WEIGHT_TYPE MAP_WEIGHT 

#                 -WEIGHT_THRESH 5

        # print swarp_opts
        swarp_cmd = "%s %s" % (sitesetup.swarp_exec, swarp_opts)

        if (add_only and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, skipping it" % (single_file))
        elif (swarp_params['reuse_singles'] and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, re-using it" % (single_file))
            single_prepared_files.append(single_file)
        else:
            logger.debug("Queuing single file %s:\n%s" % (
                prepared_file, " ".join(swarp_cmd.split()))
            )

            # print (swarp_cmd, prepared_file, single_file, swarp_params, nonsidereal_dradec, True)

            sgl_queue.put( (swarp_cmd, prepared_file, single_file, swarp_params, nonsidereal_dradec, True) )
            single_prepared_files.append(single_file)
            # time.sleep(2)

    #
    # Execute all swarps to create the single files
    #

    # Now with all swarp-runs queued, start a number of processes
    worker_args = (sgl_queue, "")
    processes = []
    for i in range(sitesetup.number_cpus):
        p = multiprocessing.Process(target=mp_swarp_single, args=worker_args)
        p.start()
        processes.append(p)

        # also add a quit-command for each process
        sgl_queue.put(None)
        
    # wait until all work is done
    sgl_queue.join()
    # join/terminate all processes
    for p in processes:
        p.join()

    #
    # If in "add" mode, rename the previous output file and add it to the list of input files
    #
    if (add_only):

        if (len(single_prepared_files) < 1):
            logger.info("No new files were added, so there's nothing to do.")
            return

        prev = 1
        while (True):
            filename = "%s.prev%02d.fits" % (outputfile, prev)
            if (not os.path.isfile(filename)):
                break
            prev += 1
            continue
                
        # Rename the current output file and its weights
        old_stacked = "%s.prev%02d.fits" % (outputfile, prev)
        old_weight = "%s.prev%02d.weight.fits" % (outputfile, prev)

        os.rename(outputfile+".fits", old_stacked)
        logger.debug("renamed old stack %s -> %s" % (outputfile+".fits", old_stacked))

        os.rename(outputfile+".weight.fits", old_weight)
        logger.debug("renamed old stack weight %s -> %s" % (outputfile+".weight.fits", old_weight))

        # Also add the new re-named old stacked file to list of input files
        single_prepared_files.append(old_stacked)
        logger.debug("Adding now old stack file to input list")
    #
    # Done re-naming the old file
    #

    # sys.exit(0)

    #############################################################################
    #
    # Now perform the background subtraction on each of the prepared single files
    # The resulting output files can then be used for the final stack without
    # additional resampling
    #
    #############################################################################

    #
    # Now use some brains to figure out the best way of setting the background 
    # subtraction to get s nice smooth background that does not over-subtract 
    # the target.
    #
    bg_opts = ""
    if (swarp_params['target_dimension'] > 0 and swarp_params['subtract_back']=='swarp'):
        dic['BACK_TYPE'] = "AUTO"
        dic['BACK_SIZE'] = 128
        dic['BACK_FILTERSIZE'] = 3

        # Rule of thum: larger objects: make both filtersize and back_size 
        # larger
        # first compute a reference size for the default settings
        ref_size = dic['BACK_SIZE'] * dic['BACK_FILTERSIZE'] \
                   * swarp_params['pixelscale'] / 60. \
                   * 0.1  # the last factor is a fudge-factor
        logger.debug("Reference size: %f" % (ref_size))
        # Now scale up the filtersize, making sure it stays between 3 and 7
        filtersize = int(math.floor(math.sqrt(swarp_params['target_dimension'] / ref_size) * dic['BACK_FILTERSIZE']))
        logger.debug("Simple filter size: %d" % (filtersize))
        if (filtersize < 3): filtersize = 3
        if (filtersize > 7): filtersize = 7

        # in a next step, modify the backsize parameter. Make sure it does not
        # become too large or too small
        backsize = (swarp_params['target_dimension'] * 60. / swarp_params['pixelscale']) / filtersize

        logger.debug("BACK-SIZE: %f" % (backsize))
        if (backsize < 64): backsize = 64
        if (backsize > 600): backsize = 600

        dic['BACK_SIZE'] = backsize
        dic['BACK_FILTERSIZE'] = filtersize

        bg_opts = """  
               -BACK_TYPE %(BACK_TYPE)s
               -BACK_SIZE %(BACK_SIZE)d
               -BACK_FILTERSIZE %(BACK_FILTERSIZE)d
        """ % dic
        logger.debug("Adding background parameters:\n\n"+bg_opts+"\n\n")
        # swarp_opts += bg_opts


    if (swarp_params['subtract_back']=='swarp'):
        logger = logging.getLogger("SwarpStack - SkySub")
        logger.info("Performing sky-subtraction on all frames")

        final_prepared_files = []

        # Prepare the worker queue
        sgl_queue = multiprocessing.JoinableQueue()

        fileid = 0
        for prepared_file in single_prepared_files:
            hdulist = pyfits.open(prepared_file)
            obsid = hdulist[0].header['OBSID']
            hdulist.close()
            fileid += 1

            # assemble all swarp options for that run
            bgsub_file = "%(singledir)s/%(obsid)s.%(fileid)d.bgsub.fits" % {
                'singledir': unique_singledir,
                'obsid': obsid,
                'fileid': fileid,}
            bgsub_weight_file = "%(singledir)s/%(obsid)s.%(fileid)d.bgsub.weight.fits" % {
                'singledir': unique_singledir,
                'obsid': obsid,
                'fileid': fileid,}

            dic = {'singledir': unique_singledir,
                   'obsid': obsid,
                   'pixelscale': swarp_params['pixelscale'],
                   'pixelscale_type': "MANUAL" if swarp_params['pixelscale'] > 0 else "MEDIAN",
                   'center_ra': out_crval1,
                   'center_dec': out_crval2,
                   'imgsizex': out_naxis1,
                   'imgsizey': out_naxis2,
                   'resample_dir': unique_singledir,
                   'inputfile': prepared_file,
                   'swarp_default': swarp_default,
                   'bgsub': "Y" if swarp_params['subtract_back']=='swarp' else "N",
                   'bgsub_file': bgsub_file,
                   'bgsub_weight_file': bgsub_weight_file,
                   'bgopts': bg_opts,
                   'fileid': fileid,
                   'delete_tmpfiles': "N" if swarp_params['keep_resamp_files'] else "Y",
               }
            
            swarp_opts = """\
                     -c %(swarp_default)s \
                     -IMAGEOUT_NAME %(bgsub_file)s \
                     -WEIGHTOUT_NAME %(bgsub_weight_file)s \
                     -PIXEL_SCALE %(pixelscale)f \
                     -PIXELSCALE_TYPE %(pixelscale_type)s \
                     -COMBINE Y \
                     -COMBINE_TYPE WEIGHTED \
                     -CENTER_TYPE MANUAL \
                     -CENTER %(center_ra)f,%(center_dec)f \
                     -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
                     -RESAMPLE_DIR %(resample_dir)s \
                     -RESAMPLE Y \
                     -SUBTRACT_BACK %(bgsub)s \
                     -WEIGHT_TYPE MAP_WEIGHT \
                     -WEIGHT_SUFFIX .weight.fits \
                     -RESCALE_WEIGHTS N \
                     -DELETE_TMPFILES %(delete_tmpfiles)s \
                     -WRITE_FILEINFO Y
                     -FSCALE_KEYWORD XXXXXXXX \
                     -FSCALE_DEFAULT 1.0 \
                     %(bgopts)s \
                     %(inputfile)s \
                     """ % dic

            
            # print swarp_opts
            swarp_cmd = "%s %s" % (sitesetup.swarp_exec, swarp_opts)
            # print swarp_cmd

            # Disable the mask, since we already prepared it earlier
            swarp_params['mask'] = None

            if (add_only and os.path.isfile(bgsub_file)):
                logger.info("This single-swarped file (%s) exist, skipping it" % (bgsub_file))
            elif (swarp_params['reuse_singles'] and os.path.isfile(bgsub_file)):
                logger.info("This single-swarped file (%s) exist, re-using it" % (bgsub_file))
                final_prepared_files.append(bgsub_file)
            else:
                logger.info("Preparing file %s, please wait ..." % (prepared_file))
                logger.debug(" ".join(swarp_cmd.split()))

                sgl_queue.put( (swarp_cmd, prepared_file, bgsub_file, swarp_params, None, False) )
                final_prepared_files.append(bgsub_file)

        #
        # Execute all swarps to create the single files
        #

        # Now with all swarp-runs queued, start a number of processes
        worker_args = (sgl_queue, "")
        processes = []
        for i in range(sitesetup.number_cpus):
            p = multiprocessing.Process(target=mp_swarp_single, args=worker_args)
            p.start()
            processes.append(p)

            # also add a quit-command for each process
            sgl_queue.put(None)

        # wait until all work is done
        sgl_queue.join()
        # join/terminate all processes
        for p in processes:
            p.join()

    else:
        # No background subtraction was requested
        final_prepared_files = single_prepared_files

    logging.debug("files to stack: %s" % (str(final_prepared_files)))

    # sys.exit(0)

    if (swarp_params['median_reject']):
        swarp_params['combine-type'].insert(0, "MEDIAN")

    #############################################################################
    #
    # Now all single files are prepared, go ahead and produce the actual stack
    # Use the background-subtracted or single files from above as direct input, 
    # i.e. do not resample these files as they are already on the right pixel grid.
    #
    #############################################################################
    for idx, combine_type in enumerate(swarp_params['combine-type']):
        dic['combine_type'] = combine_type #swarp_params['combine-type'] #"AVERAGE"
        dic['imageout'] = "%s.%s.fits" % (outputfile, combine_type)
        dic['weightout'] = "%s.%s.weight.fits" % (outputfile, combine_type) #outputfile+".weight.fits"
        dic['prepared_files'] = " ".join(single_prepared_files)
        dic['bgsub'] = "N" # as this was done before if swarp_params['subtract_back'] else "N"
        dic['clip-sigma'] = swarp_params['clip-sigma']
        dic['clip-ampfrac'] = swarp_params['clip-ampfrac']

        swarp_opts = """
                     -c %(swarp_default)s 
                     -IMAGEOUT_NAME %(imageout)s 
                     -WEIGHTOUT_NAME %(weightout)s 
                     -COMBINE_TYPE %(combine_type)s 
                     -PIXEL_SCALE %(pixelscale)f 
                     -PIXELSCALE_TYPE %(pixelscale_type)s 
                     -COMBINE Y 
                     -COMBINE_TYPE %(combine_type)s 
                     -CLIP_AMPFRAC %(clip-ampfrac)f 
                     -CLIP_SIGMA %(clip-sigma)f 
                     -CENTER_TYPE MANUAL 
                     -CENTER %(center_ra)f,%(center_dec)f 
                     -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d 
                     -RESAMPLE N 
                     -RESAMPLE_DIR %(singledir)s 
                     -SUBTRACT_BACK %(bgsub)s 
                     -WEIGHT_TYPE MAP_WEIGHT 
                     -WEIGHT_SUFFIX .weight.fits 
                     -RESCALE_WEIGHTS N 
                     -DELETE_TMPFILES N 
                     -WRITE_FILEINFO Y
                     """ % dic

        logger = logging.getLogger("SwarpStack - FinalStack")
        logger.info("Starting final stacking (%s) ..." % (combine_type))
        # print swarp_opts

        swarp_cmd = "%s %s %s" % (sitesetup.swarp_exec, swarp_opts, " ".join(final_prepared_files))
        logger.debug("swarp-options:%s"  %(swarp_cmd))
        logger.debug("\n"+" ".join(swarp_cmd.split()))
        try:
            ret = subprocess.Popen(swarp_cmd.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
            (swarp_stdout, swarp_stderr) = ret.communicate()
            if (sitesetup.log_shell_output):
                logger.debug("\nCommand:\n%s\n--> Returncode: %d\n---\nStd.Out:\n%s\n---\nStd.Err:\n%s\n---", 
                             swarp_cmd, ret.returncode, swarp_stdout, swarp_stderr)

            # logger.debug("swarp stdout:\n"+swarp_stdout)
            # if (len(swarp_stderr) > 0 and ret.returncode != 0):
            #     logger.warning("swarp stderr:\n"+swarp_stderr)
            # else:
            #     logger.debug("swarp stderr:\n"+swarp_stderr)
            #print "\n".join(swarp_stderr)
            logger.info("done, swarp returned (ret-code: %d)!" % ret.returncode)
        except OSError as e:
            podi_logging.log_exception()
            print("Execution failed:", e, file=sys.stderr)

        logger.info("Stack (%s) complete, adding headers" % (dic['imageout']))


        # Finally, open the output file and copy a bunch of headers into it
        try:
            hdustack = pyfits.open(dic['imageout'], mode='update')
        except:
            logger.error("Unable to open %s" % (dic['imageout']))
            continue

        # Also open the first frame in the stack to be used as data source

        #
        # Add some exposure-specific data from the reference header back to the 
        # output file
        #
        #firsthdu = pyfits.open(inputlist[0])
        if ('FILE0001' in hdustack[0].header):
            add_fits_header_title(hdustack[0].header, "Input files as written by swarp", 'FILE0001')

        for hdrkey in [
                'TARGRA', 'TARGDEC',
                'FILTER', 'FILTERID', 'FILTDSCR', 
                'OBSID', 'OBJECT', 
                'DATE-OBS', 'TIME-OBS', 'MJD-OBS']:
            if (hdrkey in reference_header):
                key, val, com = reference_header.cards[hdrkey]
                hdustack[0].header[key] = (val, com)
        add_fits_header_title(hdustack[0].header, "Headers inherited from input frames", 'TARGRA')

        del hdustack[0].header['EXPTIME']
        hdustack[0].header['EXPTIME'] = (stack_total_exptime, "total exposure time in stack")
        hdustack[0].header['MJD-STRT'] = (stack_start_time, "MJD at start of earliest exposure")
        hdustack[0].header['MJD-END'] = (stack_end_time, "MJD at end of last exposure")
        hdustack[0].header['NCOMBINE'] = (stack_framecount, "number of exposures in stack")

        min_saturation_level = numpy.min(all_saturation)
        hdustack[0].header['SATURATE'] = (min_saturation_level,
                                          "min saturation level")
        hdustack[0].header['SATR8AVG'] = (numpy.mean(all_saturation),
                                          "mean saturation level")
        hdustack[0].header['SATR8MAX'] = (numpy.max(all_saturation),
                                          "max saturation level")
        hdustack[0].header['SATR8MED'] = (numpy.median(all_saturation),
                                          "median saturation level")

        add_fits_header_title(hdustack[0].header, "Computed timing information", 'EXPTIME')

        # Add some additional headers
        hdustack[0].header['MAGZERO']  = (swarp_params['target_magzero'],
                                          "after flux-scaling each input exposure")
        hdustack[0].header['BACKGSUB'] = (swarp_params['subtract_back'] if swarp_params['subtract_back'] else "none",
                                          "was background subtracted?")
        hdustack[0].header['PIXLSCAL'] = (swarp_params['pixelscale'],
                                          "user-selected pixelscale")
        hdustack[0].header['REUSESGL'] = ("yes" if swarp_params['reuse_singles'] else "no",
                                          "reuse singles?")

        # Set RA/DEC keywords to make WCS work in IRAF
        pos = ephem.Equatorial(numpy.radians(hdustack[0].header['CRVAL1']),
                               numpy.radians(hdustack[0].header['CRVAL2']))
        hdustack[0].header['RA'] = (str(pos.ra), "pointing RA")
        hdustack[0].header['DEC'] = (str(pos.dec), "pointing DEC")

        #
        #
        # Add all photometry-relevant derived keywords in the output stack
        #
        #
        if (combine_type in ['SUM']):
            stack_gain = numpy.mean(gain_list)
            stack_skylevel = numpy.sum(skylevel_list)
        else:
            stack_gain = numpy.sum(gain_list)
            if (combine_type in ['WEIGHTED']):
                # compute weighted sky-level
                stack_skylevel = numpy.sum(skylevel_list*weight_list) / numpy.sum(weight_list)
            else:
                stack_skylevel = numpy.mean(skylevel_list)
        for hk in ['GAIN', 'SKYLEVEL']:
            if (hk in hdustack[0].header):
                del hdustack[0].header[hk]
        hdustack[0].header['GAIN'] = (stack_gain, "combined GAIN")
        hdustack[0].header['SKYLEVEL'] = (stack_skylevel, "combined skylevel")
        logger.info("Stack: GAIN=%f - SKY=%f" % (stack_gain, stack_skylevel))

        # Check if we are to add the sky-level back into the stack
        hdustack[0].header['NORMSKY'] = (False, "re-normalize sky")
        if (params['normalize_sky'] and not params['subtract_back'] == False):
            logger.info("Adding back skylevel (%f) to stack" % (stack_skylevel))
            hdustack[0].data += stack_skylevel
            hdustack[0].header['NORMSKY'] = True
            
        #
        # Store all configuration about non-sidereal / ephemerides corrections 
        # that were applied during execution
        #
        nonsid_mode = 'none'
        if (swarp_params['use_nonsidereal']): nonsid_mode = 'nonsidereal'
        if ('ephemerides' in swarp_params and swarp_params['ephemerides'] is not None): nonsid_mode = "ephemerides"
        hdustack[0].header['NONSIDRL'] = (nonsid_mode,
                                          "Non-sidereal correction")
        valid_nonsidereal_reference = None
        try:
            hdustack[0].header['NSID_RA']  = (options['nonsidereal']['dra'],
                                              "non-sidereal rate dRA*cosDec [arcsec/hr]")
            hdustack[0].header['NSID_DEC'] = (options['nonsidereal']['ddec'],
                                              "non-sidereal rate dDec [arcsec/hr]")
            hdustack[0].header['NSID_RT']  = (options['nonsidereal']['ref_mjd'],
                                              "non-sidereal reference MJD")
            hdustack[0].header['NSID_RF']  = (options['nonsidereal']['ref'],
                                              "non-sidereal ref. from cmd line")

            hdustack[0].header['NSREFMJD'] = options['nonsidereal']['ref_mjd']
            hdustack[0].header['NSREFILE'] = os.path.abspath(options['nonsidereal']['ref'])
            hdustack[0].header['NSREFOBS'] = options['nonsidereal']['ref_obsid']

            if (os.path.isfile(options['nonsidereal']['ref'])):
                valid_nonsidereal_reference = options['nonsidereal']['ref']
        except:
            pass

        #
        # Check if we use the ephem mode; if so add some more headers
        #
        try:
            eph = swarp_params['ephemerides']
            hdustack[0].header['EPHM-MJD'] = eph['ref-mjd']
            hdustack[0].header['EPHM-REF'] = eph['ref']
            hdustack[0].header['EPHMMODE'] = eph['mode']
            hdustack[0].header['EPHMTRGT'] = eph['target']
            hdustack[0].header['EPHMFILE'] = eph['datafile']
            hdustack[0].header['NSREFMJD'] = eph['ref-mjd']
            hdustack[0].header['NSREFILE'] = os.path.abspath(eph['ref'])
            hdustack[0].header['NSREFOBS'] = eph['ref-obsid']

            if (os.path.isfile(eph['ref'])):
                valid_nonsidereal_reference = eph['ref']
        except:
            pass
            
        add_fits_header_title(hdustack[0].header, "Non-sidereal/ephemerides configuration", 'NONSIDRL')


        add_fits_header_title(hdustack[0].header, "SwarpStack parameters supplied by user", 'MAGZERO')

        if (valid_nonsidereal_reference is not None):
            master_reduction_files_used = \
                podi_associations.collect_reduction_files_used(master_reduction_files_used, 
                                             {"mjd-reference": valid_nonsidereal_reference})
            firsthdu = pyfits.open(valid_nonsidereal_reference)
            try:
                hdustack[0].header['NSR-OBST'] = (firsthdu[0].header['TIME-OBS'], 
                                                  "reference TIME-OBS")
                hdustack[0].header['NSR-OBSD'] = (firsthdu[0].header['DATE-OBS'], 
                                                  "reference DATE-OBS")
                hdustack[0].header['NSR-OBSM'] = (firsthdu[0].header['MJD-OBS'], 
                                                  "reference MJD-OBS")

                hdustack[0].header['NSR-MIDT'] = (firsthdu[0].header['TIME-MID'], 
                                                  "reference TIME-MID")
                hdustack[0].header['NSR-MIDD'] = (firsthdu[0].header['DATE-MID'], 
                                                  "reference DATE-MID")
                hdustack[0].header['NSR-MIDM'] = (firsthdu[0].header['MJD-MID'], 
                                                  "reference MJD-MID")

                hdustack[0].header['NSR-ENDT'] = (firsthdu[0].header['TIME-END'], 
                                                  "reference TIME-END")
                hdustack[0].header['NSR-ENDD'] = (firsthdu[0].header['DATE-END'], 
                                                  "reference DATE-END")
                hdustack[0].header['NSR-ENDM'] = (firsthdu[0].header['MJD-END'], 
                                                  "reference MJD-END")
                add_fits_header_title(hdustack[0].header, 
                                      "Timing information from reference frame", 'NSR-OBST')
            except:
                pass
            firsthdu.close()

        # Add the user-defined keywords to the stacked file. This is required for
        # proper integration with the PPA framework.
        # print options['additional_fits_headers']
        first_useradded_key = None
        for key, value in options['additional_fits_headers'].items():
            if (key in hdustack[0].header):
                # if key exists, leave comment unchanged
                hdustack[0].header[key] = value
            else:
                if (first_useradded_key is None): first_useradded_key = key
                hdustack[0].header[key] = (value, "user-added keyword")
        if (len(options['additional_fits_headers']) > 0 and first_useradded_key is not None):
            add_fits_header_title(hdustack[0].header, 
                                  "user-added keywords", first_useradded_key)

        #
        # Create an association table from the master reduction files used.
        # 
        # print master_reduction_files_used
        assoc_table = podi_associations.create_association_table(master_reduction_files_used)
        hdustack.append(assoc_table)

        hdustack.flush()
        hdustack.close()
    
        if (swarp_params['median_reject'] and idx==0):
            from swarp_filter import mask_outliers_in_stack
            print("**\n"*3,"median reject starting!")
            mask_outliers_in_stack(median_frame=dic['imageout'],
                          singles=final_prepared_files)
            # rename the MEDIAN file to reference file
            fn_base = dic['imageout'][:-12]
            try:
                shutil.move(fn_base+".MEDIAN.fits", fn_base+".REFERENCE.fits")
                shutil.move(fn_base+".MEDIAN.weight.fits", fn_base+".REFERENCE.weight.fits")
            except:
                pass
            print("rejectng done!", "\n**"*3)
            # print "**\n"*3,"median reject comes now","\n**"*3

    # 
    # Delete all temporary files and the temp-directory
    # 
    if (not keep_intermediates):
        logger.info("Deleting all intermediate files")
        try:
            shutil.rmtree(unique_singledir)
        except:
            logger.error("There was a problem with recursively deleting the temp directory")
            podi_logging.log_exception()
            pass
 
    logger.info("All done!")

    return modified_files, single_prepared_files, final_prepared_files, unique_singledir


def load_horizons_ephems(object_name, ref_file, ref_mjd, filelist, params):

    if ((object_name.startswith("'") and object_name.endswith("'")) or
        (object_name.startswith('"') and object_name.endswith('"'))):
        object_name = object_name[1:-1]
    object_name = object_name.replace(".", " ")

    results = podi_ephemerides.get_ephemerides_for_object_from_filelist(
        object_name=object_name, 
        filelist=filelist,
        session_log_file="swarpstack_horizon.session",
        verbose=False
    )

    params['ephemerides'] = {
        'ref': ref_file,
        'ref-mjd': ref_mjd,
        'datafile': "NASA-TELNET",
        'data': results['data'],
        'ra': results['ra'],
        'dec': results['dec'],
        'target': object_name,
    }
    return

def get_reference_mjd(item):

    ref_obsid = None
    try:
        ref_mjd = float(item)
    except:
        if (os.path.isfile(item)):
            hdulist = pyfits.open(item)
            ref_obsid = hdulist[0].header['OBSID'] if 'OBSID' in hdulist[0].header else None
            for ext in hdulist:
                if ('MJD-OBS' in ext.header):
                    ref_mjd = ext.header['MJD-OBS']
                    # correct reference MJD to the mid-point of the exposure
                    if ('EXPTIME' in ext.header):
                        ref_mjd += ext.header['EXPTIME']/2./86400
                    break
            hdulist.close()
        else:
            ref_mjd = None

    return ref_mjd, ref_obsid


def read_swarp_params(filelist):

    params = {}

    logger = logging.getLogger("SwarpStack - Config")

    params['pixelscale'] = float(cmdline_arg_set_or_default("-pixelscale", 0))

    # params['subtract_back'] = cmdline_arg_isset("-bgsub")
    params['subtract_back'] = False
    params['sky_regions'] = None
    params['sky_regions_file'] = None
    if (cmdline_arg_isset("-bgsub")):
        params['subtract_back'] = cmdline_arg_set_or_default("-bgsub", 'swarp')

        if (params['subtract_back'] is not None and
            os.path.isfile(params['subtract_back'])):
            params['sky_regions_file'] = params['subtract_back']
            params['subtract_back'] = "_REGIONS_"
            params['sky_regions'] = read_sky_regions_file(params['sky_regions_file'])

    params['reuse_singles'] = cmdline_arg_isset("-reusesingles")
    params['use_nonsidereal'] = cmdline_arg_isset("-nonsidereal")
    params['target_dimension'] = float(cmdline_arg_set_or_default('-dimension', -1))
    params['add'] = cmdline_arg_isset("-add")
    params['reference_file'] = cmdline_arg_set_or_default("-reference", None)
    params['no-fluxscale'] = cmdline_arg_isset('-nofluxscale')

    params['ignore_quality_checks'] = cmdline_arg_isset("-ignorechecks")

    combine_methods = cmdline_arg_set_or_default('-combine', 'weighted')
    params['combine-type'] = []
    for combine_method in combine_methods.split(","):
        if (not combine_method.lower() in ['average', 'median', 'sum', 'min', 'max', 'weighted', 'chi2',
                                           'chi-old', 'chi-mode', 'chi-mean', 'clipped',
                                           'weighted_weight', 'median_weight', 'and', 'nand', 'or', 'nor',
                                           'sigmaclipmean', 'sigmaclipwmean', 'sigmaclipmedian']):
            logger = logging.getLogger("Setup")
            logger.error("The specified combine method (%s) is not supported, using average instead" % (combine_method))
            continue
        else:
            params['combine-type'].append(combine_method.upper())
    if (params['combine-type'] == []):
        params['combine-type'].append('average'.upper())
    

    params['use_ephemerides'] = cmdline_arg_isset("-ephemerides")
    if (params['use_ephemerides']):
        opts = cmdline_arg_set_or_default('-ephemerides', None)
        # print opts
        if (opts is None):
            params['use_ephemerides'] = False
        else:
            items = opts.split(',')
            # print items

            # See if the first parameter is a number, if so it's the reference MJD,
            # if not, assume it's a FITS file and we are to read MJD-OBS from the header
            ref_mjd = None
            if (items[0] == 'NASA'):
                if (len(items) >= 3):
                    ref_file = items[2]
                else:
                    # Use the first input frame as reference frame
                    ref_file = get_clean_cmdline()[2]
                
                ref_mjd, ref_obsid = get_reference_mjd(ref_file)

                if (ref_mjd is not None):
                    object_name = items[1]
                    full_filelist = list(filelist)
                    full_filelist.append(ref_file)
                    load_horizons_ephems(object_name, ref_file, ref_mjd, full_filelist, params)
                    params['ephemerides']['mode'] = "telnet:horizons"
                    params['ephemerides']['ref-obsid'] = ref_obsid
                else:
                    logger.critical("Unable to find MJD reference file (%s)" % (ref_file))
                    raise RuntimeError("Unable to find MJD reference file (%s)" % (ref_file))

            elif (items[0] == 'file'):
                ref_file = get_clean_cmdline()[2] if len(items) <=3 else items[2]
                ref_mjd, ref_obsid = get_reference_mjd(ref_file)
                if (ref_mjd is not None):

                    # Now read and process the datafile
                    import ephemerides
                    ra, dec, data = ephemerides.load_ephemerides(
                        items[1], plot=False
                    )
                    logger.info("Using ephemerides from file %s" % (items[1]))

                    params['ephemerides'] = {
                        'ref': items[0],
                        'ref-mjd': ref_mjd,
                        'datafile': os.path.abspath(items[1]),
                        'data': data,
                        'ra': ra,
                        'dec': dec,
                        'mode': 'file:'+items[1],
                        'target': "user-defined",
                        'ref-obsid': ref_obsid,
                    }
            if (ref_mjd is None):
                params['use_ephemerides'] = False


    params['clip-ampfrac'] = 0.3
    params['clip-sigma'] = 4.0
    if (cmdline_arg_isset('-clip')):
        vals = get_cmdline_arg('-clip')
        if (len(vals)>0):
            items = vals.split(',')
            if (len(items) > 0):
                params['clip-sigma'] = float(items[0])
            if (len(items) >= 2):
                params['clip-ampfrac'] = float(items[1])

    params['mask-fits'] = None
    params['mask-npix'] = float(cmdline_arg_set_or_default("-masknpix", 5))
    params['mask-nsigma'] = float(cmdline_arg_set_or_default("-masknsigma", 2))
    params['mask-list'] = None
    params['mask-mjds'] = None
    if (cmdline_arg_isset("-maskframe")):
        vals = get_cmdline_arg("-maskframe").split(",")
        mask_list = []
        for filename in vals:
            if (os.path.isfile(filename)):
                mask_list.append(filename)
        if (mask_list): params['mask-fits'] = mask_list
    
    params['mask'] = None
    if (cmdline_arg_isset("-mask")):
        mask = get_cmdline_arg("-mask")
        if (os.path.isfile(mask)):
            params['mask'] = mask

    params['huge_frame_allowed'] = cmdline_arg_isset("-huge")

    params['target_magzero'] = cmdline_arg_set_or_default("-targetzp", 25.0)

    params['illumcorr'] = cmdline_arg_set_or_default("-illumcorr", None)
    params['illumcorrfiles'] = cmdline_arg_set_or_default("-illumcorrfiles", None)
    params['un_illumcorr'] = cmdline_arg_isset("-unic")


    params['normalize_sky'] = cmdline_arg_isset("-normsky")

    params['skip_guide_otas'] = cmdline_arg_isset("-rmguideota")

    params['wipe_cells'] = read_wipecells_list()

    params['aggressive_clean'] = cmdline_arg_isset('-ocdclean')
    params['keep_resamp_files'] = cmdline_arg_isset('-keepresamp')

    params['preswarp_only'] = cmdline_arg_isset('-preswarponly')

    params['median_reject'] = cmdline_arg_isset("-medianreject")

    params['mask_saturated'] = None
    if (cmdline_arg_isset("-masksaturated")):
        params['mask_saturated'] = int(cmdline_arg_set_or_default("-masksaturated", 58000))

    params['cutout_list'] = None
    if (cmdline_arg_isset("-cutout")):
        txt = get_cmdline_arg("-cutout")
        items = txt.split(",")
        cutout = {
            'ra':  float(items[0]),
            'dec': float(items[1]),
            'size': float(items[2]),
            }
        params['cutout_list'] = cutout

    params['autophotflat'] = False
    params['apf_timeframe'] = None
    params['apf_nframes'] = None
    params['apf_noauto'] = False
    params['apf_resolution'] = None

    if(cmdline_arg_isset('-autophotflat')):
        params['autophotflat'] = True
        opt = cmdline_arg_set_or_default2("-autophotflat", None)
        if (opt is not None):
            print(opt)

            apf_options = opt.split(":")
            for o in apf_options:
                _o = o.split("=")
                opt_name = _o[0]
                opt_value = _o[1]
                if (opt_name == "dt"): #o.startswith("dt=")):
                    interval = opt_value[-1]
                    if (interval == 'd'):
                        time = float(opt_value[:-1])
                    elif (interval == "h"):
                        time = float(opt_value[:-1])/24.
                    else:
                        time = float(opt_value)
                    params['apf_timeframe'] = time

                elif (opt_name == "n"):
                    params['apf_nframes'] = int(opt_value)

                elif (opt_name == "no_auto"):
                    params['apf_noauto'] = True

                elif (opt_name == "res"):
                    params['apf_resolution'] = float(opt_value)

    cheese_grate = cmdline_arg_set_or_default("-cheesegrate", None)
    if (cheese_grate is None):
        params['cheese_grate'] = None
    elif (cheese_grate == "auto"):
        params['cheese_grate'] = 'auto'
    else:
        try:
            params['cheese_grate'] = int(cheese_grate)
        except:
            params['cheese_grate'] = None

    return params

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        import podi_swarpstack as me
        print(me.__doc__)
        sys.exit(0)

    # Setup everything we need for logging
    options = set_default_options()
    podi_logging.setup_logging(options)

    logger = logging.getLogger("SwarpStack-Startup")
    while (cmdline_arg_isset("-fromfile")):
        configfile = get_cmdline_arg("-fromfile")
        if (os.path.isfile(configfile)):
            logger.info("Reading additional command line parameters from file (%s)" % (configfile))
            conf = open(configfile, "r")
            lines = conf.readlines()
            all_items = []
            for line in lines:
                line = line.strip()
                if (len(line) <= 0):
                    continue
                elif (line.startswith("#")):
                    continue
                elif (line.startswith("-")):
                    all_items.append(line)
                else:
                    items = line.split()
                    all_items.append(items[0])
            conf.close()
            for i in range(len(sys.argv)):
                if (sys.argv[i].startswith("-fromfile")):
                    del sys.argv[i]
                    for to_add in all_items[::-1]:
                        sys.argv.insert(i, to_add)
        else:
            logger.error("Can't open the configfile (%s)" % (configfile))

    logger.debug("Reading options from command line")
    options = read_options_from_commandline(options)

    keep_intermediates = cmdline_arg_isset('-keep')
    unique_dir = None
    if (cmdline_arg_isset("-uniquedir")):
        unique_dir = get_cmdline_arg("-uniquedir")

    # print "non-sid",options['nonsidereal']
    try:

        # Read command line and store all results in params dictionary
        outputfile = get_clean_cmdline()[1]

        raw_filelist = get_clean_cmdline()[2:]
        inputlist = read_list(raw_filelist, )
        # inputlist = []
        # for f in get_clean_cmdline()[2:]:
        #     if (os.path.isfile(f)):
        #         inputlist.append(f)
        params = read_swarp_params(inputlist)
        
        # print params
        # print inputlist

        logger.debug("Commanding output: %s" % (outputfile))
        for i in inputlist:        
            logger.debug("Commanding input: %s" % (i))

        logger.debug("Starting processing")
        swarpstack(outputfile=outputfile, 
                   inputlist=inputlist, 
                   swarp_params=params, 
                   options=options,
                   keep_intermediates=keep_intermediates,
                   unique_dir=unique_dir)

    except (KeyboardInterrupt, SystemExit) as e:
        pass
    except RuntimeError as e:
        logger.critical("Encountered critical error, terminating execution! Please check command!")
        logger.debug("Terminating after runtime error!")
    except:
        if (len(get_clean_cmdline()) < 3):
            logger.error("Not enough parameters have been specified, need at least output file and 1 input file")
        else:
            podi_logging.log_exception()
        pass
    finally:
        podi_logging.shutdown_logging(options)
