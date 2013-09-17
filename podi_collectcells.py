#!/usr/local/bin/python
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
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
import scipy
import scipy.optimize
import ephem
import traceback

import Queue
import threading
import multiprocessing
import ctypes
import time

from podi_plotting import *

gain_correct_frames = False
from podi_definitions import *
import podi_findstars
import podi_search_ipprefcat
import podi_fixwcs
import podi_fixwcs_rotation
import podi_sitesetup as sitesetup
import podi_crosstalk
import podi_persistency
import podi_asyncfitswrite
import podi_fitskybackground
import podi_matchcatalogs
import podi_matchpupilghost
import podi_fringing
import podi_photcalib
import podi_nonlinearity

from astLib import astWCS

fix_cpu_count = False

if (sitesetup.number_cpus == "auto"):
    try:
        number_cpus = multiprocessing.cpu_count()
        print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
        if (number_cpus > sitesetup.max_cpu_count and sitesetup.max_cpu_count > 1):
            number_cpus = sitesetup.max_cpu_count
            print "... but using only %d of them!" % (number_cpus)
    except:
        pass
else:
    number_cpus = sitesetup.number_cpus



def collect_reduce_ota(filename,
                       verbose=False,
                       options=None):

    data_products = {
        "hdu": None,
        "wcsdata": None,
        "sky-samples": None,
        "sky": None,
        }

    if (not os.path.isfile(filename)):
        stdout_write("Couldn't find file %s ..." % (filename))
    else:
        # Create an fits extension to hold the output
        hdu = pyfits.ImageHDU()
        log_svn_version(hdu.header)

        try:
            hdulist = pyfits.open(filename, memmap=False)
        except:
            # This happed with corrupt files
            return data_products

        detsize = break_region_string(hdulist[0].header['DETSIZE'])
        det_x1, det_x2, det_y1, det_y2 = detsize
        #print det_x1, det_x2, det_y1, det_y2

        size_x, size_y = det_x2 - det_x1 + 1, det_y2 - det_y1 + 1
        #print size_x, size_y
        size_x, size_y = 4096, 4096
        #print size_x, size_y

        obsid = hdulist[0].header["OBSID"]
        ota = int(hdulist[0].header['FPPOS'][2:])
        ota_c_x, ota_c_y = int(math.floor(ota/10)), int(math.fmod(ota,10))

        # Save the fppos as name for this extension
        ota_name = "OTA%02d" % ota
        extname = "OTA%02d.SCI" % ota
        hdu.update_ext_name(extname)
        
        # Now copy the headers from the original file into the new one
        cards = hdulist[0].header.ascardlist()
        for c in cards:
            hdu.header.update(c.key, c.value, c.comment)

        # Also read the MJD for this frame. This will be needed later for the correction
        mjd = hdu.header['MJD-OBS']

        # 
        # Perform cross-talk correction, using coefficients found in the 
        # podi_crosstalk package.
        #
        # Procedure: 
        # Correct all frames for the crosstalk, then write them back into the 
        # original OTA position so we can perform the overscan subtraction etc.
        #

        # Allocate some memory for the cells in one row
        xtalk_corr = [None] * 8

        # Now go through each of the 8 lines
        for row in range(8):
            for column in range(8):
                
                # Allocate some more memory to hold the output of the cross-talk 
                # correction of this frame
                xtalk_corr[column] = numpy.zeros(hdulist[1].data.shape)

                for xtalk_column in range(8):
                    # Construct the name of each cell that's causing the crosstalk
                    xy_name = "xy%d%d" % (xtalk_column, row)

                    # Now go through the list of all cells in this row and add them to 
                    # the corrected cell content output
                    #print "Adding ",xy_name,"to ",extname, column, row, "(scaling",podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column],")"

                    correction = hdulist[xy_name].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                    if (column != xtalk_column):
                        saturated = hdulist[xy_name].data >= podi_crosstalk.xtalk_saturation_limit
                        correction[saturated] = -1 * podi_crosstalk.xtalk_saturated_correction

                    xtalk_corr[column] += correction #hdulist[xy_name].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                    #print xtalk_corr[column][100,100]

            for column in range(8):
                # Now all cells in this row have been corrected, let's write them 
                # back into the hdulist so can can continue with the overscan subtraction etc.
                xy_name = "xy%d%d" % (column, row)
                hdulist[xy_name].data = xtalk_corr[column]

        #
        # Allocate memory for the merged frame, and set all pixels by default to NaN.
        # Valid pixels will subsequently be overwritten with real numbers
        #
        merged = numpy.ones(shape=(size_x, size_y), dtype=numpy.float32)
        merged[:,:] = options['indef_pixelvalue']
        
        if (options['bgmode']):
            stdout_write("\r%s: Starting work on OTA %02d ...\n" % (obsid, ota))

        nonlin_data = None
        if (not options['nonlinearity'] == None):
            print "Using non-linearity coefficients from",options['nonlinearity']
            nonlin_data = podi_nonlinearity.load_nonlinearity_correction_table(options['nonlinearity'], ota)

        for cell in range(1,65):
            if (not options['bgmode']):
                stdout_write("\r%s:   OTA %02d, cell %s ..." % (obsid, ota, hdulist[cell].header['EXTNAME']))
            wm_cellx, wm_celly = hdulist[cell].header['WN_CELLX'], hdulist[cell].header['WN_CELLY']

            #
            # Special case for cell 0,7 (the one in the bottom left corner):
            # Copy the CRPIX values into the merged image header 
            #
            if (hdulist[cell].header['EXTNAME'] == "XY07"):
                # print "Setting CRPIXs", hdulist[cell].header['CRPIX1'], hdulist[cell].header['CRPIX2']
                hdu.header.update("CRPIX1", hdulist[cell].header['CRPIX1'], "Ref. pixel RA")
                hdu.header.update("CRPIX2", hdulist[cell].header['CRPIX2'], "Ref. pixel DEC")

            # Check if this is one of the broken cells
            cellmode_id = get_cellmode(hdulist[0].header, hdulist[cell].header)
            if (not cellmode_id == 0):
                # This means it either broken (id=-1) or in video-mode (id=1)
                continue

            #
            # Now extract just the data section.
            # Values are hard-coded as some of the DATASEC keywords are wrong
            #
            datasec = hdulist[cell].data[0:494, 0:480] 

            # Now overscan subtract and insert into large frame
            overscan_region = extract_region(hdulist[cell].data, '[500:530,1:494]')
            overscan_level = numpy.median(overscan_region)

            datasec -= overscan_level

            if (not options['nonlinearity'] == None):
                nonlin_correction = podi_nonlinearity.compute_cell_nonlinearity_correction(
                    datasec, wm_cellx, wm_celly, nonlin_data)
                datasec += nonlin_correction

            if (gain_correct_frames):
            # Correct for the gain variations in each cell
                try:
                    gain = float(hdulist[cell].header['GAIN'])
                    datasec *= gain
                except:
                    print "Couldn't find the GAIN header!"
                    pass


            #
            # Insert the reduced data-section of this cell into the large OTA frame
            #
            bx, tx, by, ty = cell2ota__get_target_region(wm_cellx, wm_celly) 
            merged[by:ty,bx:tx] = datasec

        # 
        # At this point we have a 4x4 Kpixel array with all cells merged
        #

        #
        # Get some information for the OTA
        #
        fppos = hdulist[0].header['FPPOS']

        try:
            filter_name = hdulist[0].header['FILTER']
        except KeyError:
            filter_name = 'UNKNOWN'
            
        try:
            exposure_time = hdulist[0].header['EXPTIME']
        except KeyError:
            exposure_time = 0

        # If we are to do some bias subtraction:
        if (not options['bias_dir'] == None):
            bias_filename = check_filename_directory(options['bias_dir'], "bias.fits")
            if (os.path.isfile(bias_filename)):
                bias = pyfits.open(bias_filename)

                # Search for the bias data for the current OTA
                for bias_ext in bias[1:]:
                    fppos_bias = bias_ext.header['FPPOS']
                    if (fppos_bias == fppos):
                        # This is the one
                        merged -= bias_ext.data
                        break

                bias.close()
                hdu.header.add_history("CC-BIAS: %s" % (os.path.abspath(bias_filename)))
                del bias
 

        # To do some dark subtraction:
        #
        # Missing here: Add treatment for frames with detectors switched on or off
        #
        if (not options['dark_dir'] == None):

            # For now assume all detectors are switched on
            detectorglow = "yes"

            dark_filename = check_filename_directory(options['dark_dir'], "dark_%s.fits" % (detectorglow))
            if (os.path.isfile(dark_filename)):
                dark = pyfits.open(dark_filename)
                darktime = dark[0].header['EXPTIME']

                # Search for the flatfield data for the current OTA
                for dark_ext in dark[1:]:
                    fppos_dark = dark_ext.header['FPPOS']
                    if (fppos_dark == fppos):
                        # This is the one
                        merged -= (dark_ext.data * exposure_time / darktime)
                        break

                dark.close()
                hdu.header.add_history("CC-DARK: %s" % (os.path.abspath(dark_filename)))
                del dark
 

        # If the third parameter points to a directory with flat-fields
        if (not options['flat_dir'] == None):
            flatfield_filename = check_filename_directory(options['flat_dir'], "flat_%s.fits" % (filter_name))
            if (os.path.isfile(flatfield_filename)):
                flatfield = pyfits.open(flatfield_filename)

                # Search for the flatfield data for the current OTA
                for ff_ext in flatfield[1:]:
                    fppos_flatfield = ff_ext.header['FPPOS']
                    if (fppos_flatfield == fppos):
                        # This is the one
                        merged /= ff_ext.data
                        break

                flatfield.close()
                hdu.header.add_history("CC-FLAT: %s" % (os.path.abspath(flatfield_filename)))
                del flatfield

        # Finally, apply bad pixel masks 
        if (not options['bpm_dir'] == None):
            # Determine which region file we need
            region_file = "%s/bpm_%s.reg" % (options['bpm_dir'], fppos)
            if (os.path.isfile(region_file)):
                # Apply the bad pixel regions to file, marking
                # all bad pixels as NaNs
                mask_broken_regions(merged, region_file)
                hdu.header.add_history("CC-BPM: %s" % (os.path.abspath(region_file)))

        #
        # If persistency correction is requested, perform it now
        #
        if (options['persistency_dir'] != None):
            if (options['verbose']): stdout_write("Applying persistency correction\n")

            # Get a list of all saturation catalogs
            full_filelist = podi_persistency.get_list_of_saturation_tables(options['persistency_dir'])
            #print full_filelist

            # Search for the saturation map for this frame and, if found,
            # mask out all saturated pixels
            saturated_thisframe = podi_persistency.select_from_saturation_tables(full_filelist, mjd, None)
            #print "this frame=",saturated_thisframe,mjd

            if (not saturated_thisframe == None):
                merged = podi_persistency.mask_saturation_defects(saturated_thisframe, ota, merged)
            hdu.header.update("SATMASK", saturated_thisframe)

            # Also pick all files within a given MJD range, and apply the 
            # appropriate correction (which, for now, is simply masking all pixels)
            filelist = podi_persistency.select_from_saturation_tables(full_filelist, mjd, [1,options['max_persistency_time']])
            if (len(filelist) > 0):
                merged = podi_persistency.correct_persistency_effects(ota, merged, mjd, filelist)
                persistency_catalog_counter = 0
                for filename in filelist:
                    persistency_catalog_counter += 1
                    keyname = "PERSIS%02d" % (persistency_catalog_counter)
                    hdu.header.update(keyname, filename)

        #
        # If requested, determine the optimal fringe scaling
        #
        fringe_scaling = None
        fringe_scaling_median, fringe_scaling_std = -1, -1
        if (options['fringe_dir'] != None):
            fringe_filename = check_filename_directory(options['fringe_dir'], "fringe__%s.fits" % (filter_name))
            fringe_vector_file = "%s/fringevectors__%s__OTA%02d.reg" % (options['fringe_vectors'], filter_name, ota)
            if (options['verbose']):
                print "fringe file:",fringe_filename, "    found:",os.path.isfile(fringe_filename)
                print "fringe vector:",fringe_vector_file, "    found:",os.path.isfile(fringe_vector_file)

            # Do not determine fringe scaling if either or the input files does 
            # not exist or any of the cells in this OTA is marked as video cell
            if (os.path.isfile(fringe_filename)
                and os.path.isfile(fringe_vector_file)
                and hdu.header['CELLMODE'].find("V") < 0
                ):
                fringe_hdu = pyfits.open(fringe_filename)
                for ext in fringe_hdu[1:]:
                    if (extname == ext.header['EXTNAME']):
                        if (options['verbose']): print "Working on fringe scaling for",extname
                        fringe_scaling = podi_fringing.get_fringe_scaling(merged, ext.data, fringe_vector_file)
                        if (not fringe_scaling == None):
                            good_scalings = three_sigma_clip(fringe_scaling[:,6], [0, 1e9])
                            fringe_scaling_median = numpy.median(good_scalings)
                            fringe_scaling_std    = numpy.std(good_scalings)
                        break
                hdu.header.add_history("fringe map: %s" % fringe_filename)
                hdu.header.add_history("fringe vector: %s" % fringe_vector_file)
                fringe_hdu.close()
            #print "FRNG_SCL", fringe_scaling_median
            #print "FRNG_STD", fringe_scaling_std
            hdu.header.update("FRNG_SCL", fringe_scaling_median)
            hdu.header.update("FRNG_STD", fringe_scaling_std)
            hdu.header.update("FRNG_OK", (fringe_scaling != None))

        # Insert the DETSEC header so IRAF understands where to put the extensions
        start_x = ota_c_x * 4096
        start_y = ota_c_y * 4096        
        end_x = start_x + det_x2 - det_x1
        end_y = start_y + det_y2 - det_y1
        detsec_str = "[%d:%d,%d:%d]" % (start_x, end_x, start_y, end_y)
        hdu.header.update("DETSEC", detsec_str, "position of OTA in focal plane")
                
        if (cmdline_arg_isset("-simplewcs") or options['wcs_distortion'] != None):
            #
            # Fudge with the WCS headers, largely undoing what's in the fits file right now,
            # and replacing it with a simpler version that hopefully works better
            #
            hdu.header['CTYPE1'] = "RA---TAN"
            hdu.header['CTYPE2'] = "DEC--TAN"
            del hdu.header['WAT0_001']
            del hdu.header['WAT1_001']
            del hdu.header['WAT1_002']
            del hdu.header['WAT1_003']
            del hdu.header['WAT1_004']
            del hdu.header['WAT1_005']
            del hdu.header['WAT2_001']
            del hdu.header['WAT2_002']
            del hdu.header['WAT2_003']
            del hdu.header['WAT2_004']
            del hdu.header['WAT2_005']
        # in any case, add the CUNIT headers that are missing by default
        hdu.header.update("CUNIT1", "deg", "")
        hdu.header.update("CUNIT2", "deg", "")

        coord_j2000 = ephem.Equatorial(hdu.header['RA'], hdu.header['DEC'], epoch=ephem.J2000)
        if (not options['target_coords'] == None):
            if (options['verbose']): print "Overwriting target positions",options['target_coords']
            ra, dec = options['target_coords']
            coord_j2000 = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
            if (options['verbose']): print "\n\nAdjusting ra/dec:", coord_j2000.ra, coord_j2000.dec,"\n\n"

        # Write the CRVALs with the pointing information
        #print numpy.degrees(coord_j2000.ra), numpy.degrees(coord_j2000.dec)  
        hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  
        hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) 

        # Compute total offsets as the sum from pointing and dither offset
        # Now add the pointing and dither offsets
        #print offset_total[0] / 3600. / numpy.cos(numpy.radians(hdu.header['CRVAL2']))
        #print hdu.header['CRVAL2'], numpy.cos(numpy.radians(hdu.header['CRVAL2']))
        #hdu.header['CRVAL1'] += offset_total[0] / 3600. / numpy.cos(numpy.radians(hdu.header['CRVAL2']))
        #hdu.header['CRVAL2'] += offset_total[1] / 3600.
        #
        # To do:
        # =========================================================
        # Check if the above still makes sense !!!!
        # In particular the addition of the telescope offsets 
        # should be included in RA/DEC already !!!
        # =========================================================
        #
        if (options['offset_pointing'] != None or options['offset_dither'] != None):
            offset_total = numpy.array([0,0])
            if (options['offset_pointing'] != None):
                offset_total += numpy.array(options['offset_pointing'])
            if (options['offset_dither'] != None):
                offset_total += numpy.array(options['offset_dither'])
            if (options['verbose']): 
                stdout_write("Adding user's WCS offset (ra: %f, dec: %f degrees)\n" % (offset_total[0], offset_total[1]))
            hdu.header['CRVAL1'] += offset_total[0] 
            hdu.header['CRVAL2'] += offset_total[1]


        # Now add the canned WCS solution
        if (options['wcs_distortion'] != None):
            #print "Adding header from WCS minifits (%s)" % (extname)
            wcs = pyfits.open(options['wcs_distortion'])
            wcs_header = wcs[extname].header

            cards = wcs_header.ascardlist()
            for card in cards:
                if (card.key == 'CRVAL1'):
                    hdu.header["CRVAL1"] -= wcs_header['CRVAL1'] / math.cos(math.radians(hdu.header['CRVAL2']))
                elif (card.key == "CRVAL2"):
                    hdu.header['CRVAL2'] -= wcs_header['CRVAL2']
                else:
                    hdu.header.update(card.key, card.value, card.comment)

            # Make sure to write RAs that are positive
            if (hdu.header["CRVAL1"] < 0):
                hdu.header['CRVAL1'] += 360.
                
        # Insert the new image data. This also makes sure that the headers
        # NAXIS, NAXIS1, NAXIS2 are set correctly
        hdu.data = merged

        source_cat = None
        if (options['fixwcs']):
            # Create source catalog
            
            if (hdu.header['CELLMODE'].find("V") > -1):
                source_cat = None
            elif (False):
                source_cat = podi_findstars.find_stars(hdu, binning=4, boxsize=24, dumpfile=None, verbose=False,
                                                   detect_threshold=1., detect_minarea=4, roundness_limit=[-0.2,+0.2],
                                                   max_starcount=150,
                                                   extension_id=ota)
            else:
                hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdu.header, data=hdu.data)])
                obsid = hdulist[0].header['OBSID']
                process_id = os.getpid()
                fitsfile = "%s/tmp.pid%d.%s_OTA%02d.fits" % (sitesetup.scratch_dir, process_id, obsid, ota)
                catfile = "%s/tmp.pid%d.%s_OTA%02d.cat" % (sitesetup.scratch_dir, process_id, obsid, ota)
                hdulist.writeto(fitsfile, clobber=True)
                full_path = os.path.abspath(sys.argv[0])
                basepath, dummy = os.path.split(full_path)
                sex_config_file = "%s/.config/wcsfix.sex" % (basepath)
                parameters_file = "%s/.config/wcsfix.sexparam" % (basepath)
                sexcmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s >&/dev/null" % (
                    sex_config_file, parameters_file, catfile, fitsfile)
                if (options['verbose']): print sexcmd
                os.system(sexcmd)
                try:
                    try:
                        source_cat = numpy.loadtxt(catfile)
                    except IOError:
                        print "The Sextractor catalog is empty, ignoring this OTA"
                        source_cat = None
                    else:
                        if (source_cat.shape[0] == 0 or source_cat.ndim < 2):
                            source_cat == None
                        else:
                            source_cat[:,8] = ota
                            flags = source_cat[:,7]
                            valid = (flags != 3)
                            if (verbose): print "source-cat:", source_cat.shape, numpy.sum(valid)
                            source_cat = source_cat[valid]
                except:
                    source_cat == None
                clobberfile(fitsfile)
                clobberfile(catfile)

            fixwcs_data = None
            if (source_cat != None):
                if (source_cat.shape[0] > 0 and source_cat.ndim == 2):
                    odi_ra = source_cat[:,0]
                    odi_dec = source_cat[:,1]
                    odi_mag = source_cat[:,14] #-2.5 * numpy.log10(source_cat[:,6]) + 30

                    # Read the reference catalog
                    center_ra, center_dec = center_coords(hdu.header)
                    search_size = (8+6) * (1./60.)
                    ipp_cat = podi_search_ipprefcat.get_reference_catalog(center_ra, center_dec, search_size, 
                                                                          basedir=sitesetup.wcs_ref_dir,
                                                                          cattype=sitesetup.wcs_ref_type)
                    ref_ra = ipp_cat[:,0]
                    ref_dec = ipp_cat[:,1]
                    ref_mag = ipp_cat[:,3]

                    # Cut down the number of stars to < 100 to save computing time
                    ota_odi_ra, ota_odi_dec, ota_odi_mag = podi_fixwcs.pick_brightest(odi_ra, odi_dec, odi_mag, 50)
                    ota_ref_ra, ota_ref_dec, ota_ref_mag = podi_fixwcs.pick_brightest(ref_ra, ref_dec, ref_mag, 50)

                    if (verbose): print "sending %s and %d to shift_align_wcs" % (ota_odi_ra.shape[0], ota_ref_ra.shape[0])
                    dx, dy, n, matchcount = podi_fixwcs.shift_align_wcs(ota_odi_ra, ota_odi_dec, ota_ref_ra, ota_ref_dec)
                    if (verbose): print "WCSFIX dx/dy =", dx, dy
                    fixwcs_data = (odi_ra, odi_dec, ref_ra, ref_dec, dx, dy, source_cat, ipp_cat, matchcount)
        else:
            fixwcs_data = None

        #
        # Sample that background at random place so we can derive a median background level later on
        # This is necessary both for the pupil ghost correction AND the fringing correction
        #
        starcat = None
        if (fixwcs_data != None):
            src_cat = fixwcs_data[6]
            ota_x, ota_y = src_cat[:,2], src_cat[:,3]
            # print "ota_x=",ota_x
            # print "ota_y=",ota_y
            starcat = (ota_x, ota_y)
        # Now sample the background, excluding regions close to known sources
        sky_samples = numpy.array(podi_fitskybackground.sample_background(data=merged, wcs=None, 
                                                                         starcat=starcat, 
                                                                         min_found=200, boxwidth=30))

        if (sky_samples.shape[0] <= 0): #not sky_samples.shape[1] > 4):
            print "Something went wrong with the sky-calculation"
            print "sky-samples =",sky_samples
            print "sky-samples.shape =",sky_samples.shape
            sky_level_median, sky_level_mean, sky_level_std = -1, -1, -1
            sky_samples = None
        else:
            sky_level_median = numpy.median(sky_samples[:,4])
            sky_level_mean   = numpy.mean(sky_samples[:,4])
            sky_level_std    = numpy.std(sky_samples[:,4])
            hdu.header.update("SKY_MEDI", sky_level_median, "sky-level median")
            hdu.header.update("SKY_MEAN", sky_level_mean, "sky-level mean")
            hdu.header.update("SKY_STD", sky_level_std, "sky-level rms")

    data_products['hdu'] = hdu
    data_products['wcsdata'] = fixwcs_data
    data_products['sky-samples'] = sky_samples
    data_products['sky'] = (sky_level_median, sky_level_mean, sky_level_std)
    data_products['fringe_scaling'] = fringe_scaling
    data_products['sourcecat'] = source_cat

    return data_products #hdu, fixwcs_data
    


#########
#
# This routine is a wrapper around the actual collect_reduce_ota routine,
# mainly dealing with the parallelization and inter-process communication
#
#########
def parallel_collect_reduce_ota(queue, return_queue,
                                options=None):

    while (True):
        cmd_quit, filename, ota_id = queue.get()
        if (cmd_quit):
            queue.task_done()
            return

        # Do the work
        data_products = collect_reduce_ota(filename, options=options)

        # Add the results to the return_queue so the master process can assemble the result file
        # print "Adding results for OTA",ota_id,"to return queue"
        # return_queue.put( (hdu, ota_id, wcsfix_data) )
        return_queue.put( (ota_id, data_products) )
        queue.task_done()
        
    return


#########
#
# collectcells:
# Handles all filename operations, ensuring all required files exist, and hands the work on
# each OTA off to the suite of worker processes. Finally assembles all results and writes the output-file.
#
#########
def collectcells(input, outputfile,
                 batchmode=False,
                 verbose=False,
                 options=None,
                 showsplash=True):

    if (showsplash):
        stdout_write("""\

    **********************************************************************
    * This is podi_collectcells                                          *
    * (c) 2012-2013: Ralf Kotulla, kotulla@uwm.edu                       *
    *                University of Wisconsin, Milwaukee                  *
    *                                                                    *
    * Please acknowledge the author when using any products generated    *
    * with this tool. For comments, questions or ideas for improvement   *
    * please send an email to kotulla@uwm.edu. Thank you!                *
    **********************************************************************

""")

    # print "Received options:", options
    if (options == None): options = set_default_options()

    if (options['verbose']):
        stdout_write("\nThese are the options we are using:\n")
        for opt_name, opt_value in options.iteritems():
            stdout_write("   %30s: %s\n" % (opt_name, opt_value))
        stdout_write("----- end options\n\n")

    # afw = podi_asyncfitswrite.async_fits_writer(1)

    if (os.path.isfile(input)):
        # Assume this is one of the fits files in the right directory
        # In that case, extract the FILENAME header and convert it into 
        # the filebase we need to construct the filenames of all OTA fits files.
        hdulist = pyfits.open(input)
        filebase = hdulist[0].header['FILENAME'][:18]
        hdulist.close()
        del hdulist

        # Split the input filename to extract the directory part
        directory, dummy = os.path.split(input)

    elif (os.path.isdir(input)):
        # As a safety precaution, if the first parameter is the directory containing 
        # the files, extract just the ID string to be used for this script
        if (input[-1] == "/"):
            input = input[:-1]

        basedir, filebase = os.path.split(input)
        directory = input

    else:
        stdout_write("Unable to open file %s, aborting!\n" % input)
        return

    #print "Merging cells for frame %s" % (basename)

    if (outputfile == None):
        outputfile = "%s/%s.fits" % (directory, filebase)

    hdulist = None
    for ota in all_otas:
        filename = "%s/%s.%02d.fits" % (directory, filebase, ota)
        if (not os.path.isfile(filename)):
            filename = "%s/%s.%02d.fits.fz" % (directory, filebase, ota)
        try:
            hdulist = pyfits.open(filename)
            break
        except:
            #stdout_write("Problem opening an existing fits-file (%s), aborting!\n" % filename)
            pass
            continue

    if (hdulist == None):
        print "Something is wrong here, can't find/open any of the files..."
        return -1

    if (outputfile.find("%") >= 0):
        # The output filename contains special tags that should 
        # be replaced by values from the file header

        header = hdulist[0].header
        while (outputfile.find("%") >= 0):
            start = outputfile.find("%") 
            if (outputfile[start:start+7] == "%FILTER"):
                outputfile = outputfile[:start] + header['FILTER'] + outputfile[start+7:]
            elif (outputfile[start:start+7] == "%OBJECT"):
                # The object name might contain spaces, replace them with underscores
                objectname = header['OBJECT'].replace(' ', '_')
                objectname = objectname.replace(',', '_')
                outputfile = outputfile[:start] + objectname  + outputfile[start+7:]
            elif (outputfile[start:start+6] == "%OBSID"):
                outputfile = outputfile[:start] + header['OBSID'] + outputfile[start+6:]
            elif (outputfile[start:start+7] == "%OBSSEQ"):
                obsid = header['OBSID']
                dot_position = obsid.find(".")
                obs_sequence = obsid[:dot_position]
                outputfile = outputfile[:start] + obs_sequence + outputfile[start+7:]
            elif (outputfile[start:start+8] == "%EXPTIME"):
                outputfile = "%s%.1f%s" % (outputfile[:start], header['EXPTIME'], outputfile[start+8:])
            else:
                stdout_write("found unknown tag in %s\n" % outputfile)
                break

        del header

        stdout_write("Replaced some keywords, new output filename: ---> %s\n" % (outputfile))

    # We know enough about the current frame, so close the file
    hdulist.close()
    del hdulist

    # Check if the output file contains a new directory. 
    chk_directory, chk_filename = os.path.split(outputfile)
    changed_outputfilename = False
    if (chk_directory != ''):
        # The output file does contain a directory
        if (not os.path.isdir(chk_directory)):
            # if the directory doesn't exist ey, create it
            stdout_write("Output directory does not exist, creating it...\n")
            os.mkdir(chk_directory)
            changed_outputfilename = True
    if (chk_filename == ''):
        # This happens if the output is just a directory
        stdout_write("Output filename is a directory, adding default filename (collectcells.fits)\n")
        outputfile += "collectcells.fits"
        changed_outputfilename = True
    if (outputfile[-5:] != ".fits"):
        # Output filenames have to end with .fits
        stdout_write("no fits extension given to filename, adding .fits\n")
        outputfile += ".fits"
        changed_outputfilename = True
    if (changed_outputfilename):
        stdout_write("After revision, output filename now is %s\n" % (outputfile))

    if (os.path.isfile(outputfile) and not options['clobber']):
        print "#####################################################"
        print "#"
        print "# File %s already exists, skipping!" % (outputfile)
        print "#"
        print "#####################################################"
        print "\n"
        return

    #
    # Start assembling the new list of HDUs
    #
    list_of_otas_to_collect = available_ota_coords
    if (options['central_only']):
        list_of_otas_to_collect = central_array_ota_coords

    ota_list = [None] * (len(list_of_otas_to_collect)+1)
    # And add the primary HDU to make the fits file a valid one
    ota_list[0] = pyfits.PrimaryHDU()
    
    #
    # Set up the parallel processing environment
    #
    #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
    queue = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.Queue()
    
    processes = []

    worker_args = (queue, return_queue, options)

    number_extensions = 0

    list_of_otas_being_reduced = []
    for ota_id in range(len(list_of_otas_to_collect)):
        ota_c_x, ota_c_y = list_of_otas_to_collect[ota_id]        
        ota = ota_c_x * 10 + ota_c_y

        if (cmdline_arg_isset('-singleota')):
            single_ota = int(get_cmdline_arg("-singleota"))
            if (ota != single_ota):
                continue

        filename = "%s/%s.%02d.fits" % (directory, filebase, ota)
        if (not os.path.isfile(filename)):
            # Check if there's a .fz compressed file
            filename = "%s/%s.%02d.fits.fz" % (directory, filebase, ota)
            if (not os.path.isfile(filename)):
                # This file doesn't seem to exist, so don't ask the worker to do anything
                continue

        #print "Commanding work for extension",ota
        list_of_otas_being_reduced.append(list_of_otas_to_collect[ota_id])

        queue.put( (False, filename, ota_id+1) )

    # Create all processes to handle the actual reduction and combination
    #print "Creating",number_cpus,"worker processes"
    if ('profile' in options and number_cpus == 1):
        # 
        # If profiling is activated, run one only one processor and in non-multiprocessing mode
        #
        # Tell all workers to shut down when no more data is left to work on
        #for i in range(len(processes)):
        if (verbose): stdout_write("Sending quit command!\n")
        queue.put((True,None,None))

        while (True):
            print "Doing work single-processed"
            cmd_quit, filename, ota_id = queue.get()
            if (cmd_quit):
                queue.task_done()
                break

            # Do the work
            data_products = collect_reduce_ota(filename, options=options)

            # Add the results to the return_queue so the master process can assemble the result file
            # print "Adding results for OTA",ota_id,"to return queue"
            # return_queue.put( (hdu, ota_id, wcsfix_data) )
            return_queue.put( (ota_id, data_products) )
            queue.task_done()
    else:
        for i in range(number_cpus):
            p = multiprocessing.Process(target=parallel_collect_reduce_ota, args=worker_args)
            p.start()
            processes.append(p)
            time.sleep(0.01)

        # Tell all workers to shut down when no more data is left to work on
        for i in range(len(processes)):
            if (verbose): stdout_write("Sending quit command!\n")
            queue.put((True,None,None))


    #
    # By now all workers have computed their HDUs or are busy doing so,
    # let's extract their results from the return queue and assemble the ota list
    #
    fixwcs_odi_ra, fixwcs_odi_dec = numpy.array([]), numpy.array([])
    fixwcs_ref_ra, fixwcs_ref_dec = numpy.array([]), numpy.array([])
    fixwcs_odi_x,  fixwcs_odi_y   = numpy.array([]), numpy.array([])
    fixwcs_extension = numpy.array([])
    reference_catalog = None

    fixwcs_odi_sourcecat = None
    fixwcs_bestguess = numpy.ones(shape=(len(list_of_otas_being_reduced),2)) * -1000

    ############################################################
    #
    # In this loop:
    # - Sort all indidually reduced OTAs into the larger frame, 
    #   putting them all in the right order (spiraling from inside-out)
    # - Combine all the WCS fitting information from all OTAs
    #
    ############################################################
    global_number_matches = None
    sky_samples = {}
    fringe_scaling = None
    global_source_cat = None
    for i in range(len(list_of_otas_being_reduced)):
        #hdu, ota_id, wcsfix_data = return_queue.get()
        ota_id, data_products = return_queue.get()

        hdu = data_products['hdu']
        wcsfix_data = data_products['wcsdata']

        #if ("persistency_map_updated" in data_products):
        #    # We also received an updated persistency map
        #    persistency[ota_id] = data_products["persistency_map_updated"]

        if (hdu == None):
            continue

        ota_list[ota_id] = hdu
        
        sky_samples[hdu.header['EXTNAME']] = data_products['sky-samples']

        if (options['fixwcs'] and not options['update_persistency_only']):
            
            if (wcsfix_data != None):
                odi_ra, odi_dec, ref_ra, ref_dec, dx, dy, source_cat, reference_cat, number_matches = wcsfix_data

                # Append all the returned ODI and reference catalog data to the 
                # previous list to make one big catalog across all OTAs.
                fixwcs_odi_ra  = numpy.append(fixwcs_odi_ra,  odi_ra,  axis=0)
                fixwcs_odi_dec = numpy.append(fixwcs_odi_dec, odi_dec, axis=0)
                fixwcs_ref_ra  = numpy.append(fixwcs_ref_ra,  ref_ra,  axis=0)
                fixwcs_ref_dec = numpy.append(fixwcs_ref_dec, ref_dec, axis=0)

                #print "Adding data to long list"
                fixwcs_odi_x = numpy.append(fixwcs_odi_x, source_cat[:,2], axis=0)
                fixwcs_odi_y = numpy.append(fixwcs_odi_y, source_cat[:,3], axis=0)
                _ota = numpy.ones(shape=(source_cat.shape[0],1)) * ota_id
                fixwcs_extension = numpy.append(fixwcs_extension, _ota[:,0], axis=0)
                #print "source_cat[:,2]=",source_cat[:,2]
                #print "source_cat[:,3]=",source_cat[:,3]

                reference_catalog = reference_cat if (reference_catalog == None) else numpy.append(reference_catalog, reference_cat, axis=0)
                #RK fixwcs_odi_sourcecat = source_cat if (fixwcs_odi_sourcecat == None) else numpy.append(fixwcs_odi_sourcecat, source_cat, axis=0)

                fixwcs_bestguess[i,:] = [dx, dy]

                if (number_matches != None):
                    # Deal with the individual matches from this OTA
                    #print hdu.header['EXTNAME'], number_matches.shape
                    tmp = numpy.ones(shape=(number_matches.shape[0], number_matches.shape[1]+1))
                    tmp[:,:-1] = number_matches
                    tmp[:,-1] *= ota_id
                    if (global_number_matches == None):
                        global_number_matches = tmp
                    else:
                        global_number_matches = numpy.append(global_number_matches, tmp, axis=0)
                #print "\n\n\n#matches was",number_matches.shape,"now is ",global_number_matches.shape,"\n\n\n"

            #add_to_bestguess = numpy.array([dx, dy]).reshape((1,2))
            #print fixwcs_bestguess.shape, add_to_bestguess.shape
            #continue
            #fixwcs_bestguess = numpy.append(fixwcs_bestguess, add_to_bestguess, axis=0)

        if (not data_products['fringe_scaling'] == None):
            #print data_products['fringe_scaling'].shape
            #if (fringe_scaling != None): print fringe_scaling.shape
            fringe_scaling = data_products['fringe_scaling'] if fringe_scaling == None else \
                numpy.append(fringe_scaling, data_products['fringe_scaling'], axis=0)

        if (not data_products['sourcecat'] == None):
            global_source_cat = data_products['sourcecat'] if (global_source_cat == None) else numpy.append(global_source_cat, data_products['sourcecat'], axis=0)

    if (options['fixwcs'] and not options['update_persistency_only']):
        x = open("fixwcs.nmatches","w")
        numpy.savetxt(x, global_number_matches)
        x.close()
        del x

    # Now all processes have returned their results, terminate them 
    # and delete all instances to free up memory
    for cur_process in processes:
        cur_process.terminate()
        del cur_process

    stdout_write(" done!\n")

    if (options['fixwcs'] and verbose):
        print fixwcs_extension
        print fixwcs_odi_x
        print fixwcs_odi_y
        print fixwcs_bestguess.shape
        print fixwcs_bestguess
    
    #
    # Now do some post-processing:
    # 1) Add or overwrite some headers with values from an external wcs minifits file
    #    to improve the wcs accuracy - remover - no done in parallel during cell collection.
    # 2) Move a couple of headers out of each individual extension and put it in the 
    #    primary extension instead (defined in headers_to_inherit, see podi_definitions)
    # 3) Delete a bunch of headers that are no longer necessary (defined in 
    #    headers_to_delete_from_otas, see podi_definitions)
    # 4) Collect all WCS information and fix pointing and rotation errors. Create some
    #    diagnostic plots
    # 5) Write the updated persistency map to file - deleted, no longer required
    # 6) compute the global sky level
    # 7) Determine the best-fit pupil ghost scaling and remove the pupil ghost contribution
    # 8) Compute global fringe scaling and remove scaled fringe template
    #

    # First step:
    # Delete all HDU entries that are set to None, indicating that there was something
    # seriously wrong with them
    #print "Post-processing start w/",len(ota_list),"extensions"
    tmp_ota_list = []
    for ext in range(len(ota_list)):
        if (ota_list[ext] != None):
            tmp_ota_list.append(ota_list[ext])
    ota_list = tmp_ota_list
    #print "cleaned up, now",len(ota_list),"extensions left"

    # Now update the headers in all OTA extensions.
    for extension in range(1, len(ota_list)):
        ota = ota_list[extension]

        if (ota == None):
            continue

        if (cmdline_arg_isset("-prep4sex")):
            continue

        for header in headers_to_inherit:
            # Make sure the header we are about to move exists in the first place
            if (not header in ota.header):
                continue

            # Check if the header already exists in the primary header. If not add it!
            if (not header in ota_list[0].header):
                card = ota.header.ascardlist()[header]
                ota_list[0].header.update(card.key, card.value, card.comment)
                #value = ota.header[header]
                #ota_list[0].header.update(header, value, "DESCRIPTION")
            
            # By now the value should exist in the primary header, 
            # so delete it from each of the extensions
            del ota.header[header]
                
        # Set the inherit keyword so that the headers removed from each 
        # extension are instead inherited from the primary
        ota.header.update("INHERIT", True, "Inherit headers from PrimaryHDU")

        for header in headers_to_delete_from_otas:
            # As above, make sure header exists
            if (not header in ota.header):
                continue
            del ota.header[header]

    # 
    # Now combine all sky-samples to compute the global background level.
    # Take care to exclude OTAs marked as guide-OTAs and those not covered 
    # by the narrow-band filters.
    # 
    sky_samples_global = None #numpy.empty(0)
    valid_ext = otas_for_photometry[get_valid_filter_name(ota_list[0].header)]
    for ext in sky_samples:
        # print ext, valid_ext, int(ext[3:5])
        ota_number = int(ext[3:5])
        ota_name = "OTA%02d.SCI" % (ota_number)
        not_a_guiding_ota = False
        for i in range(1, len(ota_list)):
            if (ota_list[i].header['EXTNAME'] == ota_name):
                not_a_guiding_ota = (ota_list[i].header['CELLMODE'].find("V") < 0)
                break
        if (ota_number in valid_ext and not_a_guiding_ota and not sky_samples[ext] == None):
            if (sky_samples_global == None):
                sky_samples_global = sky_samples[ext]
            else:
                sky_samples_global = numpy.append(sky_samples_global, sky_samples[ext], axis=0)

    sky_global_median = numpy.median(sky_samples_global[:,4])
    ota_list[0].header.update("SKYLEVEL", sky_global_median)

    #
    # Now that we have the global sky-level, subtract the 
    # contribution of the pupil ghost to the science frame.
    # Different from the case of the calibration frames, use the radial
    # profile here, and ignore rotation (not needed anyway) to speed things up.
    #
    if (options['pupilghost_dir'] != None):
        filter_level = get_filter_level(ota_list[0].header)
        filter_name = get_valid_filter_name(ota_list[0].header)
        pg_template = "%s/pupilghost_radial___level_%d.fits" % (options['pupilghost_dir'], filter_level)
        stdout_write("looking for radial pupil ghost template %s...\n" % (pg_template))
        # If we have a template for this level
        if (os.path.isfile(pg_template)):
            stdout_write("\n   Using pupilghost template %s, filter %s ... " % (pg_template, filter_name))
            pg_hdu = pyfits.open(pg_template)

            # Find the optimal scaling factor
            scaling, scaling_std = podi_matchpupilghost.get_pupilghost_scaling(ota_list, pg_hdu)

            # And subtract the scaled pupilghost templates.
            podi_matchpupilghost.subtract_pupilghost(ota_list, pg_hdu, scaling, rotate=False)

            ota_list[0].header.update("PUPLGOST", pg_template, "p.g. template")
            ota_list[0].header.update("PUPLGFAC", scaling, "pupilghost scaling")
            bg_scaled = podi_matchpupilghost.scaling_factors[filter_name]*sky_global_median
            ota_list[0].header.update("PUPLGFA2", bg_scaled, "analytical pupilghost scaling")
            stdout_write(" done!\n")

            pg_hdu.close()

    #
    # Fix the WCS if requested
    #
    if (options['fixwcs'] and False):

        scaling_xxx = 1800.
        # New method using external match program
        source_cat_file = outputfile[:-5]+".wcs.src.cat"
        file = open(source_cat_file, "w")
        # extract only the ra,dec,magnitude columns
        cat_src = numpy.empty(shape=(fixwcs_odi_sourcecat.shape[0],3))
        cat_src[:,0:2] = fixwcs_odi_sourcecat[:,0:2]*scaling_xxx
        cat_src[:,2] = fixwcs_odi_sourcecat[:,13]
        numpy.savetxt(file, cat_src, fmt="%.6f %.6f %.3f")
        file.close()

        reference_cat_file = outputfile[:-5]+".wcs.2mass.cat"
        file = open(reference_cat_file, "w")
        cat_ref = numpy.empty(shape=(reference_catalog.shape[0],3))
        cat_ref[:,0:2] = reference_catalog[:,0:2] * scaling_xxx
        cat_ref[:,2] = reference_catalog[:,2]
        numpy.savetxt(file, cat_ref, fmt="%.6f %.6f %.3f")
        file.close()

        import matchcat_externalmatch as mc
        #transf = mc.run_match(source_cat_file, reference_cat_file)
        transf = mc.run_match(reference_cat_file, source_cat_file, "identity rotangle=0 rottol=5 match=1")
        a,b,c,d,e,f = transf

        #Apply correction to WCS headers
        for extension in range(1, len(ota_list)):
            if (ota_list[extension] == None):
                continue

            try:
                cd11 = ota_list[extension].header['CD1_1']
                cd12 = ota_list[extension].header['CD1_2']
                cd21 = ota_list[extension].header['CD2_1']
                cd22 = ota_list[extension].header['CD2_2']

                ota_list[extension].header['CD1_1'] = b*cd11 + e*cd12
                ota_list[extension].header['CD1_2'] = c*cd11 + f*cd12
                ota_list[extension].header['CD2_1'] = b*cd21 + e*cd22
                ota_list[extension].header['CD2_2'] = c*cd21 + f*cdf2

                ota_list[extension].header['CRVAL1'] += (a/scaling_xxx)
                ota_list[extension].header['CRVAL2'] += (d/scaling_xxx)
            except:
                pass

    if (options['fixwcs']):
        debuglog = outputfile+".wcsdebug"
        declination = ota_list[1].header['CRVAL2']
        #
        # Combine the most-matches shifts from all available OTAs, and
        # determine a rough first shift position.
        #
        numpy.savetxt("number_matches", global_number_matches)
        best_guess, contrast, drxy = podi_fixwcs.optimize_shift(global_number_matches, 
                                                                declination=declination,
                                                                debuglogfile=debuglog)
        stdout_write("Found offset: %.2f', %.2f' (+/- %.1f''), contrast %.1f (%d)\n" % (
                best_guess[0]*60., best_guess[1]*60., drxy[0]*3600*0.5, contrast, best_guess[2]))

        wcs_shift = best_guess[0:2]

        # Now pickle some data for further development
        if (podi_fixwcs_rotation.output_debug_catalogs):
            numpy.savetxt("numsave.fixwcs_ref_ra.txt", fixwcs_ref_ra)
            numpy.savetxt("numsave.fixwcs_ref_dec.txt", fixwcs_ref_dec)
            numpy.savetxt("numsave.fixwcs_odi_ra.txt", fixwcs_odi_ra)
            numpy.savetxt("numsave.fixwcs_odi_dec.txt", fixwcs_odi_dec)
            numpy.savetxt("numsave.fixwcs_odi_y.txt", fixwcs_odi_y)
            numpy.savetxt("numsave.fixwcs_odi_x.txt", fixwcs_odi_x)
            numpy.savetxt("numsave.wcs_shift_guess.txt", wcs_shift_guess)
            numpy.savetxt("numsave.wcs_shift_refinement.txt", wcs_shift_refinement)

        fixrot_trans = podi_fixwcs_rotation.improve_match_and_rotation(
            fixwcs_ref_ra, fixwcs_ref_dec,
            fixwcs_odi_ra, fixwcs_odi_dec,
            wcs_shift,
            matching_radius=[5,2,2], n_repeats=3,
            verbose=False)

        # Now apply all shifts and rotations to the ODI source catalog.
        # catalog is fixwcs_odi_sourcecat
        # columns are: 0/1: ra/dec
        #              2/3: x/y
        #RK odi_sourcecat_modified = fixwcs_odi_sourcecat.copy()
        #RK if (verbose): print "wcs-shift=",wcs_shift
        #RK odi_sourcecat_modified[:,0:2] += wcs_shift
#       #RK  odi_sourcecat_modified[:,1] -= wcs_shift[1]
        #RK odi_sourcecat_modified[:,0:2] = podi_fixwcs_rotation.apply_transformation(fixrot_trans, odi_sourcecat_modified[:,0:2])

        raw_radec = global_source_cat[:,0:2] + wcs_shift
        global_source_cat[:,0:2] = podi_fixwcs_rotation.apply_transformation(fixrot_trans, raw_radec)
        
        # Now we have the corrected catalog, match again with the full 2mass reference catalog
        # 2mass catalog in variable reference_catalog
        #RK odi_2mass_matched = podi_matchcatalogs.match_catalogs(reference_catalog[:,0:2], odi_sourcecat_modified)
        odi_2mass_matched = podi_matchcatalogs.match_catalogs(reference_catalog[:,0:2], global_source_cat)

        count = numpy.sum(odi_2mass_matched[:,2] > 0)
        print "Found ",count," matched odi+2mass pairs"

        numpy.savetxt("odi+2mass.matched", odi_2mass_matched)

        ota_outlines = derive_ota_outlines(ota_list)
        if (True):
            # print "Creating some diagnostic plots"
            import podi_diagnosticplots
            podi_diagnosticplots.wcsdiag_scatter(odi_2mass_matched, outputfile[:-5]+".wcs1", 
                                                 options=options)
            podi_diagnosticplots.wcsdiag_shift(odi_2mass_matched, outputfile[:-5]+".wcs2", 
                                                 options=options, ota_outlines=ota_outlines)
        
        source_cat_file = outputfile+".src.cat"
        file = open(source_cat_file, "w")
        #RK numpy.savetxt(file, fixwcs_odi_sourcecat)
        numpy.savetxt(file, global_source_cat)
        file.close()

        reference_cat_file = outputfile+".2mass.cat"
        file = open(reference_cat_file, "w")
        numpy.savetxt(file, reference_catalog)
        file.close()

        checkalign = outputfile+".cadat"
        x = open(checkalign, "w")
        dummy = numpy.ones(shape=(fixwcs_ref_ra.shape[0],2))
        dummy[:,0] = fixwcs_ref_ra[:]
        dummy[:,1] = fixwcs_ref_dec[:]
        numpy.savetxt(x, dummy)

        y = open("2mass.cat", "w")
        numpy.savetxt(y, dummy)
        y.close()

        print >>x, "\n\n\n\n\n"
        dummy = numpy.ones(shape=(fixwcs_odi_ra.shape[0],2))
        dummy[:,0] = fixwcs_odi_ra[:]
        dummy[:,1] = fixwcs_odi_dec[:]
        numpy.savetxt(x, dummy)


        #
        # Apply all WCS shifts to the WCS information in the ouput file header
        #
        #print "\n\n\n\nbefore any shifts =",ota_list[1].header['CRVAL1'], ota_list[1].header['CRVAL2'],"\n\n\n"
        for extension in range(1, len(ota_list)):
            if (ota_list[extension] == None):
                continue
            podi_fixwcs.apply_wcs_shift(wcs_shift, ota_list[extension].header,
                                        fixrot_trans)

            
        #
        # Save the two catalogs to the output file
        #

        # Save a simple copy of the 2MASS reference catalog
        ota_list.append( sexcat_to_tablehdu(fixwcs_ref_ra, fixwcs_ref_dec) )

        # Re-compute all ODI star positions from their pixel positions to match the new WCS solution.
        # Also return the catalog as TableHDU so we can add it to the output file
        #RK odi_cat_hdu, odi_source_catalog = odi_sources_to_tablehdu(ota_list, fixwcs_odi_sourcecat)
        odi_cat_hdu, odi_source_catalog = odi_sources_to_tablehdu(ota_list, global_source_cat)
        ota_list.append(odi_cat_hdu)

        print >>x, "\n\n\n\n\n"
        dummy = numpy.ones(shape=(fixwcs_odi_ra.shape[0],2))
        dummy[:,0] = fixwcs_odi_ra[:] + wcs_shift[0]
        dummy[:,1] = fixwcs_odi_dec[:] + wcs_shift[1]
        numpy.savetxt(x, dummy)

        print >>x, "\n\n\n\n\n"
        numpy.savetxt(x, odi_2mass_matched)

        #print >>x, "\n\n\n\n\n"
        #numpy.savetxt(x, fixwcs_odi_sourcecat)

        x.close()

        numpy.savetxt("matched_cat.cat", odi_2mass_matched)

        zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = 99., 99., None, 99.
        if (options['photcalib']):
            stdout_write("\n\n\nExecuting photcalib...\n\n\n")
            exptime = ota_list[0].header['EXPTIME']
            titlestring = "%s\n(obsid: %s - filter: %s- exptime: %ds)" % (
                ota_list[0].header['OBJECT'],
                ota_list[0].header['OBSID'],
                ota_list[0].header['FILTER'],
                int(ota_list[0].header['EXPTIME']),
                )
            filter_name = get_valid_filter_name(ota_list[0].header)
            zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = \
                podi_photcalib.photcalib(global_source_cat, outputfile, filter_name, 
                                         exptime=exptime,
                                         diagplots=True,
                                         plottitle=titlestring,
                                         otalist=ota_list,
                                         options=options)

        ota_list[0].header.update("PHOTZP", zeropoint_median, "phot. zeropoint corr for exptime")
        ota_list[0].header.update("PHOTZPE", zeropoint_std, "zeropoint std.dev.")
        ota_list[0].header.update("PHOTZP_X", zeropoint_exptime, "phot zeropoint for this frame")


    #
    # If requested by user via command line:
    # Execute the fringing correction
    #
    if (not options['fringe_dir'] == None and not fringe_scaling == None): #.shape[0] > 0):
        
        # Determine the global scaling factor
        good_scalings = three_sigma_clip(fringe_scaling[:,6], [0, 1e9])
        fringe_scaling_median = numpy.median(good_scalings)
        fringe_scaling_std    = numpy.std(good_scalings)

        # and log the values in the primary header
        ota_list[0].header.update("FRNG_SCL", fringe_scaling_median)
        ota_list[0].header.update("FRNG_STD", fringe_scaling_std)

        # Construct the name of the fringe map
        filter_name = ota_list[0].header['FILTER']
        fringe_filename = check_filename_directory(options['fringe_dir'], "fringe__%s.fits" % (filter_name))
        # print fringe_filename
        if (os.path.isfile(fringe_filename)):
            fringe_hdulist = pyfits.open(fringe_filename)

            # Now do the correction
            for ext in range(1, len(ota_list)):
                extname = ota_list[ext].header['EXTNAME']
                if (type(ota_list[ext]) != pyfits.hdu.image.ImageHDU):
                    continue
                # print "subtracting",extname
                ota_list[ext].data -= (fringe_hdulist[extname].data * fringe_scaling_median)

    #print "Waiting for a bit"
    #afw.wait()
    #print "done waiting, writing output file"
    #print ota_list
    hdulist = pyfits.HDUList(ota_list)
    for i in range(1, len(hdulist)):
        if 'SIMPLE' in hdulist[i].header:
            del hdulist[i].header['SIMPLE']
    hdulist.verify()

    #print "hdulist=",hdulist

    if (not batchmode):
        stdout_write("writing output file (%s)..." % (outputfile))
        clobberfile(outputfile)
        hdulist.writeto(outputfile, clobber=True)
        # afw.write(hdulist, outputfile)
        stdout_write(" done!\n")
    else:
        stdout_write(" continuing ...")
        return hdulist

    # afw.finish(userinfo=True)

    return 0



def odi2wcs(xy, wcs):
    wcs.updateFromHeader()
    radec = numpy.ones((xy.shape[0],2)) * 1e9
    for i in range(xy.shape[0]):
        x, y = xy[i,0], xy[i,1]
        radec[i,0], radec[i,1] = wcs.pix2wcs(x, y)
    return radec

def wcs2odi(radec, wcs):
    wcs.updateFromHeader()
    xy = numpy.ones((radec.shape[0],2))
    for i in range(radec.shape[0]):
        xy[i,0], xy[i,1] = wcs.wcs2pix(radec[i,0], radec[i,1])
    return xy



def sexcat_to_tablehdu(fixwcs_ref_ra, fixwcs_ref_dec):
    
    columns = [\
        pyfits.Column(name='ra',      format='E', array=fixwcs_ref_ra),
        pyfits.Column(name='dec',     format='E', array=fixwcs_ref_dec),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("CAT.2MASS", comment=None)
    return tbhdu


def odi_sources_to_tablehdu(ota_list, fixwcs_odi_sourcecat):

    # # First of all update all Ra/Dec values
    # final_cat = None
    # for ext in range(1, len(ota_list)):
    #     extname = ota_list[ext].header['EXTNAME']
    #     if (extname[0:3] == "OTA" and extname[5:] == ".SCI"):
    #         ota = int(ota_list[ext].header['EXTNAME'][3:5])
    #         wcs = astWCS.WCS(ota_list[ext].header, mode="pyfits")

    #         in_this_ota = fixwcs_odi_sourcecat[:,12] == ota
    #         sources_ota = fixwcs_odi_sourcecat[in_this_ota]

    #         xy = sources_ota[:,2:4]
            
    #         sources_ota[:,0:2] = odi2wcs(xy, wcs)
            
    #         if (final_cat == None):
    #             final_cat = sources_ota
    #         else:
    #             final_cat = numpy.append(final_cat, sources_ota, axis=0)

    # source_info = [ ra, dec, center_x, center_y, fwhm_x, fwhm_y, amplitude, peak, 
    #            bg_level, bg_variance, s_n, area, extension_id]

    source_cat = fixwcs_odi_sourcecat

    columns = [\
        pyfits.Column(name='RA',             format='E', unit='degrees', array=source_cat[:, 0], disp='right ascension'),
        pyfits.Column(name='DEC',            format='E', unit='degrees', array=source_cat[:, 1], disp='declination'),
        pyfits.Column(name='X',              format='E', unit='pixel',   array=source_cat[:, 2], disp='center x'),
        pyfits.Column(name='Y',              format='E', unit='pixel',   array=source_cat[:, 3], disp='center y'),
        pyfits.Column(name='FWHM_IMAGE',     format='E', unit='pixel',   array=source_cat[:, 4], disp='FWHM in pixels'),
        pyfits.Column(name='FWHM_WORLD',     format='E', unit='deg',     array=source_cat[:, 5], disp='FWHM in degrees'),
        pyfits.Column(name='BACKGROUND',     format='E', unit='counts',  array=source_cat[:, 6], disp='background level'),
        pyfits.Column(name='FLAGS',          format='E', unit='',        array=source_cat[:, 7], disp='SExtractor flags'),
        pyfits.Column(name='OTA',            format='E', unit='',        array=source_cat[:, 8], disp='source OTA'),
        pyfits.Column(name='MAG_D05',        format='E', unit='mag',     array=source_cat[:, 9], disp='0.5 arcsec, 5 pixels'),
        pyfits.Column(name='MAGERR_D05',     format='E', unit='mag',     array=source_cat[:,17], disp=''),
        pyfits.Column(name='MAG_D08',        format='E', unit='mag',     array=source_cat[:,10], disp='0.8 arcsec, 7 pixels'),
        pyfits.Column(name='MAGERR_D08',     format='E', unit='mag',     array=source_cat[:,18], disp=''),
        pyfits.Column(name='MAG_D10',        format='E', unit='mag',     array=source_cat[:,11], disp='1.0 arcsec, 9 pixels'),
        pyfits.Column(name='MAGERR_D10',     format='E', unit='mag',     array=source_cat[:,19], disp=''),
        pyfits.Column(name='MAG_D15',        format='E', unit='mag',     array=source_cat[:,12], disp='1.5 arcsec, 14 pixels'),
        pyfits.Column(name='MAGERR_D15',     format='E', unit='mag',     array=source_cat[:,20], disp=''),
        pyfits.Column(name='MAG_D20',        format='E', unit='mag',     array=source_cat[:,13], disp='2.0 arcsec, 18 pixels'),
        pyfits.Column(name='MAGERR_D20',     format='E', unit='mag',     array=source_cat[:,21], disp=''),
        pyfits.Column(name='MAG_D25',        format='E', unit='mag',     array=source_cat[:,14], disp='2.5 arcsec, 23 pixels'),
        pyfits.Column(name='MAGERR_D25',     format='E', unit='mag',     array=source_cat[:,22], disp=''),
        pyfits.Column(name='MAG_D30',        format='E', unit='mag',     array=source_cat[:,15], disp='3.0 arcsec, 27 pixels'),
        pyfits.Column(name='MAGERR_D30',     format='E', unit='mag',     array=source_cat[:,23], disp=''),
        pyfits.Column(name='MAG_D35',        format='E', unit='mag',     array=source_cat[:,16], disp='3.5 arcsec, 32 pixels'),
        pyfits.Column(name='MAGERR_D35',     format='E', unit='mag',     array=source_cat[:,24], disp=''),
        pyfits.Column(name='FLUX_MAX',       format='E', unit='counts',  array=source_cat[:,25], disp='max count rate'),
        pyfits.Column(name='AWIN_IMAGE',     format='E', unit='pixel',   array=source_cat[:,26], disp='major semi-axis'),
        pyfits.Column(name='BWIN_IMAGE',     format='E', unit='pixel',   array=source_cat[:,27], disp='minor semi-axis'),
        pyfits.Column(name='THETAWIN_IMAGE', format='E', unit='degrees', array=source_cat[:,28], disp='position angle'),
        pyfits.Column(name='ELONGATION',     format='E', unit='',        array=source_cat[:,29], disp='elongation'),
        pyfits.Column(name='ELlPTICITY',     format='E', unit='',        array=source_cat[:,30], disp='ellipticity'),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("CAT.ODI", comment=None)
    return tbhdu, source_cat


def set_default_options(options_in=None):

    options = {}
    if (options_in != None):
        options = options_in

    # Get directory of the executable. This also serves as the 
    # fall-back directory for some data
    full_path = os.path.abspath(sys.argv[0])
    options['exec_dir'], dummy = os.path.split(full_path)

    options['update_persistency_only'] = False
    options['persistency_dir'] = None
    options["persistency_map"] = None
    options['max_persistency_time'] = 600

    options['fringe_dir'] = None
    options['fringe_vectors'] = "%s/.fringevectors/" % (options['exec_dir'])

    options['pupilghost_dir'] = None

    options['bias_dir'] = None
    options['dark_dir'] = None
    options['flat_dir'] = None
    options['bpm_dir']  = None

    options['nonlinearity'] = None

    options['fixwcs'] = False
    options['wcs_distortion'] = None

    options['indef_pixelvalue'] = numpy.NaN

    options['offset_pointing'] = [0,0]
    options['offset_dither'] = [0,0]
    options['target_coords'] = None

    options['verbose'] = False

    options['central_only'] = False

    options['bgmode'] = False

    options['photcalib'] = False

    options['plotformat'] = ['png']

    return options


def check_filename_directory(given, default_filename):
    if (given == None):
        return None
    elif (os.path.isfile(given)):
        return given
    elif (os.path.isdir(given)):
        return "%s/%s" % (given, default_filename)

    return ""

def read_options_from_commandline(options=None):

    if (options == None):
        options = set_default_options()

    options['verbose'] = cmdline_arg_isset("-verbose")

    # Handle all reduction flags from command line
    if (cmdline_arg_isset("-cals")):
        cals_dir = get_cmdline_arg("-cals")
        options['bias_dir'] = cals_dir
        options['dark_dir'] = cals_dir
        options['flat_dir'] = cals_dir

    options['bias_dir'] = cmdline_arg_set_or_default("-bias", options['bias_dir'])
    options['dark_dir'] = cmdline_arg_set_or_default("-dark", options['dark_dir'])
    options['flat_dir'] = cmdline_arg_set_or_default("-flat", options['flat_dir'])

    options['bpm_dir']  = cmdline_arg_set_or_default("-bpm", options['bpm_dir'])
    if (options['bpm_dir'] == "auto"):
        options['bpm_dir'] = options['exec_dir']
        
    if (options['verbose']):
        print """
Calibration data:
            Bias: %s
            Dark: %s
      Flatfields: %s
  Bad pixel mask: %s
""" % (options['bias_dir'], options['dark_dir'], options['flat_dir'], options['bpm_dir'])

    options['persistency_dir'] = cmdline_arg_set_or_default('-persistency', None)
    if (not options['persistency_dir'] == None):
        # This automatically creates the index.cat file
        mjd_catalog_list = podi_persistency.get_list_of_saturation_tables(options['persistency_dir'])

    options["update_persistency_only"] = cmdline_arg_isset("-update_persistency_only")

    options['fringe_dir'] = cmdline_arg_set_or_default('-fringe', None)
    options['fringe_vectors'] = cmdline_arg_set_or_default("-fringevectors", options['fringe_vectors'])

    options['pupilghost_dir'] = cmdline_arg_set_or_default('-pupilghost', None)

    options['fixwcs'] = cmdline_arg_isset("-fixwcs")
    # For now assume that the WCS template file is located in the same directory as the executable
    root_dir, py_exe = os.path.split(os.path.abspath(sys.argv[0]))
    options['wcs_distortion'] = root_dir + "/wcs_distort2.fits"
    options['wcs_distortion'] = cmdline_arg_set_or_default("-wcs", options['wcs_distortion'])
    if (not os.path.isfile(options['wcs_distortion'])):
        options['wcs_distortion'] = None

    options['clobber'] = not cmdline_arg_isset("-noclobber")
    
    # Set the fallback value for undefined pixels (mostly the gaps between the OTA cells)
    # This works perfectly fine in ds9, but not SExtractor
    if (cmdline_arg_isset("-prep4sex")):
        # Check if the user requested us to prepare the frame for SExtractor
        # SExtractor doesn't like NaNs, so replace all of them with something
        # more negative than -1e30 (that's -1 times SEx's BIG variable)
        options['indef_pixelvalue'] = -1e31
    
    try:
        tmp = float(cmdline_arg_set_or_default('-indefval', numpy.NaN))
        options['indef_pixelvalue'] = tmp
    except:
        pass

    if (cmdline_arg_isset("-wcsoffset")):
        tmp = get_cmdline_arg("-wcsoffset")
        items = tmp.split(',')
        options['offset_pointing'] = [float(items[0]), float(items[1])]
        stdout_write("Applying a user-defined WCS offset of %.3f, %.3f degrees\n" % (options['offset_pointing'][0], options['offset_pointing'][1]))

    #
    # Read all offsets from command line
    # For convenience, there are two sets of offset parameters, that internally simply 
    # get added up. The reason for this is to make specifying them on the command line 
    # easier, since the pointing offsets stay constant across a dither pattern, while 
    # the dither offsets change.
    #
    options['target_coords'] = None
    options['offset_pointing'] = None
    options['offset_dither'] = None
    # -target: overwrites the pointing information from the wcs header
    if (cmdline_arg_isset("-target")):
        _target_coords = cmdline_arg_set_or_default("-target", "0,0")
        ra,dummy,dec = _target_coords.partition(",")
        options['target_coords'] = (ra, dec)
    # -pointing: applies a given offset to the pointing position
    if (cmdline_arg_isset("-pointing")):
        _offset_pointing = cmdline_arg_set_or_default("-pointing", "0,0")
        dx,dummy,dy = _offset_pointing.partition(",")
        options['offset_pointing'] = [float(dx), float(dy)]
    # -dither: identical to -pointing
    if (cmdline_arg_isset("-dither")):
        _offset_dither = cmdline_arg_set_or_default("-dither", "0,0")
        dx,dummy,dy = _offset_dither.partition(",")
        options['offset_dither'] = [float(dx), float(dy)]
    #  .
    # /-\
    #  |   This section is likely outdated 
    #

    options['central_only'] = cmdline_arg_isset("-centralonly")

    options['bgmode'] = cmdline_arg_isset("-bgmode")

    options['photcalib'] = cmdline_arg_isset("-photcalib")

    options['nonlinearity'] = cmdline_arg_set_or_default("-nonlinearity", None)

    if (cmdline_arg_isset('-plotformat')):
        inputstr = cmdline_arg_set_or_default("-plotformat", "png")
        options['plotformat'] = inputstr.split(",")
        print "writing plots as ",options['plotformat']
        


    return options

if __name__ == "__main__":

    # Read the input directory that contains the individual OTA files
    input = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    if (len(get_clean_cmdline())>2):
        outputfile = get_clean_cmdline()[2]
    else:
        print "No output filename has been given, setting to default mergedcells.fits"
        outputfile = "mergedcells.fits"
    print "Writing results into",outputfile

    # Set the options for collectcells to some reasonable start values
    options = set_default_options()

    # Then read the actual given parameters from the command line
    options = read_options_from_commandline(options)

    # Collect all cells, perform reduction and write result file
    try:
        if (cmdline_arg_isset('-profile')):
            options['profile'] = True
            import cProfile, pstats
            cProfile.run("""collectcells(input, outputfile,
                     options=options)""", "profiler")
            p = pstats.Stats("profiler")
            p.strip_dirs().sort_stats('time').print_stats()
            p.sort_stats('time').print_stats()
        else:
            collectcells(input, outputfile, options=options)
    except:
        stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
        etype, error, stackpos = sys.exc_info()
        stdout_write("# Exception report:")
        stdout_write("#  ==> %s\n" % (error))
        print traceback.format_exc()
        stdout_write("#\n##############################\n")

