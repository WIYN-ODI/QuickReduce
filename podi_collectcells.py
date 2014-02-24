#!/usr/bin/env python
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
Function
--------------------------------------------------------------------------------
**podi_collectcells** is the main script of the pODI QuickReduce pipeline. As
such, it handles all reduction steps, combining all cells from all available 
OTAs into a single FITS-file. If invoked with the respective command line flags,
it also performs the astrometric and photometric calibration.

How to run podi_collectcells
--------------------------------------------------------------------------------
``podi_collectcells input output.fits (flags)``
Input here can be either
- a directory containing all OTA fits files. In this case, the directory has to
  contain the basename of the files, for example ``o20130101T123456.0/``.
- the filename of a single OTA FITS file. It is then assumed that all other OTA
  files reside in the same directory

Available flags
--------------------------------------------------------------------------------
* **-cals=/dir/to/cals**

  Specifies the directory that holds the calibration data. -cals is processed 
  first, but you can change the path to individual calibration frames with any 
  of the next options.


* **-bias=/dir/to/bias**

  Overwrites the -cals path to specify a different directory for the bias frames


* **-dark=/dir/to/darks**

  Same for darks

* **-flat=/dir/to/flats**

  Same for flats -

* **-bpm=/dir/to/bpms**

  Same as above. Alternatively, you can specify "auto" (without the " ) as 
  directory. In that case the script looks for the bad pixel masks in the 
  directory of the script.

* **-pupilghost=/dir**

  same as above for the pupil ghost file

* **-fringe=/dir**

  same as above for the fringing templates

* **-fringevectors**

  Fringe correction needs special "vector" frames marking dark and bright 
  regions. These are stored in a set of ds9 region files, with one file for each
  detector, and this option allows to specify the directory holding all these 
  files.

* **-gain**

  Apply a gain correction to all cells. If combined with the -nonlinearity 
  option, the relative gain values from the non-linearity are used, otherwise it
  uses the constants defined in each cell header. 

  Note that for correct reduction, this flag needs to be applied consistently,
  in particular when creating the bias, dark, and flat-field frames.

* **-persistency=/dir**

  same as above. However, since persistency files are both read (from the 
  earlier exposures) AND written (to keep track of saturated stars in the 
  current exposure) make sure you have permission to write into the specified 
  directory as well as sufficient disk-space.

* **-wcs=some/file.fits**

  Give the filename of the canned WCS solution that specifies the detector 
  distortion. See scamp2minifits on a way to create the necessary input file.

* **-fixwcs**

  Activates the automatic WCS correction routine. If enabled, collectcells run 
  a source detection (podi_findstars), obtains a catalog of reference star 
  positions (podi_ipprefcat) and corrects telescope pointing errors (using 
  functions from podi_fixwcs).

* **-simplewcs**

  Modifies the WCS solution that remains in the fits-frame. With this option, 
  the RA-ZPX solution is changed into a simple RA-TAN system without any 
  distortion factors. Mostly for debugging and testing.

* **-pointing**

* **-dither**

* **-target**

  Now defunct

* **-prep4sex**

  Changes the value of undefined pixels (those in the cell-gaps and the detector
  edges) from the default NaN (Not a Number) value to -1e31, a value exceeding 
  the SExtractor value for bad pixels, telling SExtractor to ignore them

* **-singleota=XX**

  Limits the output to only a specific OTA. For example, to only extract OTA42, 
  give -singleota=42.

* **-noclobber**

  Skip the data reduction if the output file already exists

* **-wcsoffset=d_ra,d_dec**

  Apply WCS offset in degrees to WCS solution before files are fed to 
  SourceExtractor. Use this option if the telescope pointing is WAY off.

* **-centralonly**

  Only collect cells from the central 3x3 array, and ignore all data from the 
  outlying OTAs

* **-bgmode**

  Disable the initial welcome splash screen

* **-photcalib**

  Enable the photometric calibration of the output frames.

* **-nonlinearity=/some/dir/nonlinearity_coefficients.fits**

  Activate the non-linearity correction and use the specified file as source 
  for all correction coefficients.

* **-plotformat=format1,format2**

  Generate the diagnostic plots in the specified formats. All formats supported
  by matplotlib are available - these are typically png, jpg, ps, pdf, svg

* **-nootalevelplots**

  By default, all diagnostic plots are created for the full focalplane and for 
  each OTA individually. Adding this option disables the OTA level, restricting
  diagnostic plots to only the focalplane views.

* **-qasubdirs**

  Creates a directory structure parallel to the output file and sorts the 
  diagnostic plots into (sub-)directories. The default directory structure is:

  output.fits
  QA/
  wcs1.png
  wcs1/
  OTA00.png
  OTA22.png
  seeing.png
  seeing/
  OTA00.png
  OTA22.png

* **-qasubdirname**

  Allows to specify the name of the "QA" directory as shown above.

* **-noqaplots**

  Disable the creation of the QA plots.

* **-addfitskey=KEY,value**

  Add a user-defined fits header entry to the primary header of the output file.
  KEY needs to be fits-conform, i.e. 8 characters or less without spaces or 
  special characters. value is stored as a string in the output file. If value 
  contains spaces, the value needs to be quoted (eg. -addfitskey=WIYN,"greatest 
  of the universe")

* **-nonsidereal=a,b,c**

  a: proper motion of the object (dRA*cos(dec)) given in arcseconds per hour.
  b: proper motion dDec in arcsec / hour
  c: Reference. This can either be a MJD directly, or a FITS file containing the 
     MJD-OBS header keyword which is then taken to be the reference MJD.

* **-fitradialZP**

  Fit a linear radial profile to the photometric zeropoint data. The center is 
  assumed to be roughly the center of OTA 3,3 as per Daniel Harbeck's demo-
  script. Fit results (RADZP_P0, RADZP_P1) and errors (RADZP_E0, RADZP_E1) are 
  stored in the primary header of the output file. If enabled, this options also
  generates another diagnostic plot showing the uncorrected photometric ZP as
  function of radius before and after the best-fit has been taken out.

Methods
--------------------------------------------------------------------------------

"""

import sys
import os
import signal
import pyfits
import numpy
import scipy
import scipy.optimize
import ephem
import traceback
#import psutil
import datetime

import Queue
import threading
import multiprocessing
import ctypes
import time
import logging
import itertools

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
import podi_logging
import podi_cosmicrays

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
    """
    collect_reduce_ota does work relating to the initial instrument 
    detrending **for a single OTA** from cross-talk correction, overscan-
    subtraction, bias & dark correction and flat-fielding. It also derives a 
    sky-background level and creates a source catalog using SourceExtractor.

    Parameters
    ----------
    filename : string

        Filename of the raw OTA FITS file to be reduced

    verbose : Bool

        Output more detailed progress information during execution. Mostly made
        for debugging and error tracking.

    options : struct/dictionary

        This dictionary contains all configuration parameters necessary for the
        reduction, such as flat-field directories and specified reduction flags
        read from the command-line.

    Returns
    -------
    data_products : dictionary

        The returned dictionary contains all information generated during the 
        reduction, such as the reduced ImageHDU, the source catalog, 
        data about background levels, etc.

    """

    data_products = {
        "hdu": None,
        "wcsdata": None,
        "sky-samples": None,
        "sky": None,
        "tech-header": None,
        }
    
    if (not os.path.isfile(filename)):
        stdout_write("Couldn't find file %s ..." % (filename))
    else:
        # Create an fits extension to hold the output
        hdu = pyfits.ImageHDU()
        log_svn_version(hdu.header)

        # also create a tech-header to keep track of values used in this frame
        tech_header = pyfits.Header()

        # Keep track of what input files we used
        reduction_files_used = {}

        try:
            hdulist = pyfits.open(filename, memmap=False)
        except:
            # This happed with corrupt files
            return data_products
            
        reduction_files_used['raw'] = filename

        # Get the binning factor
        binning = get_binning(hdulist[1].header)
        
        # Get the image output size (depends on binning)
        size_x, size_y = get_collected_image_dimensions(binning)
        #print size_x, size_y

        obsid = hdulist[0].header["OBSID"]
        ota = int(hdulist[0].header['FPPOS'][2:])
        ota_c_x, ota_c_y = int(math.floor(ota/10)), int(math.fmod(ota,10))

        # Save the fppos as name for this extension
        ota_name = "OTA%02d" % ota
        extname = "OTA%02d.SCI" % ota
        hdu.update_ext_name(extname)
        hdu.header['OTA'] = (ota, "OTA designation")

        # podi_logging.podi_getlogger("%s - %s" % (obsid, extname), options['log_setup'])
        logger = logging.getLogger("%s - %s" % (obsid, extname))

        # Now copy the headers from the original file into the new one
        cards = hdulist[0].header.cards
        for (keyword, value, comment) in cards:
            hdu.header[keyword] = (value, comment)

        # Also read the MJD for this frame. This will be needed later for the correction
        mjd = hdu.header['MJD-OBS']

        #
        # in the old data, some header EXTNAMES are duplicated, so we can't look
        # up extensions by name in this case. Create a dictionary with EXTNAMES 
        # and extension IDs to work around this problem.
        #
        extname2id = {}
        for i in range(1, len(hdulist)):
            _extname = hdulist[i].header['EXTNAME']
            extname2id[_extname.lower()] = i
        logger.debug("Setting up extension lookup directory:\n"+str(extname2id))

        # 
        # Perform cross-talk correction, using coefficients found in the 
        # podi_crosstalk package.
        #
        # Procedure: 
        # Correct all frames for the crosstalk, then write them back into the 
        # original OTA position so we can perform the overscan subtraction etc.
        #

        logger.info("Starting work on OTA %02d of %s ..." % (ota, obsid))

        # Allocate some memory for the cells in one row
        xtalk_corr = [None] * 8

        # Now go through each of the 8 lines
        logger.debug("Starting crosstalk correction (%s)" % (extname))
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

                    correction = hdulist[extname2id[xy_name]].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                    if (column != xtalk_column):
                        saturated = hdulist[extname2id[xy_name]].data >= podi_crosstalk.xtalk_saturation_limit
                        correction[saturated] = -1 * podi_crosstalk.xtalk_saturated_correction

                    xtalk_corr[column] += correction #hdulist[xy_name].data * podi_crosstalk.xtalk_matrix[extname][row][xtalk_column][column]
                    #print xtalk_corr[column][100,100]

            for column in range(8):
                # Now all cells in this row have been corrected, let's write them 
                # back into the hdulist so can can continue with the overscan subtraction etc.
                xy_name = "xy%d%d" % (column, row)
                hdulist[extname2id[xy_name]].data = xtalk_corr[column]
        logger.debug("Done with crosstalk correction")

        #
        # Allocate memory for the merged frame, and set all pixels by default to NaN.
        # Valid pixels will subsequently be overwritten with real numbers
        #
        merged = numpy.ones(shape=(size_x, size_y), dtype=numpy.float32)
        merged[:,:] = options['indef_pixelvalue']
        
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

        nonlin_data = None
        if (options['nonlinearity-set'] or options['gain_method'] == "relative"):
            nonlinearity_file = options['nonlinearity']
            if (options['nonlinearity'] == None or 
                options['nonlinearity'] == "" or
                not os.path.isfile(nonlinearity_file)):
                nonlinearity_file = podi_nonlinearity.find_nonlinearity_coefficient_file(mjd, options)
            if (options['verbose']):
                print "Using non-linearity coefficients from",nonlinearity_file
            logger.debug("Using non-linearity coefficients from file %s"  % (nonlinearity_file))
            nonlin_data = podi_nonlinearity.load_nonlinearity_correction_table(nonlinearity_file, ota)
            if (options['nonlinearity-set']):
                reduction_files_used['nonlinearity'] = nonlinearity_file

        #
        # Search, first in the flat-field, then in the bias-frame for a 
        # TECHDATA extension that holds GAINS and readout noises for each cell. 
        #
        techdata = None
        if (not options['techdata'] == None):
            techfile = None

            # Set the filename for the TECHDATA extension based on the user-option
            if (options['techdata'] == "from_flat" and not options['flat_dir'] == None):
                techfile = check_filename_directory(options['flat_dir'], "flat_%s_bin%d.fits" % (filter_name, binning))
                
            elif (options['techdata'] == "from_bias" and not options['bias_dir'] == None):
                techfile = check_filename_directory(options['bias_dir'], "bias_bin%s.fits" % (binning))

            else:
                techfile= check_filename_directory(options['techdata'], "techdata_%s_bin%d.fits" % (filter_name, binning))
                
            # Check if the specified file exists and read the data if possible
            if (os.path.isfile(techfile)):
                logger.debug("Reading techdata from file %s" % (techfile))
                techhdulist = pyfits.open(techfile)
                try:
                    techdata = techhdulist['TECHDATA'].header
                    reduction_files_used['techdata'] = techfile
                except:
                    pass
                techhdulist.close()
            else:
                basedir, _ = os.path.split(os.path.abspath(sys.argv[0]))
                techfile = "%s/techdata.fits" % (basedir)
                if (os.path.isfile(techfile)):
                    logger.debug("Reading techdata from file %s" % (techfile))
                    techhdulist = pyfits.open(techfile)
                    try:
                        techdata = techhdulist['TECHDATA'].header
                        reduction_files_used['techdata'] = techfile
                    except:
                        pass
                    techhdulist.close()
                else:
                    logger.debug("Was looking for techfile %s but couldn't find it" % (techfile))

           
        all_gains = numpy.ones(shape=(8,8)) * -99
        all_readnoise = numpy.ones(shape=(8,8)) * -99
        all_readnoise_electrons = numpy.ones(shape=(8,8)) * -99

        for wm_cellx, wm_celly in itertools.product(range(8),repeat=2):
            #if (not options['bgmode']):
            #    stdout_write("\r%s:   OTA %02d, cell %s ..." % (obsid, ota, hdulist[cell].header['EXTNAME']))
            cellname = "xy%d%d" % (wm_cellx, wm_celly)
            cell = extname2id[cellname]

            #
            # Special case for cell 0,7 (the one in the bottom left corner):
            # Copy the CRPIX values into the merged image header 
            #
            if (hdulist[cell].header['EXTNAME'].lower() == "xy07"):
                # print "Setting CRPIXs", hdulist[cell].header['CRPIX1'], hdulist[cell].header['CRPIX2']
                hdu.header["CRPIX1"] = (hdulist[cell].header['CRPIX1'], "Ref. pixel RA")
                hdu.header["CRPIX2"] = (hdulist[cell].header['CRPIX2'], "Ref. pixel DEC")



            # Check if this is one of the broken cells
            cellmode_id = get_cellmode(hdulist[0].header, hdulist[cell].header)
            if (not cellmode_id == 0):
                # This means it either broken (id=-1) or in video-mode (id=1)
                continue

            # logger.debug("ota %02d, cell %d,%d: gain=%f, ron=%f, ron(e-)=%f" % (
            #     ota, wm_cellx, wm_celly, gain, readnoise, readnoise_electrons))
            # all_gains[wm_cellx, wm_celly] = gain
            # all_readnoise[wm_cellx, wm_celly] = readnoise
            # all_readnoise_electrons[wm_cellx, wm_celly] = readnoise_electrons

            #
            # Now extract just the data section.
            # Values are hard-coded as some of the DATASEC keywords are wrong
            #
            datasec = extract_datasec_from_cell(hdulist[cell].data, binning) #[0:494, 0:480] 
            # print "datasec.shape=",datasec.shape
            # Now overscan subtract and insert into large frame
            overscan_region = extract_biassec_from_cell(hdulist[cell].data, binning)
            #extract_region(hdulist[cell].data, '[500:530,1:494]')
            overscan_level = numpy.median(overscan_region)

            datasec -= overscan_level

            # logger.debug("cell %s: gain=%.2f overscan=%6.1f" % (hdulist[cell].header['EXTNAME'], gain, overscan_level))

            if (options['nonlinearity-set']):
                nonlin_correction = podi_nonlinearity.compute_cell_nonlinearity_correction(
                    datasec, wm_cellx, wm_celly, nonlin_data)
                datasec += nonlin_correction

            #
            # Insert the reduced data-section of this cell into the large OTA frame
            #
            bx, tx, by, ty = cell2ota__get_target_region(wm_cellx, wm_celly, binning)
            merged[by:ty,bx:tx] = datasec

            #
            # Now that the cell is pre-reduced (merged and overscan-subtracted), assemble 
            # the tech-header by copying the information from the input techdata 
            # to the output techdata
            #
            if (not techdata == None):
                ids = "%02d%d%d" % (ota, wm_cellx, wm_celly)
                for keybase in techdata_keywords:
                    keyword = keybase+ids
                    if (keyword in techdata):
                        key, val, com = techdata.cards[keyword]
                        tech_header[key] = (val, com)
                # Also set the gain and readnoise values to be used later for the gain correction
                all_gains[wm_cellx, wm_celly] = \
                    techdata['GN__'+ids] if ('GN__'+ids in techdata) else \
                        hdulist[cell].header['GAIN'] if 'GAIN' in hdulist[cell].header else backup_gain
                all_readnoise[wm_cellx, wm_celly] = \
                    techdata['RN__'+ids] if 'RN__'+ids in techdata else backup_readnoise
                all_readnoise_electrons[wm_cellx, wm_celly] = \
                    techdata['RNE_'+ids] if 'RNE_'+ids in techdata else backup_readnoise_electrons
            else:
                all_gains[wm_cellx, wm_celly] = \
                    hdulist[cell].header['GAIN'] if 'GAIN' in hdulist[cell].header else backup_gain
                all_readnoise[wm_cellx, wm_celly] = backup_readnoise
                all_readnoise_electrons[wm_cellx, wm_celly] = backup_readnoise_electrons
                
            # work on next cell

        logger.debug("Collected all cells for OTA %02d of %s" % (ota, obsid))
        # for c in tech_header: print tech_header.cards[c]

        # 
        # At this point we have a 4x4 Kpixel array with all cells merged
        #

        # If we are to do some bias subtraction:
        if (not options['bias_dir'] == None):
            bias_filename = check_filename_directory(options['bias_dir'], "bias_bin%s.fits" % (binning))
            if (os.path.isfile(bias_filename)):
                bias = pyfits.open(bias_filename)
                reduction_files_used['bias'] = bias_filename

                # Search for the bias data for the current OTA
                for bias_ext in bias[1:]:
                    if (not is_image_extension(bias_ext)):
                        continue
                    fppos_bias = bias_ext.header['FPPOS']
                    if (fppos_bias == fppos):
                        # This is the one
                        merged -= bias_ext.data
                        logger.debug("Subtracting bias: %s" % (bias_filename))
                        break

                bias.close()
                hdu.header.add_history("CC-BIAS: %s" % (os.path.abspath(bias_filename)))
                del bias
 
        #
        # Do some dark subtraction:
        # Add at some point: use different darks for all detectors switched on 
        # to minimize the integration glow in guide-OTAs
        #
        if (not options['dark_dir'] == None):

            # For now assume all detectors are switched on
            detectorglow = "yes"

            dark_filename = check_filename_directory(options['dark_dir'], "dark_%s_bin%d.fits" % (detectorglow, binning))
            if (os.path.isfile(dark_filename)):
                dark = pyfits.open(dark_filename)
                reduction_files_used['dark'] = dark_filename
                darktime = dark[0].header['EXPMEAS']

                # Search for the flatfield data for the current OTA
                for dark_ext in dark[1:]:
                    if (not is_image_extension(dark_ext)):
                        continue
                    fppos_dark = dark_ext.header['FPPOS']
                    if (fppos_dark == fppos):
                        # This is the one
                        dark_scaling = exposure_time / darktime
                        logger.debug("Subtracting dark: %s (scaling=%.2f)" % (dark_filename, dark_scaling))
                        merged -= (dark_ext.data * dark_scaling)
                        break

                dark.close()
                hdu.header.add_history("CC-DARK: %s" % (os.path.abspath(dark_filename)))
                del dark
 
        # By default, mark the frame as not affected by the pupilghost. This
        # might be overwritten if the flat-field has PG specifications.
        hdu.header['PGAFCTD'] = False

        # If the third parameter points to a directory with flat-fields
        if (not options['flat_dir'] == None):
            flatfield_filename = check_filename_directory(options['flat_dir'], "flat_%s_bin%d.fits" % (filter_name, binning))
            if (os.path.isfile(flatfield_filename)):
                flatfield = pyfits.open(flatfield_filename)
                reduction_files_used['flat'] = flatfield_filename

                # Search for the flatfield data for the current OTA
                for ff_ext in flatfield[1:]:
                    if (not is_image_extension(ff_ext)):
                        continue
                    fppos_flatfield = ff_ext.header['FPPOS']
                    if (fppos_flatfield == fppos):
                        # This is the one
                        merged /= ff_ext.data
                        logger.debug("Dividing by flatfield: %s" % (flatfield_filename))

                        # If normalizing with the flat-field, overwrite the gain
                        # keyword with the average gain value of the flatfield.
                        ff_gain = flatfield[0].header['GAIN'] \
                                  if 'GAIN' in flatfield[0].header else 1.3
                        hdu.header['GAIN'] = ff_gain
                        logger.debug("Overwriting GAIN keyword: %f" % (ff_gain))
                        
                        logger.debug("Checking if extension has PGAFCTD keyword: %s" % (str('PGAFCTD' in ff_ext.header)))
                        if ('PGAFCTD' in ff_ext.header):
                            logger.debug("Value of PGAFCTD header keyword: %s" % (str(ff_ext.header['PGAFCTD'])))
                        if ('PGAFCTD' in ff_ext.header and ff_ext.header['PGAFCTD']):
                            # Mark this extension as having a pupilghost problem.
                            hdu.header['PGAFCTD'] = True

                            # Also copy the center position of the pupilghost
                            # If this is not stored in the flat-field, assume some
                            # standard coordinates
                            logger.debug("Found an extension affected by pupilghost: %s" % (extname))
                            if (extname in pupilghost_centers):
                                pupilghost_center_x, pupilghost_center_y = pupilghost_centers[extname]
                                hdu.header['PGCNTR_X'] = (pupilghost_center_x, "pupil ghost center position X")
                                hdu.header['PGCNTR_Y'] = (pupilghost_center_y, "pupil ghost center position Y")
                            for pgheader in (
                                    'PGCENTER', 'PGCNTR_X', 'PGCNTR_Y',
                                    'PGTMPL_X', 'PGTMPL_Y', 
                                    'PGREG_X1', 'PGREG_X2', 'PGREG_Y1', 'PGREG_Y2',
                                    'PGEFCTVX', 'PGEFCTVY', 
                                    'PGSCALNG',
                                    'PGROTANG',
                            ):
                                if (pgheader in ff_ext.header):
                                    hdu.header[pgheader] = ff_ext.header[pgheader]

                            # if ('PGCNTR_X' in ff_ext.header):
                            #     pupilghost_center_x = ff_ext.header['PGCNTR_X']
                            #     hdu.header['PGCNTR_X'] = (pupilghost_center_x, "pupil ghost center position X")
                            # if ('PGCNTR_Y' in ff_ext.header):
                            #     pupilghost_center_y = ff_ext.header['PGCNTR_Y']
                            #     hdu.header['PGCNTR_Y'] = (pupilghost_center_y, "pupil ghost center position Y")
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
                logger.debug("Applying BPM file: %s" % (region_file))
                mask_broken_regions(merged, region_file)
                hdu.header.add_history("CC-BPM: %s" % (os.path.abspath(region_file)))
                reduction_files_used['bpm'] = region_file

        #
        # Now apply the gain correction
        #

        # Compute the gain across this OTA
        valid_gain = (all_gains > 0) & (numpy.isfinite(all_gains))
        gain_avg = numpy.mean(all_gains[valid_gain])
        hdu.header['GAIN_RAW'] = (gain_avg, 'gain averaged over all cells')
        hdu.header['NGAINRAW'] = (numpy.sum(valid_gain), 'number of cells contrib. to avg. gain')

        logger.debug("GAIN setting:"+str(options['gain_correct']))
        if (options['gain_correct']):
            # Correct for the gain variations in each cell
            logger.debug("Applying gain correction (OTA %02d) - method: %s" % (ota, 
                options['gain_method'] if not options['gain_method'] == None else "default:techdata"))

            if (options['gain_method'] == 'relative'):
                reduction_files_used['gain'] = nonlinearity_file
                # Find the relative gain correction factor based on the non-linearity correction data
                logger.debug("Apply gain correction from nonlinearity data")

                for cx, cy in itertools.product(range(8), repeat=2):
                    x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy, binning)
                    merged[y1:y2, x1:x2], gain = podi_nonlinearity.apply_gain_correction(
                        merged[y1:y2, x1:x2], cx, cy, nonlin_data, return_gain=True)
                    all_gains[cx,cy] /= gain
                
            elif (options['gain_method'] == 'header'):
                logger.debug("Applying gain correction  with GAINS from header")
                reduction_files_used['gain'] = filename
                for cx, cy in itertools.product(range(8), repeat=2):
                    x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy, binning)

                    cell = extname2id["xy%d%d" % (cx,cy)]
                    # Use what's in the header
                    gain = float(hdulist[cell].header['GAIN']) if 'GAIN' in hdulist[cell].header else backup_gain
                    merged[y1:y2, x1:x2] *= gain
                    all_gains[cx,cy] /= gain

            else:
                if (not techdata == None):
                    logger.debug("Using GAINs from tech-data file %s" % (techfile))
                    reduction_files_used['gain'] = techfile
                    for cx, cy in itertools.product(range(8), repeat=2):
                        x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy, binning)

                        gain = all_gains[cx,cy]
                        if (gain <= 0):
                            all_gains[cx,cy] = -99
                            continue
                        # print "Applying gain", gain, "to cell", cx, cy
                        merged[y1:y2, x1:x2] *= gain
                        all_gains[cx,cy] /= gain
                else:
                    logger.warning("GAIN correction using TECHDATA requested, but can't find a tech-data file")


        #
        # Compute the average gain, read-noise, etc and store how many cells contributed
        #
        # Set the GAIN to be used for photometry. This gain is 1.0 if we 
        # already applied the gain correction to the data
        valid_gain = (all_gains > 0) & (numpy.isfinite(all_gains))
        gain_avg = numpy.mean(all_gains[valid_gain])
        hdu.header['GAIN'] = (gain_avg, 'gain averaged over all cells')
        hdu.header['NGAIN'] = (numpy.sum(valid_gain), 'number of cells contrib. to avg. gain')

        valid_ron = (all_readnoise > 0) & (numpy.isfinite(all_readnoise))
        ron_avg = numpy.mean(all_readnoise[valid_ron])
        hdu.header['RDNOISE'] = (ron_avg, 'readnoise averaged over all cells')
        hdu.header['NRDNOISE'] = (numpy.sum(valid_ron), 'number of cells contrib. to avg. readnoise')

        valid_rone = (all_readnoise_electrons > 0) & (numpy.isfinite(all_readnoise_electrons))
        rone_avg = numpy.mean(all_readnoise_electrons[valid_rone])
        hdu.header['RDNOISEE'] = (rone_avg, 'readnoise in e- averaged over all cells')
        hdu.header['NRDNOSEE'] = (numpy.sum(valid_rone), 'number cells contrib. to avg. readnoise in e-')

        
        
                
        #
        # Optionally, normalize the frame by some header value
        #
        if ("normalize" in options):
            norm_factor = 1
            if (options['normalize'] == "EXPTIME" or
                options['normalize'] == "EXPMEAS"):
                # get the factor from the header
                norm_factor = hdu.header[options['normalize']]

                if (norm_factor > 0):
                    logger.debug("normalizing data with constant %f (%s)" % (
                            norm_factor, options['normalize']))
                    # normalize the data
                    merged /= norm_factor
                    # and fix the EXPTIME/EXPMEAS header as well
                    hdu.header['EXPTIME'] /= norm_factor
                    hdu.header['EXPMEAS'] /= norm_factor
            hdu.header['NORMALIZ'] = (norm_factor, "normalization constant")
                
        #
        # If persistency correction is requested, perform it now
        #
        if (options['persistency_dir'] != None):
            logger.debug("Applying persistency correction")
            if (options['verbose']): stdout_write("Applying persistency correction\n")

            # Get a list of all saturation catalogs
            full_filelist = podi_persistency.get_list_of_saturation_tables(options['persistency_dir'])
            # print full_filelist

            # Search for the saturation map for this frame and, if found,
            # mask out all saturated pixels
            #if (verbose): 
            logger.debug("MJD of this frame: %f" % (mjd))

            saturated_thisframe = podi_persistency.select_from_saturation_tables(full_filelist, mjd, None)
            reduction_files_used['saturation'] = saturated_thisframe
            #print "this frame=",saturated_thisframe,mjd

            if (not saturated_thisframe == None):
                merged = podi_persistency.mask_saturation_defects(saturated_thisframe, ota, merged)
            hdu.header["SATMASK"] = saturated_thisframe

            # Also pick all files within a given MJD range, and apply the 
            # appropriate correction (which, for now, is simply masking all pixels)
            filelist = podi_persistency.select_from_saturation_tables(full_filelist, mjd, [1,options['max_persistency_time']])
            if (len(filelist) > 0):
                # Extract only the filenames from the filelist dictionary
                persistency_files_used = []
                for assoc_persistency_file, assoc_mjd in filelist.iteritems():
                    persistency_files_used.append(assoc_persistency_file)
                reduction_files_used['persistency'] = persistency_files_used

                merged = podi_persistency.correct_persistency_effects(ota, merged, mjd, filelist)
                persistency_catalog_counter = 0
                for filename in filelist:
                    persistency_catalog_counter += 1
                    keyname = "PERSIS%02d" % (persistency_catalog_counter)
                    hdu.header[keyname] = filename


        #
        # If requested, determine the optimal fringe scaling
        #
        fringe_scaling = None
        fringe_scaling_median, fringe_scaling_std = -1, -1
        if (options['fringe_dir'] != None):
            fringe_filename = check_filename_directory(options['fringe_dir'], "fringe__%s.fits" % (filter_name))
            fringe_vector_file = "%s/fringevectors__%s__OTA%02d.reg" % (options['fringe_vectors'], filter_name, ota)
            reduction_files_used['fringemap'] = fringe_filename
            reduction_files_used['fringevector'] = fringe_vector_file
            logger.debug("fringe file %s found: %s" % (fringe_filename, os.path.isfile(fringe_filename)))
            logger.debug("fringe vector %s found: %s" % (fringe_vector_file, os.path.isfile(fringe_vector_file)))
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
            hdu.header["FRNG_SCL"] = fringe_scaling_median
            hdu.header["FRNG_STD"] = fringe_scaling_std
            hdu.header["FRNG_OK"] = (fringe_scaling != None)

        # Insert the DETSEC header so IRAF understands where to put the extensions
        start_x = ota_c_x * size_x #4096
        start_y = ota_c_y * size_y #4096        
        end_x = start_x + size_x #det_x2 - det_x1
        end_y = start_y + size_y #det_y2 - det_y1
        detsec_str = "[%d:%d,%d:%d]" % (start_x, end_x, start_y, end_y)
        hdu.header["DETSEC"] = (detsec_str, "position of OTA in focal plane")
                
        if (cmdline_arg_isset("-simplewcs") or options['wcs_distortion'] != None):
            logger.debug("Clearing existing WCS solution")
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
        hdu.header["CUNIT1"] = ("deg", "")
        hdu.header["CUNIT2"] = ("deg", "")

        logger.debug("Converting coordinates to J2000")
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

            # Modify the WCS solution to properly represents the WCS solution of binned data
            # print "Correcting the WCS solution for binning factor",binning,"..."
            for hdr_name in ('CRPIX1', 'CRPIX2'):
                wcs_header[hdr_name] /= binning
            for hdr_name in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
                wcs_header[hdr_name] *= binning

            reduction_files_used['wcs'] = options['wcs_distortion']

            cards = wcs_header.cards
            for (keyword, value, comment) in cards:
                if (keyword == 'CRVAL1'):
                    hdu.header["CRVAL1"] -= wcs_header['CRVAL1'] / math.cos(math.radians(hdu.header['CRVAL2']))
                elif (keyword == "CRVAL2"):
                    hdu.header['CRVAL2'] -= wcs_header['CRVAL2']
                else:
                    hdu.header[keyword] = (value, comment)

            # Make sure to write RAs that are positive
            if (hdu.header["CRVAL1"] < 0):
                hdu.header['CRVAL1'] += 360.
                
        # Insert the new image data. This also makes sure that the headers
        # NAXIS, NAXIS1, NAXIS2 are set correctly
        hdu.data = merged

        #
        # If requested, perform cosmic ray rejection
        #
        if (options['crj'] > 0):
            logger.debug("Starting cosmic ray removal with %d iterations" % (options['crj']))
            corrected, mask = podi_cosmicrays.remove_cosmics(hdu=hdu, 
                                                             n_iterations=options['crj'])
            hdu.data = corrected
            logger.debug("Done with cosmic ray removal")

        hdu.header['CRJ_ITER'] = (options['crj'], "cosmic ray removal iterations")


        source_cat = None
        if (options['fixwcs']):
            # Create source catalog
            
            if (options['skip_guide_ota'] and hdu.header['CELLMODE'].find("V") > -1):
                source_cat = None
            elif (False):
                source_cat = podi_findstars.find_stars(hdu, binning=4, boxsize=24, dumpfile=None, verbose=False,
                                                   detect_threshold=1., detect_minarea=4, roundness_limit=[-0.2,+0.2],
                                                   max_starcount=150,
                                                   extension_id=ota)
            else:
                logger.debug("Running SourceExtractor")
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
                sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
                    sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
                    fitsfile)
                if (options['verbose']): print sexcmd

                start_time = time.time()
                try:
                    ret = subprocess.Popen(sexcmd.split(), 
                                           stdout=subprocess.PIPE, 
                                           stderr=subprocess.PIPE)
                    (sex_stdout, sex_stderr) = ret.communicate()
                    #os.system(sexcmd)
                    if (ret.returncode != 0):
                        logger.warning("Sextractor might have a problem, check the log")
                        logger.debug("Stdout=\n"+sex_stdout)
                        logger.debug("Stderr=\n"+sex_stderr)
                except OSError as e:
                    podi_logging.log_exception()
                    print >>sys.stderr, "Execution failed:", e
                end_time = time.time()
                logger.debug("SourceExtractor returned after %.3f seconds" % (end_time - start_time))

                try:
                    try:
                        source_cat = numpy.loadtxt(catfile)
                    except IOError:
                        print "The Sextractor catalog is empty, ignoring this OTA"
                        source_cat = None
                    else:
                        if (source_cat.shape[0] == 0 or source_cat.ndim < 2):
                            source_cat = None
                        else:
                            source_cat[:,8] = ota
                            flags = source_cat[:,7]
                            no_flags = (flags == 0)
                            logger.debug("Found %d sources, %d with no flags" % (source_cat.shape[0], numpy.sum(no_flags)))
                except:
                    source_cat = None
                if (sitesetup.sex_delete_tmps):
                    clobberfile(fitsfile)
                    clobberfile(catfile)

            fixwcs_data = None

        #
        # Sample that background at random place so we can derive a median background level later on
        # This is necessary both for the pupil ghost correction AND the fringing correction
        #
        starcat = None
        if (source_cat != None):
            ota_x, ota_y = source_cat[:,2], source_cat[:,3]
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
            hdu.header["SKY_MEDI"] = (sky_level_median, "sky-level median")
            hdu.header["SKY_MEAN"] = (sky_level_mean, "sky-level mean")
            hdu.header["SKY_STD"] = (sky_level_std, "sky-level rms")


    data_products['hdu'] = hdu
    data_products['wcsdata'] = None #fixwcs_data
    data_products['sky-samples'] = sky_samples
    data_products['sky'] = (sky_level_median, sky_level_mean, sky_level_std)
    data_products['fringe_scaling'] = fringe_scaling
    data_products['sourcecat'] = source_cat
    data_products['tech-header'] = tech_header

    data_products['reduction_files_used'] = reduction_files_used

    return data_products #hdu, fixwcs_data
    













def collect_reduction_files_used(masterlist, files_this_frame):
    """
    Keeps track of all files used during the reduction. This function maintains 
    a list of files and their corresponding reduction step, and properly handles
    individual files as well as list of filenames.
    
    This routine is called during reduction to keep track of all file 
    associations.

    Parameters
    ----------
    masterlist : dictionary

        Dictionary of all reduction steps that have external files associated. 
        For each step it maintains a list of files.

    files_this_frame : dictionary

        Like the master_list, just for only one step. 


    Returns
    -------
    the updated master_list

    """

    for key, value in files_this_frame.iteritems():
        if (key in masterlist):
            existing_keys = masterlist[key]
            if (type(value) == list):
                for val1 in value:
                    masterlist[key].append(val1)
            else:
                masterlist[key].append(value)
        else:
            # This is a new key, so just copy it
            if (type(value) == list):
                masterlist[key] = value
            else:
                masterlist[key] = [value]

    # print masterlist
    return masterlist










def create_association_table(master, verbose=False):
    """

    Convert the association dictionary maintained and updated by 
    :proc:collect_reduction_files_used and creates a FITS-compatible 
    associations table that will be stored in the resulting output FITS file.

    Parameters
    ----------
    master : dictionary

        the master associations dictionary 

    verbose : Bool

        Activate some debugging output

    Returns
    -------
    tbhdu : TableHDU

        A FITS TableHDU containing the assocations table. Each entry in this 
        table contains the reduction step, the name of the associated file, and 
        the full path of this file.

    """

    reduction_step = []
    full_filename = []
    short_filename = []

    for key, value in master.iteritems():
        # print key,":",value
        for filename in set(value):
            reduction_step.append(key)
            if (filename == None):
                continue
            full_filename.append(os.path.abspath(filename))
            
            dirname, filebase = os.path.split(filename)
            short_filename.append(filebase)

            if (verbose):
                print "% 15s : %s" % (key, filename)

    # print reduction_step
    # print short_filename

    columns = [\
        pyfits.Column(name='correction',    format='A25',  array=reduction_step),
        pyfits.Column(name='filename_full', format='A375', array=full_filename),
        pyfits.Column(name='filename',      format='A100', array=short_filename),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("ASSOCIATIONS", comment=None)

    return tbhdu












def parallel_collect_reduce_ota(queue, return_queue,
                                options=None):
    """
    A minimal wrapper handling the parallel execution of collectcells.

    Parameters
    ----------
    queue : Queue

        Input Queue containing names of raw OTA frames to be reduced

    return_queue : Queue

        Queue to report the reduced data back to the main process

    options : dictionary

        containing all reduction parameters and settings

    Returns
    -------
    no return values, all returns are handled via the return_queue

    """

    # Setup the multi-processing-safe logging
    podi_logging.podi_logger_setup(options['log_setup'])

    while (True):
        cmd_quit, filename, ota_id = queue.get()
        if (cmd_quit):
            queue.task_done()
            return

        # Do the work
        try:
            data_products = collect_reduce_ota(filename, options=options)
        except (KeyboardInterrupt, SystemExit):
            queue.task_done()
            while (not queue.empty()):
                queue.get()
                queue.task_done()
            break
        except:
            # All other problems:
            # Note them in the debug log, and raise the error for further 
            # processing
            podi_logging.log_exception("parallel_collectcells")
            raise

        # Add the results to the return_queue so the master process can assemble the result file
        # print "Adding results for OTA",ota_id,"to return queue"
        # return_queue.put( (hdu, ota_id, wcsfix_data) )
        return_queue.put( (ota_id, data_products) )
        queue.task_done()
        
    return












def kill_all_child_processes(process_tracker):
    """
    Small function to clean up after collectcells timeouts. It receives a list 
    of processes IDs of all child processes initiated during execution and kills
    them one after the other.

    Parameters
    ----------
    process_tracker : Queue

        A queue of process-IDs of all child processes

    """

    try:
        while (True):
            child_pid = process_tracker.get(True, 1.)
            stdout_write("\n   terminating child process (process-ID %d) ..." % (child_pid))
            try:
                os.kill(child_pid, signal.SIGKILL)
            except:
                pass
            stdout_write(" done!")
    except Queue.Empty:
        stdout_write("\n   all child processes terminated!\n")

    return










def collectcells_with_timeout(input, outputfile,
                              batchmode=False,
                              verbose=False,
                              options=None,
                              showsplash=True,
                              timeout=None,
                              process_tracker=None):

    """ 
    Minimal wrapper to enable a timeout feature for collectcells. The actual 
    collectcells is started as a sub-process that can be joined for a specified
    time. If this time is up but the process is still running, it has exceeded
    its lifetime and will be terminated.

    Parameters
    ----------
    timeout : float

        Allowed maximum execution time before the task will be terminated

    process_tracker : Queue

        Forwarded from the main process, this queue will contain process IDs of 
        all sub-processes to allow for cleaning up process children after the
        timeout

    See collectcells for information about parameters.
    
    """

    # Start a collectcells as subprocess
    
    cc_args = (input, outputfile,
               process_tracker,
               batchmode,
               verbose,
               options,
               showsplash)

    p = multiprocessing.Process(target=collectcells, args=cc_args)
    if (verbose): print "Starting collectcells with timeout"
    p.start()
    if (verbose): print "collectcells started!"

    timeout = timeout if timeout > 0 else None
    p.join(timeout)
    if (p.is_alive()):
        stdout_write("\n\nTimeout event triggered, shutting things down ...")
        kill_all_child_processes(process_tracker)

        stdout_write("\n   Killing collectcells after timeout...")
        p.terminate()
        stdout_write(" all done!\n\n")

    return








def format_filename(input_filename_or_header, outputfile):
    """

    Format the input string and replace special keywords with information from
    the FITS header.

    """

    if (type(input_filename_or_header) == str):
        # This is a string, interpret it as a filename
        if (os.path.isfile(input_filename_or_header)):
            filename = input_filename_or_header
        elif (os.path.isdir(input_filename_or_header)):
            dirname = input_filename_or_header
            if (dirname.endswith("/")): dirname = dirname[:-1]
            # For now, assume that OTA 33 always exists, and that there is no .fz ending
            directory, obsid = os.path.split(dirname)
            filename = "%s/%s.33.fits" % (dirname, obsid)
            if (not os.path.isfile(filename)):
                filename = filename+".fz"
                if (not os.path.isfile(filename)):
                    return outputfile

        hdulist = pyfits.open(filename)
        header = hdulist[0].header
    else:
        # This is a header
        header = input_filename_or_header

    while (outputfile.find("%") >= 0):
        start = outputfile.find("%") 
        if (outputfile[start:start+7] == "%FILTER"):
            outputfile = outputfile[:start] + header['FILTER'] + outputfile[start+7:]
        elif (outputfile[start:start+7] == "%OBJECT"):
            # The object name might contain spaces, replace them with underscores
            # Also replace all other special characters with underscores
            objectname = header['OBJECT'].replace(' ', '_')
            objectname = objectname.replace(',', '_')
            objectname = objectname.replace('(', '_')
            objectname = objectname.replace(')', '_')
            objectname = objectname.replace('/', '_')
            objectname = objectname.replace('\\', '_')
            objectname = objectname.replace('`', '_')
            objectname = objectname.replace('"', '_')
            objectname = objectname.replace('\'', '_')
            objectname = objectname.replace(')', '_')
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

    return outputfile





def collectcells(input, outputfile,
                 process_tracker,
                 batchmode=False,
                 verbose=False,
                 options=None,
                 showsplash=True,
                 ):

    """
    collectcells handles all filename operations, ensuring all required files 
    exist, and hands the work on each OTA off to the suite of worker processes. 
    Finally assembles all results and writes the output-file.

    The actual basic reduction steps are done in collect_reduce_ota, but 
    collectcells take the output of collect_reduce_ota and finishes the 
    reduction, creating diagnostic plots, as well as performing reduction steps
    that require data from multiple OTAs, such as for the astrometric or 
    photometric calibration, removal of the pupil ghost, fringe removal with a 
    global scaling factor, etc.

    Parameters
    ----------
    input : string

        Either the name of a directory or filename that allows to construct the 
        name of all OTA FITS files of this exposure.

    outputfile : string

        Name of the output file

    batchmode : Bool

    verbose : Bool

        write extra progress updates to console

    options : dictionary

        Contains all reduction parameters and settings.

    showsplash : Bool

        Write a small splash screen with the name of the program, copyright
        notice and author information to console before starting the actual work.

    Returns
    -------
    ImageHDU

        If batchmode is set to True, collectcells will return the output HDU in 
        memory without writing it to disk, enabling post-processing without
        I/O overhead.

    """

    if (options == None): options = set_default_options()

    logger = logging.getLogger('CollectCells')
    logger.debug("Starting --collectcells--")

    if (showsplash):
        splash = """\

    **********************************************************************
    * This is podi_collectcells                                          *
    * (c) 2012-2013: Ralf Kotulla, kotulla@uwm.edu                       *
    *                University of Wisconsin (Milwaukee & Madison)       *
    *                WIYN Observatory, Inc                               *
    *                                                                    *
    * Please acknowledge the author when using any products generated    *
    * with this tool. For comments, questions or ideas for improvement   *
    * please send an email to kotulla@uwm.edu. Thank you!                *
    **********************************************************************

"""
        stdout_write(splash)

    # print "Received options:", options

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

        filebase = hdulist[0].header['FILENAME'][:-8]
        hdulist.close()
        del hdulist

        # Split the input filename to extract the directory part
        directory, dummy = os.path.split(input)
        if (directory == None or directory == ''):
            directory = "."

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
        stdout_write("Something is wrong here, can't find/open any of the files...")
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
                # Also replace all other special characters with underscores
                objectname = header['OBJECT'].replace(' ', '_')
                objectname = objectname.replace(',', '_')
                objectname = objectname.replace('(', '_')
                objectname = objectname.replace(')', '_')
                objectname = objectname.replace('/', '_')
                objectname = objectname.replace('\\', '_')
                objectname = objectname.replace('`', '_')
                objectname = objectname.replace('"', '_')
                objectname = objectname.replace('\'', '_')
                objectname = objectname.replace(')', '_')
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

    input_header = hdulist[0].header

    binning = get_binning(hdulist[0].header)

    # Check if the output file contains a new directory. 
    chk_directory, chk_filename = os.path.split(outputfile)
    changed_outputfilename = False
    if (chk_directory != ''):
        # The output file does contain a directory
        if (not os.path.isdir(chk_directory)):
            # if the directory doesn't exist ey, create it
            stdout_write("Output directory does not exist, creating it...\n")
            os.makedirs(chk_directory)
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
    ota_list[0].header["PLVER"] = (pipeline_plver, "name and version")
    ota_list[0].header["PIPELINE"] = (pipeline_name, "pipeline name")
    ota_list[0].header["PLVERSIO"] = (pipeline_version, "pipeline version")
    ota_list[0].header["PLAUTHOR"] = ("Ralf Kotulla", "pipeline author")
    ota_list[0].header["PLEMAIL"] = ("kotulla@uwm.edu", "contact email")
    for key, value in options['additional_fits_headers'].iteritems():
        ota_list[0].header[key] = (value, "user-added keyword")

    ota_list[0].header['BINNING'] = (binning, "binning factor")

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

        if (ota in options['skip_otas']):
            continue

        #print "Commanding work for extension",ota
        list_of_otas_being_reduced.append(list_of_otas_to_collect[ota_id])

        queue.put( (False, filename, ota_id+1) )

    logger.info("Performing instrumental detrending")
    # Create all processes to handle the actual reduction and combination
    #print "Creating",number_cpus,"worker processes"
    if ('profile' in options or number_cpus == 1):
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
            if (not process_tracker == None):
                if (verbose): print "Adding current slave process to process-tracker...", i
                process_tracker.put(p.pid)
                if (verbose): print "done adding to process-tracker"

            time.sleep(0.01)

        # Tell all workers to shut down when no more data is left to work on
        for i in range(len(processes)):
            if (verbose): stdout_write("Sending quit command!\n")
            queue.put((True,None,None))


    ############################################################
    #
    # Here is a good point to do some processing that does not depend on
    # individual OTAs. All slave-processes are busy, so the main process 
    # has some free time at hand before the first OTAs return
    #
    ############################################################

    #
    # Create a master list that keeps track of all additional files used for this frame
    #
    master_reduction_files_used = {}

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

    #
    # Compute some additional timing keywords for the mid- and end-point
    # of the exposure (important for observations of moving things)
    #

    # Read the date string from the header, making sure to trim off the 
    # fractions of a second to make ot python compatible
    date_format = "%Y-%m-%dT%H:%M:%S.%f"
    time_format = "%H:%M:%S.%f"
    date_obs = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], date_format)
    date_mid = date_obs + datetime.timedelta(seconds=(0.5*hdulist[0].header['EXPMEAS']))
    date_end = date_obs + datetime.timedelta(seconds=hdulist[0].header['EXPMEAS'])
    mjd = hdulist[0].header['MJD-OBS']

    #
    # write the mid-exposure and end-exopsure times back to FITS header
    #
    ota_list[0].header['DATE-MID'] = (datetime.datetime.strftime(date_mid, date_format)[:-3], 
                                      "Date at exposure mid-point")
    ota_list[0].header['TIME-MID'] = (datetime.datetime.strftime(date_mid, time_format)[:-3], 
                                      "Time at exposure mid-point")
    ota_list[0].header['DATE-END'] = (datetime.datetime.strftime(date_end, date_format)[:-3], 
                                      "Date at exposure end-point")
    ota_list[0].header['TIME-END'] = (datetime.datetime.strftime(date_end, time_format)[:-3], 
                                      "Time at exposure end-point")
    ota_list[0].header['MJD-MID'] = (mjd + (0.5 * hdulist[0].header['EXPMEAS'])/86400., 
                                     "MJD at exposure mid-point")
    ota_list[0].header['MJD-END'] = (mjd + (hdulist[0].header['EXPMEAS']/86400.), 
                                     "MJD at exposure end-point")

    # We know enough about the current frame, so close the file
    hdulist.close()
    del hdulist

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
    global_gain_sum, global_gain_count = 0, 0
    all_tech_headers = []

    for i in range(len(list_of_otas_being_reduced)):
        #hdu, ota_id, wcsfix_data = return_queue.get()
        try:
            ota_id, data_products = return_queue.get()
        except (KeyboardInterrupt, SystemExit):
            while (not return_queue.empty()):
                return_queue.get()
            raise
            return


        hdu = data_products['hdu']
        wcsfix_data = data_products['wcsdata']

        if (hdu == None):
            continue

        if ('reduction_files_used' in data_products):
            files_this_frame = data_products['reduction_files_used']
            # print "\n\n\n\n\nfiles_this_frame=\n\n",files_this_frame,"\n\n\n"
            collect_reduction_files_used(master_reduction_files_used, files_this_frame)

        global_gain_sum += (hdu.header['GAIN'] * hdu.header['NGAIN'])
        global_gain_count += hdu.header['NGAIN']

        ota_list[ota_id] = hdu
        
        sky_samples[hdu.header['EXTNAME']] = data_products['sky-samples']

        if (not data_products['fringe_scaling'] == None):
            #print data_products['fringe_scaling'].shape
            #if (fringe_scaling != None): print fringe_scaling.shape
            fringe_scaling = data_products['fringe_scaling'] if fringe_scaling == None else \
                numpy.append(fringe_scaling, data_products['fringe_scaling'], axis=0)

        if (not data_products['sourcecat'] == None):
            global_source_cat = data_products['sourcecat'] if (global_source_cat == None) \
                else numpy.append(global_source_cat, data_products['sourcecat'], axis=0)

        if ('tech-header' in data_products):
            all_tech_headers.append(data_products['tech-header'])
            # print "techdata for ota",ota_id,"\n",data_products['tech-header']

    #
    # Update the global gain variables
    #
    ota_list[0].header['GAIN'] = (global_gain_sum / global_gain_count)
    ota_list[0].header['NGAIN'] = global_gain_count

    # Now all processes have returned their results, terminate them 
    # and delete all instances to free up memory
    for cur_process in processes:
        cur_process.join()
        #cur_process.terminate()
        #del cur_process

    logger.debug("all data received from worker processes!")
    logger.info("Starting post-processing")
    additional_reduction_files = {}


    if (options['fixwcs'] and verbose):
        print fixwcs_extension
        print fixwcs_odi_x
        print fixwcs_odi_y
        print fixwcs_bestguess.shape
        print fixwcs_bestguess
        
    if(verbose):
        print master_reduction_files_used
        
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

        logger.debug("Updating header for extension %s..." % (ota.header['EXTNAME']))

        for header in headers_to_inherit:
            # Make sure the header we are about to move exists in the first place
            if (not header in ota.header):
                continue

            # Check if the header already exists in the primary header. If not add it!
            if (not header in ota_list[0].header):
                (keyword, value, comment) = ota.header.cards[header]
                ota_list[0].header[keyword] = (value, comment)
            
            # By now the value should exist in the primary header, 
            # so delete it from each of the extensions
            del ota.header[header]
                
        # Set the inherit keyword so that the headers removed from each 
        # extension are instead inherited from the primary
        ota.header["INHERIT"] = (True, "Inherit headers from PrimaryHDU")

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
    logger.debug("Combining sky-samples from all OTAs into global sky value")
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
    ota_list[0].header["SKYLEVEL"] = sky_global_median
    ota_list[0].header["SKYBG"] = sky_global_median
    ota_list[0].header["SKYNOISE"] = math.sqrt( 8.**2 + (sky_global_median*ota_list[0].header['GAIN'])**2 )
    logger.debug("Found global median sky-value = %.1f" % (sky_global_median))

    #
    # Now that we have the global sky-level, subtract the 
    # contribution of the pupil ghost to the science frame.
    # Different from the case of the calibration frames, use the radial
    # profile here, and ignore rotation (not needed anyway) to speed things up.
    #
    if (options['pupilghost_dir'] != None):
        logger.debug("Getting ready to subtract pupil ghost from science frame")
        filter_level = get_filter_level(ota_list[0].header)
        filter_name = get_valid_filter_name(ota_list[0].header)
        binning = ota_list[0].header['BINNING']
#        pg_template = "%s/pupilghost_radial___level_%d__bin%d.fits" % (options['pupilghost_dir'], filter_level, binning)
        pg_template = "%s/pupilghost_template___level_%d__bin%d.fits" % (options['pupilghost_dir'], filter_level, binning)
        logger.debug("looking for radial pupil ghost template %s...\n" % (pg_template))
        # If we have a template for this level
        if (os.path.isfile(pg_template)):
            logger.debug("\n   Using pupilghost template %s, filter %s ... " % (pg_template, filter_name))
            pg_hdu = pyfits.open(pg_template)

            # Find the optimal scaling factor
            logger.debug("Searching for optimal pupilghost scaling factor")
            any_affected, scaling, scaling_std = podi_matchpupilghost.get_pupilghost_scaling(ota_list, pg_hdu)
            logger.debug("Check if any OTA is affected: %s" % ("yes" if any_affected else "no"))
            logger.debug("Optimal scaling factor found: %.2f +/- %.2f" % (scaling, scaling_std))

            if (any_affected):
                # And subtract the scaled pupilghost templates.
                logger.debug("Commencing pupilghost subtraction")
                podi_matchpupilghost.subtract_pupilghost(ota_list, pg_hdu, scaling, 
                                                         # rotate=False,
                                                         rotate=True,
                                                         source_center_coords='header')

                ota_list[0].header["PUPLGOST"] = (pg_template, "p.g. template")
                ota_list[0].header["PUPLGFAC"] = (scaling, "pupilghost scaling")
                bg_scaled = podi_matchpupilghost.scaling_factors[filter_name]*sky_global_median
                ota_list[0].header["PUPLGFA2"] = (bg_scaled, "analytical pupilghost scaling")
                stdout_write(" done!\n")
                additional_reduction_files['pupilghost'] = pg_template

            pg_hdu.close()
        else:
            logger.info("Pupilghost correction requested, but no template found:")
            logger.info("  was looking for filename %s" % (pg_template))

    #
    # Fix the WCS if requested
    # New WCS matching using CCMatch
    #
    ota_list[0].header['WCSFIXED'] = False
    ota_list[0].header['WCSCAL'] = (False, "Was astrometric solution found")
    ota_list[0].header['WCSXRMS'] = (-1., "RMS in x-dir of astrometric solution")
    ota_list[0].header['WCSYRMS'] = (-1., "RMS in y-dir of astrometric solution")
    ota_list[0].header['ASTRMCAT'] = ("", "astrometric reference catalog")

    enough_stars_for_fixwcs = not global_source_cat == None \
                              and global_source_cat.shape[0]>3
    if (options['fixwcs']):
        if (enough_stars_for_fixwcs):
            logger.debug("Found enough stars for astrometric calibration")
        else:
            logger.info("Couldn't find enough stars for astrometric calibration!")

    if (options['fixwcs'] 
        and enough_stars_for_fixwcs):

        logger.info("Performing astrometric calibration")
        # The entire source catalog is located in: --> global_source_cat
        # Current HDUList is in: --> ota_list

        import dev_ccmatch
        if (dev_ccmatch.create_debug_files): numpy.savetxt("debug.wcs_raw", global_source_cat)
        ccmatched = dev_ccmatch.ccmatch(source_catalog=global_source_cat,
                                       reference_catalog=None, # meaning ccmtch will obtain it
                                       input_hdu=ota_list, 
                                       mode="otashear")

        # Use the fixed HDUList
        ota_list = ccmatched['hdulist']
        ota_list[0].header['WCSFIXED'] = True
        ota_list[0].header['ASTRMCAT'] = "2MASS"
        if ("WCS_QUAL" in ota_list[0].header):
            ota_list[0].header['WCSCAL'] = ota_list[0].header['WCS_QUAL'] > 1.5
        
        # Also extract some of the matched/calibrated catalogs
        odi_2mass_matched = ccmatched['matched_src+2mass']
        global_source_cat = ccmatched['calibrated_src_cat']

        # Append the 2MASS reference catalog to output frame
        logger.debug("Creating a FITS table for the full 2MASS reference catalog")
        twomass_hdu = twomasscat_to_tablehdu(ccmatched['2mass-catalog']) #fixwcs_ref_ra, fixwcs_ref_dec)
        ota_list.append(twomass_hdu)

        logger.debug("Creating a FITS table for the full ODI catalog")
        src_tbhdu = odi_sources_to_tablehdu(ccmatched['calibrated_src_cat'])
        #logger.debug(src_tbhdu)
        ota_list.append(src_tbhdu)


        #
        # Also create a TableHDU for the matched ODI+2MASS catalog
        # 
        logger.debug("Creating a FITS table for the matched ODI+2MASS catalog")
        odi_2mass_cat = ccmatched['matched_src+2mass']
        columns = [pyfits.Column(name='ODI_RA', format='D', unit='degrees', 
                                 array=odi_2mass_cat[:,  0], disp='ODI right ascension'),
                   pyfits.Column(name='ODI_DEC', format='D', unit='degrees', 
                                 array=odi_2mass_cat[:,  1], disp='ODI declination'),
                   pyfits.Column(name='TWOMASS_RA', format='D', unit='degrees', 
                                 array=odi_2mass_cat[:, -2], disp='2MASS right ascension'),
                   pyfits.Column(name='TWOMASS_DEC', format='D', unit='degrees', 
                                 array=odi_2mass_cat[:, -1], disp='2MASS declination'),
                   pyfits.Column(name='OTA', format='E', unit='',
                                 array=odi_2mass_cat[:, 8], disp='source OTA'),
        ]
        coldefs = pyfits.ColDefs(columns)
        matchedhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
        matchedhdu.update_ext_name("CAT.ODI+2MASS", comment=None)
        matchedhdu.header['MATCHRAD'] = (2., "matching radius in arcsec")
        ota_list.append(matchedhdu)

        # Compute the WCS quality statistics
        # This also writes to the Primary and OTA-level headers
        wcs_quality = dev_ccmatch.global_wcs_quality(odi_2mass_cat, ota_list)
        # print "WCS quality =",wcs_quality
        ota_list[0].header['WCSXRMS'] = wcs_quality['full']['RMS-RA']
        ota_list[0].header['WCSYRMS'] = wcs_quality['full']['RMS-DEC']

        # Compute the image quality using all detected sources
        # Make sure to only include round source (elongation < 1.3) and decent 
        # photometry (i.e. high enoug S/N)
        good_fwhm_values = (global_source_cat[:, SXcolumn['flags']] == 0) & \
                           (global_source_cat[:, SXcolumn['elongation']] < 1.3) & \
                           (global_source_cat[:, SXcolumn['mag_err_3.0']] <= 0.2)
        seeing = global_source_cat[:, SXcolumn['fwhm_world']] * 3600. # convert to arcsec
        seeing_all = numpy.median(seeing)
        seeing_clipped = three_sigma_clip(seeing)
        seeing_filtered = numpy.median(seeing_clipped)
        ota_list[0].header['FWHM_ALL'] = (seeing_all, "median FWHM of all sources")
        ota_list[0].header['FWHM_FLT'] = (seeing_filtered, "median FWHM, 3sigma clipped")
        ota_list[0].header['FWHMNFLT'] = (seeing_clipped.shape[0], "number of src in FWHM comp")

        ota_list[0].header['FWHMSTAR'] = (seeing_all, "median FWHM of SDSS-matched stars")
        ota_list[0].header['SEEING'] = (seeing_all, "Seeing [arcsec]")
        ota_list[0].header['SEEING_N'] = (seeing_clipped.shape[0], "number of stars in seeing comp")

        
    if (options['fixwcs'] 
        and options['create_qaplots']
        and enough_stars_for_fixwcs):
        ota_outlines = derive_ota_outlines(ota_list)
            # print "Creating some diagnostic plots"

        diagnostic_plot_title = "%s\n(obsid: %s - filter: %s- exptime: %ds)" % (
            ota_list[0].header['OBJECT'],
            ota_list[0].header['OBSID'],
            ota_list[0].header['FILTER'],
            int(ota_list[0].header['EXPTIME']),
            )

        ota_wcs_stats = wcs_quality #None #{}
        # for extension in range(1, len(ota_list)):
        #     if (ota_list[extension] == None):
        #         continue

        #     # Compute the typical RMS of the alignment
        #     ota = int(ota_list[extension].header['EXTNAME'][3:5])
        #     results = podi_fixwcs.compute_wcs_quality(odi_2mass_matched, ota, ota_list[extension].header)
        #     # print results
        #     ota_wcs_stats[ota_list[extension].header['EXTNAME']] = results
        # results = podi_fixwcs.compute_wcs_quality(odi_2mass_matched, None, ota_list[0].header)
        # # print results
        # ota_wcs_stats['full'] = results
        
        import podi_diagnosticplots
        title_info = ota_list[0].header.copy()

        # Create the WCS scatter plot
        plotfilename = create_qa_filename(outputfile, "wcs1", options)
        podi_diagnosticplots.wcsdiag_scatter(matched_radec_odi=odi_2mass_matched[:,0:2], 
                                             matched_radec_2mass=odi_2mass_matched[:,-2:],
                                             matched_ota=odi_2mass_matched[:,SXcolumn['ota']],
                                             filename=plotfilename, 
                                             options=options,
                                             ota_wcs_stats=ota_wcs_stats,
                                             also_plot_singleOTAs=options['otalevelplots'],
                                             title_info=title_info)

        # Create the WCS shift plot
        plotfilename = create_qa_filename(outputfile, "wcs2", options)
        podi_diagnosticplots.wcsdiag_shift(matched_radec_odi=odi_2mass_matched[:,0:2],
                                           matched_radec_2mass=odi_2mass_matched[:,-2:],
                                           matched_ota=odi_2mass_matched[:,SXcolumn['ota']],
                                           filename=plotfilename, #outputfile[:-5]+".wcs2", 
                                           options=options,
                                           ota_wcs_stats=ota_wcs_stats,
                                           ota_outlines=ota_outlines,
                                           also_plot_singleOTAs=options['otalevelplots'],
                                           title_info=title_info)

    #
    # Add some default photometric calibration keywords
    #
    ota_list[0].header['PHOTMCAT'] = (None, "catalog used for photometric calibration")
    
    ota_list[0].header["PHOTZP"]   = (-99., "phot. zeropoint corr for exptime")
    ota_list[0].header["PHOTZPSD"] = (-99., "zeropoint std.dev.")
    ota_list[0].header["PHOTZP_X"] = (-99., "phot zeropoint for this frame")
    ota_list[0].header["PHOTZPSP"] = (-99., "phot ZP upper 1sigma limit")
    ota_list[0].header["PHOTZPSM"] = (-99., "phot ZP lower 1sigma limit")
    ota_list[0].header["PHOTZPER"] = (-99., "phot ZP std.err of the mean")
    ota_list[0].header["PHOTZP_N"] = (-1, "number stars in clipped distrib.")
    ota_list[0].header["PHOTZPN0"] = (-1, "total number of matched ref stars")

    # AuCaP compatible photometric zeropoint
    ota_list[0].header["MAGZERO"]  = (-99, "phot. zeropoint corr for exptime")
    ota_list[0].header["MAGZSIG"]  = (-99, "phot ZP dispersion")
    ota_list[0].header["MAGZERR"]  = (-99, "phot ZP uncertainty")

    # Add some information on what apertures were used for the photometric calibration
    ota_list[0].header['MAG0MODE'] = ("none", "how was aperture determined")
    ota_list[0].header['MAG0SIZE'] = (-99, "what aperture size was used")
    ota_list[0].header['MAG0_MAG'] = ("none", "id string for magnitude")
    ota_list[0].header['MAG0_ERR'] = ("none", "is string for mag error")

    # Compute a flux scaling keyword to be used with swarp
    ota_list[0].header["FLXSCALE"] = (1.0, "flux scaling factor for ZP=25")

    # Add default values forthe radial ZP trend
    ota_list[0].header['RADZPFIT'] = (False, "was a radial ZP (ZP=p0 + p1 x r) fit done?")
    ota_list[0].header['RADZP_P0'] = (-99, "radial ZP fit, ZP offset at r=0")
    ota_list[0].header['RADZP_P1'] = (-99, "radial ZP fit, slope in mag/degree")
    ota_list[0].header['RADZP_E0'] = (-99, "radial ZP fit, ZP offset error")
    ota_list[0].header['RADZP_E1'] = (-99, "radial ZP fit, slope error")

    # Add keywords for the magnitude/count/error restricted zeropoint calculation
    ota_list[0].header['ZPRESMED'] = (-99., "phot ZP of restricted catalog")
    ota_list[0].header['ZPRESSTD'] = (-99., "phot. ZP std.dev. of rest catalog")
    ota_list[0].header['ZPRES_SP'] = (-99., "phot ZP upper 1 sigma of rest catalog")
    ota_list[0].header['ZPRES_SM'] = (-99., "phot ZP lower 1 sigma of rest catalog")
    ota_list[0].header['ZPRES__N'] = (  -1, "phot ZP restricted catalog size")
    ota_list[0].header['ZPRES_MD'] = ( -99, "phot ZP restricted cat. median ODI mag")
    ota_list[0].header['ZPRES_MX'] = ( -99, "phot ZP restricted cat. max ODI mag")
    ota_list[0].header['ZPRES_MN'] = ( -99, "phot ZP restricted cat. min ODI mag")

    # Add more keywords for the ZP-ODI_mag relation
    ota_list[0].header['ZPSLP_P0'] = (0., "phot ZP - magnitude slope - y0")
    ota_list[0].header['ZPSLP_P1'] = (0., "phot ZP - magnitude slope - y1")
    ota_list[0].header['ZPSLP_E0'] = (0., "phot ZP - magnitude slope - dy0")
    ota_list[0].header['ZPSLP_E1'] = (0., "phot ZP - magnitude slope - dy1")

    # Theoretical photometric ZP
    ota_list[0].header['MAGZREF'] = (-99, "reference photometric zeropoint")
    
    # ZP corrected for atmospheric extinction
    ota_list[0].header['MAGZ_AM1'] = (-99, "phot zeropoint corrected to airmass = 1.0")
            
    # color-term corrections
    ota_list[0].header['MAGZ_CT'] = (False, "was a color-term correction used")
    ota_list[0].header['MAGZ_COL'] = ("", "color used in color-term correction")
    ota_list[0].header['MAGZ_CTC'] = (-99, "slope of color-term")

    # sky-brightness during observation
    ota_list[0].header['SKYMAG'] = (-99, "sky brightness in mag/arcsec^2")

    #
    # If requested, perform photometric calibration
    #
    if (options['photcalib'] 
        and options['fixwcs']
        and enough_stars_for_fixwcs):
        zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = 99., 99., None, 99.
        logger.info("Starting photometric calibration")

        exptime = ota_list[0].header['EXPTIME']
        titlestring = "%s\n(obsid: %s - filter: %s- exptime: %ds)" % (
            ota_list[0].header['OBJECT'],
            ota_list[0].header['OBSID'],
            ota_list[0].header['FILTER'],
            int(ota_list[0].header['EXPTIME']),
            )
        filter_name = get_valid_filter_name(ota_list[0].header)

        photcalib_details = {}
        zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = \
            podi_photcalib.photcalib(global_source_cat, outputfile, filter_name, 
                                     exptime=exptime,
                                     diagplots=True,
                                     plottitle=titlestring,
                                     otalist=ota_list,
                                     options=options,
                                     detailed_return=photcalib_details)

        ota_list[0].header['PHOTMCAT'] = ("SDSS")

        ota_list[0].header["PHOTZP"] = (zeropoint_median, "phot. zeropoint corr for exptime")
        ota_list[0].header["PHOTZPSD"] = (zeropoint_std, "zeropoint std.dev.")
        ota_list[0].header["PHOTZP_X"] = (zeropoint_exptime, "phot zeropoint for this frame")
        ota_list[0].header["PHOTZPSP"] = (photcalib_details['zp_upper1sigma'], "phot ZP upper 1sigma limit")
        ota_list[0].header["PHOTZPSM"] = (photcalib_details['zp_lower1sigma'], "phot ZP lower 1sigma limit")
        ota_list[0].header["PHOTZPER"] = (photcalib_details['stderrofmean'], "phot ZP std.err of the mean")
        ota_list[0].header["PHOTZP_N"] = (photcalib_details['n_clipped'], "number stars in clipped distrib.")
        ota_list[0].header["PHOTZPN0"] = (photcalib_details['n_raw'], "total number of matched ref stars")

        ota_list[0].header["MAGZERO"] = (photcalib_details['median'], "phot. zeropoint corr for exptime")
        ota_list[0].header["MAGZSIG"] = (photcalib_details['std'], "phot ZP dispersion")
        ota_list[0].header["MAGZERR"] = (photcalib_details['stderrofmean'], "phot ZP uncertainty")
        flux_scaling = math.pow(10., -0.4*(zeropoint_exptime - 25.0))
        ota_list[0].header["FLXSCALE"] = flux_scaling
        # For swarp to work properly, also copy the FLXSCALE keyword into every image extension
        for i in range(len(ota_list)):
            if (is_image_extension(ota_list[i])):
                ota_list[i].header["FLXSCALE"] = flux_scaling

        # Add some information on what apertures were used for the photometric calibration
        ota_list[0].header['MAG0MODE'] = (photcalib_details['aperture_mode'], "how was aperture determined")
        ota_list[0].header['MAG0SIZE'] = (photcalib_details['aperture_size'], "what aperture size was used")
        ota_list[0].header['MAG0_MAG'] = (photcalib_details['aperture_mag'], "id string for magnitude")
        ota_list[0].header['MAG0_ERR'] = (photcalib_details['aperture_magerr'], "is string for mag error")

        if (not photcalib_details['radialZPfit'] == None):
            ota_list[0].header['RADZPFIT'] = True
            ota_list[0].header['RADZP_P0'] = photcalib_details['radialZPfit'][0]
            ota_list[0].header['RADZP_P1'] = photcalib_details['radialZPfit'][1]
            ota_list[0].header['RADZP_E0'] = photcalib_details['radialZPfit_error'][0]
            ota_list[0].header['RADZP_E1'] = photcalib_details['radialZPfit_error'][1]

        if (not photcalib_details['zp_restricted'] == None):
            (sel_median, sel_std, sel_psigma, sel_msigma, sel_n, sel_medodimag, sel_maxodimag, sel_minodimag) = photcalib_details['zp_restricted']
            ota_list[0].header['ZPRESMED'] = sel_median
            ota_list[0].header['ZPRESSTD'] = sel_std
            ota_list[0].header['ZPRES_SP'] = sel_psigma
            ota_list[0].header['ZPRES_SM'] = sel_msigma
            ota_list[0].header['ZPRES__N'] = sel_n
            ota_list[0].header['ZPRES_MD'] = sel_medodimag
            ota_list[0].header['ZPRES_MX'] = sel_maxodimag
            ota_list[0].header['ZPRES_MN'] = sel_minodimag
            
        if (not photcalib_details['zp_magnitude_slope'] == None):
            fit, uncert = photcalib_details['zp_magnitude_slope']
            ota_list[0].header['ZPSLP_P0'] = fit[0]
            ota_list[0].header['ZPSLP_P1'] = fit[1]
            ota_list[0].header['ZPSLP_E0'] = uncert[0]
            ota_list[0].header['ZPSLP_E1'] = uncert[1]

        ref_ZP = -99. if not filter_name in reference_zeropoint else reference_zeropoint[filter_name][0]
        ota_list[0].header['MAGZREF'] = (ref_ZP, "reference photometric zeropoint")

        # Also compute the zeropoint after correction for airmass
        zp_airmass1 = -99.
        if (filter_name in atm_extinction):
            zp_airmass1 = zeropoint_median + (ota_list[0].header['AIRMASS']-1) * atm_extinction[filter_name]
        ota_list[0].header['MAGZ_AM1'] = (zp_airmass1, "phot Zeropoint corrected for airmass")
            
        # Add some information whether or not we performed a color-term correction
        colorterm_correction = (not photcalib_details['colorterm'] == None)
        ota_list[0].header['MAGZ_CT'] = colorterm_correction
        ota_list[0].header['MAGZ_COL'] = photcalib_details['colorcorrection'] if colorterm_correction else ""
        ota_list[0].header['MAGZ_CTC'] = photcalib_details['colorterm'] if colorterm_correction else 0.0

        # Compute the sky-brightness 
        sky_arcsec = sky_global_median / (0.11**2) # convert counts/pixel to counts/arcsec*2
        sky_mag = -2.5 * math.log10(sky_arcsec) + zeropoint_exptime
        ota_list[0].header['SKYMAG'] = sky_mag



        # Now convert the matched source catalog into a valid FITS extension 
        # and add it to the output.
        if (not odi_sdss_matched == None and odi_sdss_matched.shape[0] > 0):
            match_tablehdu = create_odi_sdss_matched_tablehdu(odi_sdss_matched)
            # Copy a bunch of headers so we can makes heads and tails of the catalog
            # even if it's separated from the rest of the file.
            for hdrname in ['AIRMASS',
                            'FILTER', 'FILTERID',
                            'MJD-OBS',
                            'OBSID',
                            'EXPTIME', 'EXPMEAS',
                            'TELFOCUS',
                            'ROTSTART', 'ROTEND',
                            'ADCMODE',
                            'GAIN',
            ]:
                if (hdrname in ota_list[0].header):
                    match_tablehdu.header[hdrname] = ota_list[0].header[hdrname]
            ota_list.append(match_tablehdu)

            # Also use the matched catalog to determine the seeing of only stars
            star_seeing = odi_sdss_matched[:, SXcolumn['fwhm_world']+2] * 3600.
            cleaned = three_sigma_clip(star_seeing)
            seeing = numpy.median(cleaned)
            logger.debug("Seeing is %.2fdeg = %.2f arcsec" % (seeing, seeing*3600.))
            ota_list[0].header['FWHMSTAR'] = (seeing, "median FWHM of SDSS-matched stars")
            ota_list[0].header['SEEING'] = (seeing, "Seeing [arcsec]")
            ota_list[0].header['SEEING_N'] = (cleaned.shape[0], "number of stars in seeing comp")

            # Compute a approximate detection limit
            # This assumes an aperture with diameter of 2x the seeing 
            # (i.e. radius = seeing).
            aperture_area = (seeing / 0.11)**2 * math.pi
            readnoise = 8.
            bgcounts = aperture_area * (sky_global_median + readnoise**2)
            counts_sn1 = math.sqrt(bgcounts)
            limiting_mag = -2.5*math.log10(counts_sn1) + zeropoint_exptime
            ota_list[0].header['PHOTDPTH'] = limiting_mag

    # Create the image quality plot
    # Also create a diagnostic plot for the Seeing.
    # choose the sdssm-matched catalog if available, otherwise use the raw source catalog.
    if (options['fixwcs'] and options['create_qaplots'] and enough_stars_for_fixwcs):

        if (options['photcalib'] and not odi_sdss_matched == None and odi_sdss_matched.shape[0] > 0):
            # Use the SDSS catalog if available
            flags = odi_sdss_matched[:,SXcolumn['flags']+2]
            valid_flags = (flags == 0)
            ra = odi_sdss_matched[:,SXcolumn['ra']][valid_flags]
            dec= odi_sdss_matched[:,SXcolumn['dec']][valid_flags]
            fwhm = odi_sdss_matched[:,SXcolumn['fwhm_world']+2][valid_flags]
            ota = odi_sdss_matched[:,SXcolumn['ota']+2][valid_flags]

        else:
            # This should be cleaned up to make the call for this plot nicer
            flags = global_source_cat[:,SXcolumn['flags']]
            valid_flags = (flags == 0) & (global_source_cat[:,SXcolumn['mag_err_auto']] < 0.2)
            ra = global_source_cat[:,SXcolumn['ra']][valid_flags]
            dec= global_source_cat[:,SXcolumn['dec']][valid_flags]
            fwhm = global_source_cat[:,SXcolumn['fwhm_world']][valid_flags]
            ota = global_source_cat[:,SXcolumn['ota']][valid_flags]

        plotfilename = create_qa_filename(outputfile, "seeing", options)
        podi_diagnosticplots.diagplot_psfsize_map(ra, dec, fwhm, ota, 
                                                  output_filename=plotfilename, #outputfile[:-5]+".seeing",
                                                  title=diagnostic_plot_title,
                                                  ota_outlines=ota_outlines, 
                                                  options=options,
                                                  also_plot_singleOTAs=options['otalevelplots'],
                                                  title_info=title_info)


    if (not options['nonsidereal'] == None):
        logger.info("Starting non-sidereal WCS modification")
        apply_nonsidereal_correction(ota_list, options, logger)

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
        ota_list[0].header["FRNG_SCL"] = fringe_scaling_median
        ota_list[0].header["FRNG_STD"] = fringe_scaling_std

        # Construct the name of the fringe map
        filter_name = ota_list[0].header['FILTER']
        fringe_filename = check_filename_directory(options['fringe_dir'], "fringe__%s.fits" % (filter_name))
        # print fringe_filename
        if (os.path.isfile(fringe_filename)):
            additional_reduction_files['fringemap'] = fringe_filename
            fringe_hdulist = pyfits.open(fringe_filename)

            # Now do the correction
            for ext in range(1, len(ota_list)):
                extname = ota_list[ext].header['EXTNAME']
                if (type(ota_list[ext]) != pyfits.hdu.image.ImageHDU):
                    continue
                # print "subtracting",extname
                ota_list[ext].data -= (fringe_hdulist[extname].data * fringe_scaling_median)


    #
    # Add some information about the filter bandpasses to the output file
    #
    filtername = ota_list[0].header['FILTER']
    bandpass = [0.,0.] if (not filtername in filter_bandpass) else filter_bandpass[filtername]
    ota_list[0].header['PHOTCLAM'] = (bandpass[0], "central wavelength of filter [nm]")
    ota_list[0].header['PHOTBW'] = (bandpass[1], "RMS width of filter [nm]")
    ota_list[0].header['PHOTFWHM'] = (bandpass[1], "FWHM of filter [nm]")

    #
    # Prepare the Tech-HDU and add it to output HDU
    #
    techhdu = pyfits.ImageHDU(name='TECHDATA')
    # Now add all tech-data to the techhdu
    for techhdr in all_tech_headers:
        for (key, value, comment) in techhdr.cards:
            techhdu.header[key] = (value, comment)
    ota_list.append(techhdu)


    #
    # Create an association table from the master reduction files used.
    # 
    
    master_reduction_files_used = collect_reduction_files_used(master_reduction_files_used, 
                                                               additional_reduction_files)
    assoc_table = create_association_table(master_reduction_files_used)
    ota_list.append(assoc_table)


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


def apply_nonsidereal_correction(ota_list, options, logger=None):
    # This WCS in this frame should be corrected for non-sidereal tracking
    # Tracking rates are given in arcseconds per hour
    # Note that d_ra is given as dRA*cos(dec)

    if (logger == None):
        logger = logging.getLogger("NonsiderealCorr")

    if (options['nonsidereal'] == None):
        logger.debug("No nonsidereal option set, skipping non-sidereal correction")
        return
    if (not ('dra' in options['nonsidereal'] and
             'ddec' in options['nonsidereal'] and
             'ref_mjd' in options['nonsidereal'])):
        logger.debug("One of the nonsidereal option missing, skipping non-sidereal correction")
        logger.debug("Available options are:"+str(options['nonsidereal']))
        return
             
    mjd = ota_list[0].header['MJD-OBS']
    delta_t_hours = (mjd - options['nonsidereal']['ref_mjd']) * 24.
    dra_t = options['nonsidereal']['dra'] * delta_t_hours / 3600.
    ddec_t = options['nonsidereal']['ddec'] * delta_t_hours / 3600.

    ota_list[0].header['NSIDPMRA'] = options['nonsidereal']['dra']
    ota_list[0].header['NSIDPMDE'] = options['nonsidereal']['ddec']
    ota_list[0].header['NSIDBASE'] = options['nonsidereal']['ref_mjd']
    logger.debug("Tracking rates are dRA=%(dra)f dDEC=%(ddec)f arcsec/hour" % options['nonsidereal'])
    logger.debug("Time-offset to reference frame: %f hours" % (delta_t_hours))

    for ext in range(len(ota_list)):
        if (not is_image_extension(ota_list[ext])):
            continue
        declination = ota_list[ext].header['CRVAL2']
        dra_corrected = dra_t / math.cos(math.radians(declination))

        ota_list[ext].header['CRVAL1'] -= dra_corrected
        ota_list[ext].header['CRVAL2'] -= ddec_t

        ota_list[ext].header['NSIDDRA'] = dra_corrected
        ota_list[ext].header['NSIDDDEC'] = ddec_t
        logger.debug("Adding offset of %f %f arcsec to OTA %s" % (
            dra_corrected*3600., ddec_t*3600, ota_list[ext].header['EXTNAME'])
        )

    return

    
def create_odi_sdss_matched_tablehdu(odi_sdss_matched):

    # Define what columns are in the SDSS catalog
    SDSScolumn_names = [
        'ra', 'dec',
        'u', 'u_err',
        'g', 'g_err',
        'r', 'r_err',
        'i', 'i_err',
        'z', 'z_err',
        ]
    SDSScolumn = {}
    # Convert names into IDs, account for the number of columns
    # (in the matched catalog) already filled by the ODI catalog
    for name in SDSScolumn_names:
        SDSScolumn[name] = len(SDSScolumn) + len(SXcolumn)

    # Note:
    # The +2 in the collumn indices accounts for the fact that the Ra/Dec 
    # from the SDSS catalog are inserted after the Ra/Dec of the ODI 
    # catalog, shifting all other columns to the right by two positions.
    columns = [\
        #
        # Ra/Dec from ODI catalog
        #
        pyfits.Column(name='ODI_RA', disp='right ascension',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 0]),
        pyfits.Column(name='ODI_DEC', disp='right ascension',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 1]),

        #
        # Ra/Dec from SDSS
        #
        pyfits.Column(name='SDSS_RA', disp='right ascension',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 2]),
        pyfits.Column(name='SDSS_DEC', disp='right ascension',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 3]),


        #
        # Now add the SDSS magnitudes
        #
        pyfits.Column(name='SDSS_MAG_U', disp='SDSS magnitude u-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['u']]),
        pyfits.Column(name='SDSS_ERR_U', disp='SDSS magnitude error u-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['u_err']]),

        pyfits.Column(name='SDSS_MAG_G', disp='SDSS magnitude g-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['g']]),
        pyfits.Column(name='SDSS_ERR_G', disp='SDSS magnitude error g-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['g_err']]),

        pyfits.Column(name='SDSS_MAG_R', disp='SDSS magnitude r-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['r']]),
        pyfits.Column(name='SDSS_ERR_R', disp='SDSS magnitude error r-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['r_err']]),

        pyfits.Column(name='SDSS_MAG_I', disp='SDSS magnitude i-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['i']]),
        pyfits.Column(name='SDSS_ERR_I', disp='SDSS magnitude error i-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['i_err']]),

        pyfits.Column(name='SDSS_MAG_Z', disp='SDSS magnitude z-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['z']]),
        pyfits.Column(name='SDSS_ERR_Z', disp='SDSS magnitude error z-band',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['z_err']]),


        pyfits.Column(name='ODI_FWHM', disp='FWHM in ODI frame',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, SXcolumn['fwhm_world']+2]),

       pyfits.Column(name='ODI_MAG_AUTO', format='E', unit='mag',
                     array=odi_sdss_matched[:,SXcolumn['mag_auto']+2],
                     disp='auto-mag'),
       pyfits.Column(name='ODI_ERR_AUTO', format='E', unit='mag',
                     array=odi_sdss_matched[:,SXcolumn['mag_err_auto']+2],
                     disp='auto-mag error'),
        
        ]


    #
    # Magnitudes and errors for all ODI apertures
        #
    for i in range(len(SXapertures)):
        apsize = SXapertures[i]
        col_mag = "mag_aper_%0.1f" % (apsize)
        col_magerr = "mag_err_%0.1f" % (apsize)
        columns.append(pyfits.Column(name='ODI_MAG_D%02d' % (int(apsize*10.)), 
                                     format='E', unit='mag',
                                     array=odi_sdss_matched[:,SXcolumn[col_mag]+2],
                                     disp='aperture mag %0.1f arcsec diameter' % (apsize)
                                 )
        )
        columns.append(pyfits.Column(name='ODI_ERR_D%02d' % (int(apsize*10.)),
                                     format='E', unit='mag',
                                     array=odi_sdss_matched[:,SXcolumn[col_magerr]+2],
                                     disp='ap. mag error %0.1f arcsec diameter' % (apsize)
                                 )
        )


    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
    tbhdu.update_ext_name("CAT.PHOTCALIB", comment=None)

    return tbhdu





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



def twomasscat_to_tablehdu(catalog): 
    
    columns = [\
        pyfits.Column(name='RA',  format='D', array=catalog[:,0]),
        pyfits.Column(name='DEC', format='D', array=catalog[:,1]),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("CAT.2MASS", comment=None)
    return tbhdu











def odi_sources_to_tablehdu(source_cat):

    """
    Create a FITS table containing the source catalog created during reduction.

    Parameters
    ----------
    source_cat : numpy array

        Global ODI source catalog collected from the source catalog of the 
        individual OTA catalogs

    Returns
    -------
    tbhdu : pyfits.TableHDU

        Binary FITS table extension

    """
    
    #
    # Start with some general source information
    #
    columns = [\
               pyfits.Column(name='RA', format='D', unit='degrees',
                             array=source_cat[:,SXcolumn['ra']],
                             disp='right ascension'),
               pyfits.Column(name='DEC',            format='D', unit='degrees',
                             array=source_cat[:,SXcolumn['dec']], 
                             disp='declination'),
               pyfits.Column(name='X',              format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['x']],
                             disp='center x'),
               pyfits.Column(name='Y',              format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['y']], 
                             disp='center y'),
               pyfits.Column(name='FWHM_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['fwhm_image']],
                             disp='FWHM in pixels'),
               pyfits.Column(name='FWHM_WORLD',     format='E', unit='deg',
                             array=source_cat[:,SXcolumn['fwhm_world']],
                             disp='FWHM in degrees'),
               pyfits.Column(name='BACKGROUND',     format='E', unit='counts',
                             array=source_cat[:,SXcolumn['background']],
                             disp='background level'),
               pyfits.Column(name='FLAGS',          format='I', unit='',
                             array=source_cat[:,SXcolumn['flags']],
                             disp='SExtractor flags'),
               pyfits.Column(name='OTA',            format='I', unit='',
                             array=source_cat[:,SXcolumn['ota']],
                             disp='source OTA'),

               pyfits.Column(name='MAG_AUTO',        format='E', unit='mag',
                             array=source_cat[:,SXcolumn['mag_auto']],
                             disp='auto-mag'),
               pyfits.Column(name='MAGERR_AUTO',     format='E', unit='mag',
                             array=source_cat[:,SXcolumn['mag_err_auto']],
                             disp=''),
        
               pyfits.Column(name='FLUX_MAX',       format='E', unit='counts',
                             array=source_cat[:,SXcolumn['flux_max']],
                             disp='max count rate'),
               pyfits.Column(name='AWIN_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['major_axis']],
                             disp='major semi-axis'),
               pyfits.Column(name='BWIN_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['minor_axis']],
                             disp='minor semi-axis'),
               pyfits.Column(name='THETAWIN_IMAGE', format='E', unit='degrees',
                             array=source_cat[:,SXcolumn['position_angle']],
                             disp='position angle'),
               pyfits.Column(name='ELONGATION',     format='E', unit='',
                             array=source_cat[:,SXcolumn['elongation']],
                             disp='elongation'),
               pyfits.Column(name='ELLIPTICITY',    format='E', unit='',
                             array=source_cat[:,SXcolumn['ellipticity']],
                             disp='ellipticity'),
           ]

    # Add all aperture photometry
    for i in range(len(SXapertures)):
        apsize = SXapertures[i]
        col_mag = "mag_aper_%0.1f" % (apsize)
        col_magerr = "mag_err_%0.1f" % (apsize)
        columns.append(pyfits.Column(name='MAG_D%02d' % (int(apsize*10.)), 
                                     format='E', unit='mag',
                                     array=source_cat[:,SXcolumn[col_mag]],
                                     disp='aperture mag %0.1f arcsec diameter' % (apsize)
                                 )
        )
        columns.append(pyfits.Column(name='MAGERR_D%02d' % (int(apsize*10.)),
                                     format='E', unit='mag',
                                     array=source_cat[:,SXcolumn[col_magerr]],
                                     disp='ap. mag error %0.1f arcsec diameter' % (apsize)
                                 )
        )

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("CAT.ODI", comment=None)
    return tbhdu












def set_default_options(options_in=None):
    """
    Initialize the options dictionary by defining all available options and
    assigning safe and reasonable default values.

    Parameters
    ----------
    options_in : dictionary

        If a directory already exists, add to it, otherwise create a new one.

    Returns
    -------
    options : dictionary

        The options dictionary with default values set.

    """
    
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

    options['gain_correct'] = False
    options['gain_method'] = None

    options['nonlinearity'] = None
    options['nonlinearity-set'] = False

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
    options['otalevelplots'] = True

    options['additional_fits_headers'] = {}

    options['structure_qa_subdirs'] = False
    options['structure_qa_subdir_name'] = "QA"
    options['create_qaplots'] = True

    options['skip_guide_ota'] = False

    options['log_setup'] = None

    options['skip_otas'] = []

    options['nonsidereal'] = None

    options['fitradialZP'] = False

    options['techdata'] = None

    options['crj'] = 0
    options['crj_sigclip'] = 5.0
    options['crj_sigfrac'] = 0.3
    options['crj_objlim'] = 5.0
    options['crj_saturation'] = -1
    options['crj_method'] = "cy"

    return options












def check_filename_directory(given, default_filename):
    """
    Some of the options support either a directory or a filename. This function
    checks if the input is a directory or a filename. In the first case, add 
    the specified default filename and return the filename. In the latter case, 
    simply return the filename.

    Parameters
    ----------
    given : string

        The value specified by the user, either a directory or filename

    default_filename : string

        In the case the user specified only a directory, append the name of this
        filename

    Returns
    -------
    filename

    Example
    -------
    This function is called e.g. during bias-subtraction. Using the -bias 
    option, the user can specify the bias file to be used during the reduction.
    However, the default way of specifying the bias file is to simply give the 
    directory with the calibration files, and collectcells adds the 'bias.fits'
    during execution.

    """

    if (given == None):
        return None
    elif (os.path.isfile(given)):
        return given
    elif (os.path.isdir(given)):
        return "%s/%s" % (given, default_filename)

    return ""











def read_options_from_commandline(options=None):
    """
    Read all command line options and store them in the options dictionary.

    """

    logger = logging.getLogger("ReadOptions")

    if (options == None):
        options = set_default_options()

    options['verbose'] = cmdline_arg_isset("-verbose")

    # Handle all reduction flags from command line
    if (cmdline_arg_isset("-cals")):
        cals_dir = get_cmdline_arg("-cals")
        if (not os.path.isdir(cals_dir)):
            logger.critical("The specified cals-directory (%s) does not exist!!!" % (cals_dir))
            stdout_write("\n\n   The specified cals-directory (%s) does not exist!!!\n\n\n" % (cals_dir))
            sys.exit(0)

        options['bias_dir'] = cals_dir
        options['dark_dir'] = cals_dir
        options['flat_dir'] = cals_dir
        options['techdata'] = cals_dir
        logger.info("tech-data: %s" % (options['techdata']))

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

    
    options['gain_correct'] = cmdline_arg_isset("-gain")
    options['gain_method'] = cmdline_arg_set_or_default("-gain", None)

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
    options['wcs_distortion'] = root_dir + "/2mass_distort5.fits"
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

    options['nonlinearity-set'] = cmdline_arg_isset("-nonlinearity")
    options['nonlinearity'] = cmdline_arg_set_or_default("-nonlinearity", None)

    if (cmdline_arg_isset('-plotformat')):
        inputstr = cmdline_arg_set_or_default("-plotformat", "png")
        options['plotformat'] = inputstr.split(",")
        print "writing plots as ",options['plotformat']
        
    options['otalevelplots'] = not cmdline_arg_isset("-nootalevelplots")

    options['structure_qa_subdirs'] = cmdline_arg_isset("-qasubdirs")
    if (cmdline_arg_isset('-qasubdirname')):
        options['structure_qa_subdir_name'] = cmdline_arg_set_or_default('-qasubdirname', "QA")
        options['structure_qa_subdirs'] = True

    options['create_qaplots'] = not cmdline_arg_isset("-noqaplots")
    
    # Now loop over all headers again and isolate the -addfitskey entries
    for entry in sys.argv[1:]:
        if (entry[:11] == "-addfitskey"):
            key_value = entry[12:]
            key, value = key_value.split(',', 2)
            print "Adding fits keyword %s = %s" % (key, value) 

            options['additional_fits_headers'][key] = value

    # Determine which, if any, OTAs are to be skipped
    if (cmdline_arg_isset("-skipota")):
        toskip = cmdline_arg_set_or_default("-skipota", "")
        items = toskip.split(',')
        for i in items:
            options['skip_otas'].append(int(i))

    if (cmdline_arg_isset("-nonsidereal")):
        logger.debug("Found -nonsidereal command line flag")
        value = cmdline_arg_set_or_default("-nonsidereal", None)
        if (not value == None):
            items = value.split(',')
            if (len(items) == 3):
                ns = {}
                ns['dra'] = float(items[0])
                ns['ddec'] = float(items[1])
                ns['ref'] = items[2]
                ns['ref_mjd'] = None
                try:
                    ns['ref_mjd'] = float(ns['ref'])
                    logger.debug("Found reference MJD (%f) on command line" % (ns['ref_mjd']))
                except:
                    if (os.path.isfile(ns['ref'])):
                        hdulist = pyfits.open(ns['ref'])
                        if ("MJD-OBS" in hdulist[0].header):
                            ns['ref_mjd'] = hdulist[0].header['MJD-OBS']
                            logger.debug("Found reference MJD (%f) in file %s" % (ns['ref_mjd'], ns['ref']))
                        hdulist.close()
                    else:
                        logger.critical("Could not find the reference file for the -nonsidereal option!")
                        logger.critical("Disabling nonsidereal WCS correction")
                    pass
                if (not ns['ref_mjd'] == None):
                    options['nonsidereal'] = ns
            else:
                logger.critical("I don't understand the -nonsidereal parameter")
        logger.debug("non-sidereal setup: %s" % (str(ns)))

    options['fitradialZP'] = cmdline_arg_isset("-fitradialZP")

    options['techdata'] = cmdline_arg_set_or_default("-techdata", options['techdata'])

    options['crj'] = int(cmdline_arg_set_or_default("-crj", 0))
    options['crj_sigclip'] = float(cmdline_arg_set_or_default("-crjsigclip", sitesetup.crj_sigclip))
    options['crj_sigfrac'] = float(cmdline_arg_set_or_default("-crjsigfrac", sitesetup.crj_sigfrac))
    options['crj_objlim'] = float(cmdline_arg_set_or_default("-crjobjlim", sitesetup.crj_objlim))
    options['crj_saturation'] = float(cmdline_arg_set_or_default("-crjsaturation", sitesetup.crj_saturation))
    options['crj_method'] = cmdline_arg_set_or_default("-crjmethod", "cy")

    return options











if __name__ == "__main__":

    if (len(sys.argv) <= 1 or sys.argv[1] == "-help"):
        #print help('podi_matchpupilghost')
        import podi_collectcells as me
        print me.__doc__
        sys.exit(0)

    m = multiprocessing.Manager()
    process_tracker = m.Queue()

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

    # Setup everything we need for logging
    podi_logging.setup_logging(options)

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
            if (cmdline_arg_isset("-timeout")):
                timeout = float(cmdline_arg_set_or_default("-timeout", 900))
                print "Setting timeout to",timeout,"seconds"
                collectcells_with_timeout(input, outputfile, options=options,
                                          timeout=timeout,
                                          process_tracker=process_tracker)
            else:
                collectcells(input, outputfile, process_tracker=process_tracker, options=options)
    except:
        print "Cleaning up left over child processes"
        podi_logging.log_exception()
        kill_all_child_processes(process_tracker)


    finally:
        podi_logging.shutdown_logging(options)
