#!/usr/bin/env python3
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

from __future__ import print_function

import sys
import os
import signal
import astropy.io.fits as pyfits
import numpy
import scipy
import scipy.optimize
import ephem
import traceback
#import psutil
import datetime
import warnings

import queue
import threading
import multiprocessing
import multiprocessing.reduction
import ctypes
import time
import logging
import itertools
import bottleneck
import errno
import psutil

from astropy.io.fits.verify import VerifyWarning
warnings.simplefilter('ignore', category=VerifyWarning)

from podi_plotting import *

#
# Un-comment this block to get warnings with tracebacks for easier locating
#
# import warnings
# #warnings.simplefilter("error")
# warnings.simplefilter("ignore", RuntimeWarning)
# import traceback
# import warnings
# import sys
# def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
#     print"\n"*3
#     traceback.print_stack()
#     log = file if hasattr(file,'write') else sys.stderr
#     log.write(warnings.formatwarning(message, category, filename, lineno, line))
#     print "\n"*3
# warnings.showwarning = warn_with_traceback

gain_correct_frames = False
from podi_definitions import *
from podi_commandline import *
from wiyn_filters import *
#import podi_findstars
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
import podi_illumcorr
import podi_almanac
import podi_focalplanelayout
import podi_guidestars
import podi_shifthistory
import podi_photflat

import psf_quality
import podi_readfitscat

import podi_associations
import podi_calibrations
from podi_calibrations import check_filename_directory

from podi_reductionlog import *
from version import record_pipeline_versioning

from astLib import astWCS

from sharedmemory import SharedMemory

fix_cpu_count = False

if (sitesetup.number_cpus == "auto"):
    try:
        number_cpus = multiprocessing.cpu_count()
        # print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
        if (number_cpus > sitesetup.max_cpu_count and sitesetup.max_cpu_count > 1):
            number_cpus = sitesetup.max_cpu_count
            # print "... but using only %d of them!" % (number_cpus)
    except:
        pass
else:
    number_cpus = sitesetup.number_cpus




def read_techdata(techdata_hdulist, ota_x, ota_y, wm_cellx, wm_celly):

    logger = logging.getLogger("ReadTechData")
    #logger.info("Reading from %s" % (str(techdata_hdulist)))

    gain, readnoise, readnoise_e = None, None, None
    cellidx_x = ota_x*8+wm_cellx
    cellidx_y = ota_y*8+(7-wm_celly)


    try:
        gain = techdata_hdulist['GAIN'].data[cellidx_y, cellidx_x]
        readnoise = techdata_hdulist['READNOISE'].data[cellidx_y, cellidx_x]
        readnoise_e = techdata_hdulist['READNOISE_E'].data[cellidx_y, cellidx_x]
    except:
        logger.error("TECHDATA file is likely corrupted!")
        pass

    return gain, readnoise, readnoise_e


def apply_wcs_distortion(filename, hdu, binning, reduction_log=None):

    logger = logging.getLogger("ApplyWCSmodel")

    if (reduction_log is not None):
        reduction_log.attempt('wcs_dist')

    try:
        wcs = pyfits.open(filename)
    except:
        if (reduction_log is not None):
            reduction_log.fail('wcs_dist')
        logger.error("Could not open WCS distortion model (%s)" % (filename))
        return False

    extname = hdu.header['EXTNAME']

    try:
        wcs_header = wcs[extname].header
    except:
        if (reduction_log is not None):
            reduction_log.fail('wcs_dist')
        logger.warning("Could not find distortion model for %s" % (extname))
        return False

    try:
        for hdr_name in ('CRPIX1', 'CRPIX2'):
            wcs_header[hdr_name] /= binning
        for hdr_name in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
            wcs_header[hdr_name] *= binning

        cards = wcs_header.cards
        for (keyword, value, comment) in cards:
            if (keyword not in ['CRVAL1', 'CRVAL2']):
                hdu.header[keyword] = (value, comment)

        d_crval1 = wcs_header['CRVAL1']
        if (d_crval1 > 180):
            d_crval1 -= 360
        hdu.header['CRVAL2'] += wcs_header['CRVAL2']
        hdu.header["CRVAL1"] += d_crval1 / math.cos(math.radians(hdu.header['CRVAL2']))
        # print "change crval1 by",wcs_header['CRVAL1'], d_crval1, wcs_header['CRVAL1'] / math.cos(math.radians(hdu.header['CRVAL2']))

        # Make sure to write RAs that are positive
        if (hdu.header["CRVAL1"] < 0):
            hdu.header['CRVAL1'] -= math.floor(hdu.header['CRVAL1'] / 360.) * 360.
        elif (hdu.header['CRVAL1'] > 360.):
            hdu.header['CRVAL1'] = math.fmod(hdu.header['CRVAL1'], 360.0)

        if ('RADESYS' in hdu.header):
            hdu.header['RADESYS'] = 'ICRS'

        # Rewrite TAN to TPV to make IRAF happy :-(
        for hdr_name in ['CTYPE1', 'CTYPE2']:
            if (hdr_name in hdu.header):
                hdu.header[hdr_name] = hdu.header[hdr_name].replace("TAN", "TPV")

    except:
        logger.critical("something went wrong while applying the WCS model")
        if (reduction_log is not None):
            reduction_log.partial_fail('wcs_dist')
        podi_logging.log_exception()

    if (reduction_log is not None):
        reduction_log.success('wcs_dist')
    return True


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

    reduction_log = ReductionLog()
    data_products = {
        "hdu": None,
        "header": None,
        "wcsdata": None,
        "sky-samples": None,
        "sky": None,
        "tech-header": None,
        "fringe_scaling": None,
        "fringe-template": None,
        "source-cat": None,
        'pupilghost-scaling': None,
        'pupilghost-template': None,
        'reduction_files_used': None,
        'reduction_log': reduction_log,
        'psf': None,
        }
   
    if (not os.path.isfile(filename)):
        stdout_write("Couldn't find file %s ..." % (filename))
    else:
        # Create an fits extension to hold the output
        hdu = pyfits.ImageHDU()
        # log_svn_version(hdu.header)

        # Set the inherit keyword so that the headers removed from each 
        # extension are instead inherited from the primary
        hdu.header["INHERIT"] = (True, "Inherit headers from PrimaryHDU")

        # also create a tech-header to keep track of values used in this frame
        tech_header = pyfits.Header()

        # Keep track of what input files we used
        reduction_files_used = {}

        try:
            hdulist = pyfits.open(filename, memmap=False)
        except:
            # This happed with corrupt files
            return data_products

        # Check if we can find all 1+64 extensions
        if (len(hdulist) < 65):
            # Something is wrong with this file
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
        hdu.name = extname
        hdu.header['OTA'] = (ota, "OTA designation")

        # podi_logging.podi_getlogger("%s - %s" % (obsid, extname), options['log_setup'])
        logger = logging.getLogger("%s - %s" % (obsid, extname))

        fpl = podi_focalplanelayout.FocalPlaneLayout(hdulist[0])
        logger.debug("Focalplane Layout: %s" % (fpl.layout))

        # Now copy the headers from the original file into the new one
        firstkey = None
        cards = hdulist[0].header.cards
        for (keyword, value, comment) in cards:
            firstkey = keyword if firstkey is None else firstkey
            hdu.header[keyword] = (value, comment)
        add_fits_header_title(hdu.header, "Instrument telemetry from raw data", firstkey)

        #
        # Add default headers here
        #
        hdu.header['SOFTBIN'] = 0
        if ('RADESYS' in hdu.header):
            del hdu.header['RADESYS']
        add_fits_header_title(hdu.header, "Pipeline added/modified metadata", 'SOFTBIN')

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

        mastercals = podi_calibrations.ODICalibrations(
            cmdline_options=options,
            hdulist=hdulist)

        # Now go through each of the 8 lines
        if (not mastercals.apply_crosstalk()):
            reduction_log.not_selected('crosstalk')
        elif (mastercals.crosstalk(mjd, ota) is None):
            logger.warning("Cross-talk correction requested, but failed!")
            reduction_log.fail('crosstalk')
        else:
            xtalk_file = mastercals.crosstalk(mjd, ota)
            logger.debug("Starting crosstalk correction (%s)" % (extname))
            reduction_files_used['crosstalk'] = xtalk_file

            outcome = podi_crosstalk.apply_crosstalk_correction(
                hdulist,
                xtalk_file = xtalk_file,
                fpl = fpl,
                extname2id = extname2id,
                options = options,
                reduction_log = reduction_log
            )
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

        # nonlin_data = None
        # if (options['nonlinearity-set'] or options['gain_method'] == "relative"):
        #     nonlinearity_file = options['nonlinearity']
        #     if (options['nonlinearity'] == None or
        #         options['nonlinearity'] == "" or
        #         not os.path.isfile(nonlinearity_file)):
        #         nonlinearity_file = podi_calibrations.find_nonlinearity_coefficient_file(mjd, options)
        #     if (options['verbose']):
        #         print "Using non-linearity coefficients from",nonlinearity_file
        #     logger.debug("Using non-linearity coefficients from file %s"  % (nonlinearity_file))
        #     nonlin_data = podi_nonlinearity.load_nonlinearity_correction_table(nonlinearity_file, ota)
        #     if (options['nonlinearity-set']):
        #         reduction_files_used['nonlinearity'] = nonlinearity_file

        nonlin_data = None
        if (not mastercals.apply_nonlinearity()):
            reduction_log.not_selected('nonlinearity')
        elif (mastercals.nonlinearity(mjd) is None):
            reduction_log.fail('nonlinearity')
        else:
            reduction_log.attempt('nonlinearity')
            nonlinearity_file = mastercals.nonlinearity(mjd=mjd)
            logger.debug("Using non-linearity coefficients from file %s" % (nonlinearity_file))
            nonlin_data = podi_nonlinearity.load_nonlinearity_correction_table(nonlinearity_file, ota)
            if (nonlin_data is not None):
                reduction_files_used['nonlinearity'] = nonlinearity_file
                reduction_log.success('nonlinearity')
            else:
                logger.warning("Unable to load non-linearity corrections (%s)" % (nonlinearity_file))
                reduction_log.fail('nonlinearity')

        relative_gains = None
        if (mastercals.apply_relative_gain()):
            # In this mode, we also require a non-linearity correction file to be loaded
            if (nonlin_data is not None):
                # all done, we already loaded the file as part of the non-linearity correction startup
                relative_gains = nonlin_data
                pass
            else:
                relative_gains_file = mastercals.nonlinearity(mjd=mjd)
                if (relative_gains_file is not None):
                    relative_gains = podi_nonlinearity.load_nonlinearity_correction_table(relative_gains_file, ota)

        #
        # Search, first in the flat-field, then in the bias-frame for a 
        # TECHDATA extension that holds GAINS and readout noises for each cell. 
        #
        techdata = None
        techhdulist = None
        if (options['techdata'] is not None):
            techfile = None

            # Set the filename for the TECHDATA extension based on the user-option
            if (options['techdata'] == "from_flat" and options['flat_dir'] is not None):

                for ft in sitesetup.flat_order:
                    techfile = check_filename_directory(options['flat_dir'], 
                        "%s_%s_bin%d.fits" % (ft, filter_name, binning))
                    if (os.path.isfile(techfile)):
                        break
                
            elif (options['techdata'] == "from_bias" and options['bias_dir'] is not None):
                techfile = check_filename_directory(options['bias_dir'], "bias_bin%s.fits" % (binning))

            else:
                for ft in sitesetup.flat_order:
                    techfile = check_filename_directory(options['techdata'], 
                        "techdata_%s_%s_bin%d.fits" % (ft, filter_name, binning))
                    if (os.path.isfile(techfile)):
                        break

            # Check if the specified file exists and read the data if possible
            if (os.path.isfile(techfile)):
                logger.debug("Reading techdata from file %s" % (techfile))
                techhdulist = pyfits.open(techfile)
                reduction_files_used['techdata'] = techfile
            else:
                techfile = "%s/techdata.fits" % (sitesetup.exec_dir)
                if (os.path.isfile(techfile)):
                    logger.debug("Reading techdata from file %s" % (techfile))
                    techhdulist = pyfits.open(techfile)
                    reduction_files_used['techdata'] = techfile
                else:
                    logger.debug("Was looking for techfile %s but couldn't find it" % (techfile))

        all_gains = numpy.ones(shape=(8,8)) * -99
        all_readnoise = numpy.ones(shape=(8,8)) * -99
        all_readnoise_electrons = numpy.ones(shape=(8,8)) * -99

        #
        # Handle the trimcell option
        #
        hdu.header['TRIMCELL'] = (-1, "trim cell edges (-1: no)")
        if (options['trimcell'] is None):
            reduction_log.not_selected('trimcell')
        else:
            reduction_log.success('trimcell')
            hdu.header['TRIMCELL'] = options['trimcell']

        reduction_log.attempt('overscan')
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
            cellmode_id = get_cellmode(hdulist[0].header, hdulist[cell].header, fpl)
            if (not cellmode_id == 0):
                # This means it either broken (id=-1) or in video-mode (id=1)
                if (options['keep_cells'] == False):
                    # no keep_cells option
                    continue
                elif (options['keep_cells'] is None):
                    # keep all cells
                    pass
                else:
                    cell_id = "%02d.%1d%1d" % (ota, wm_cellx, wm_celly)
                    if (cell_id in options['keep_cells']):
                        pass
                    else:
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
            datasec = extract_datasec_from_cell(
                hdulist[cell].data, binning,
                trimcell=options['trimcell'])
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
            bx, tx, by, ty = cell2ota__get_target_region(
                wm_cellx, wm_celly, binning,
                trimcell=options['trimcell'])
            merged[by:ty,bx:tx] = datasec

            #
            # Now that the cell is pre-reduced (merged and overscan-subtracted), assemble 
            # the tech-header by copying the information from the input techdata 
            # to the output techdata
            #
            if (techhdulist is not None):
                gain, readnoise, readnoise_e = read_techdata(
                    techhdulist, ota_c_x, ota_c_y, wm_cellx, wm_celly)
                if (gain is not None and readnoise is not None):
                    logger.debug("Using tech-data for gain & readnoise: %f %f" % (gain, readnoise))
                else:
                    logger.debug("No gain/readnoise data available")

                # Also set the gain and readnoise values to be used later for the gain correction
                all_gains[wm_cellx, wm_celly] = \
                    gain if gain is not None else \
                    hdulist[cell].header['GAIN'] if 'GAIN' in hdulist[cell].header else backup_gain
                all_readnoise[wm_cellx, wm_celly] = \
                    readnoise if gain is not None else backup_readnoise
                all_readnoise_electrons[wm_cellx, wm_celly] = \
                    readnoise_e if readnoise_e is not None else backup_readnoise_electrons
            else:
                all_gains[wm_cellx, wm_celly] = \
                    hdulist[cell].header['GAIN'] if 'GAIN' in hdulist[cell].header else backup_gain
                all_readnoise[wm_cellx, wm_celly] = backup_readnoise
                all_readnoise_electrons[wm_cellx, wm_celly] = backup_readnoise_electrons
                
            # work on next cell
        reduction_log.success('overscan')
        logger.debug("Collected all cells for OTA %02d of %s" % (ota, obsid))

        # 
        # At this point we have a 4x4 Kpixel array with all cells merged
        #

        # If we are to do some bias subtraction:
        if (not mastercals.apply_bias()):
            reduction_log.not_selected('bias')
        elif (mastercals.bias() is None):
            reduction_log.fail('bias')
        else:
            bias_filename = mastercals.bias()
            bias = pyfits.open(bias_filename)
            reduction_files_used['bias'] = bias_filename

            # Search for the bias data for the current OTA
            for bias_ext in bias[1:]:
                if (not is_image_extension(bias_ext)):
                    continue
                fppos_bias = bias_ext.header['FPPOS']
                if (fppos_bias == fppos):
                    # This is the one
                    try:
                        if (not options['keep_cells'] == False):
                            bias_ext.data[numpy.isnan(bias_ext.data)] = 0.
                        merged -= bias_ext.data
                        logger.debug("Subtracting bias: %s" % (bias_filename))
                    except:
                        logger.warning("Unable to subtract bias, dimensions don't match (data: %s, bias: %s)" % (
                            str(merged.shape), str(bias_ext.data.shape)))
                        reduction_log.fail('bias')
                        pass
                    break

            bias.close()
            hdu.header.add_history("CC-BIAS: %s" % (os.path.abspath(bias_filename)))
            del bias
            reduction_log.success('bias')
 
        #
        # Do some dark subtraction:
        # Add at some point: use different darks for all detectors switched on 
        # to minimize the integration glow in guide-OTAs
        #
        if (not mastercals.apply_dark()):
            reduction_log.not_selected('dark')
        elif (mastercals.dark() is None):
            reduction_log.fail('dark')
        else:
            # For now assume all detectors are switched on
            dark_filename = mastercals.dark()
            dark = pyfits.open(dark_filename)
            reduction_files_used['dark'] = dark_filename
            if ('EXPMEAS' in dark[0].header):
                darktime = dark[0].header['EXPMEAS']
            elif ('EXPTIME' in dark[0].header):
                darktime = dark[0].header['EXPTIME']
            else:
                darktime = 1.

            # Search for the flatfield data for the current OTA
            for dark_ext in dark[1:]:
                if (not is_image_extension(dark_ext)):
                    continue
                fppos_dark = dark_ext.header['FPPOS']
                if (fppos_dark == fppos):
                    # This is the one
                    dark_scaling = exposure_time / darktime
                    try:
                        if (not options['keep_cells'] == False):
                            dark_ext.data[numpy.isnan(dark_ext.data)] = 0.
                        merged -= (dark_ext.data * dark_scaling)
                        logger.debug("Subtracting dark: %s (scaling=%.2f)" % (dark_filename, dark_scaling))
                    except:
                        logger.warning("Unable to subtract dark, dimensions don't match (data: %s, dark: %s)" % (
                            str(merged.shape), str(dark_ext.data.shape)))
                        reduction_log.fail('dark')
                        pass
                    break

            dark.close()
            hdu.header.add_history("CC-DARK: %s" % (os.path.abspath(dark_filename)))
            del dark
            reduction_log.success('dark')

        # By default, mark the frame as not affected by the pupilghost. This
        # might be overwritten if the flat-field has PG specifications.
        hdu.header['PGAFCTD'] = False

        # If the third parameter points to a directory with flat-fields
        hdu.header['GAIN'] = 1.3
        gain_from_flatfield = None
        pupilghost_center_x = pupilghost_center_y = None
        if (not mastercals.apply_flat()):
            reduction_log.not_selected('flat')
        elif (mastercals.flat(sitesetup.flat_order) is None):
            reduction_log.fail('flat')
        else:
            flatfield_filename = mastercals.flat(sitesetup.flat_order)
            flatfield = pyfits.open(flatfield_filename)
            reduction_files_used['flat'] = flatfield_filename

            # Search for the flatfield data for the current OTA
            for ff_ext in flatfield[1:]:
                if (not is_image_extension(ff_ext)):
                    continue
                fppos_flatfield = ff_ext.header['FPPOS']
                if (fppos_flatfield == fppos):
                    # This is the one
                    try:
                        if (not options['keep_cells'] == False):
                            ff_ext.data[numpy.isnan(ff_ext.data)] = 1.
                        merged /= ff_ext.data
                        logger.debug("Dividing by flatfield: %s" % (flatfield_filename))
                    except:
                        logger.warning("Unable to apply flat-field, dimensions don't match (data: %s, flat: %s)" % (
                            str(merged.shape), str(ff_ext.data.shape)))
                        reduction_log.partial_fail('flat')
                        pass

                    # If normalizing with the flat-field, overwrite the gain
                    # keyword with the average gain value of the flatfield.
                    ff_gain = flatfield[0].header['GAIN'] \
            if ('GAIN' in flatfield[0].header and flatfield[0].header['GAIN'] > 0) else -1.
                    gain_from_flatfield = ff_gain

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

                        for pgheader in (
                                'PGCENTER', 'PGCNTR_X', 'PGCNTR_Y',
                                'PGTMPL_X', 'PGTMPL_Y',
                                'PGREG_X1', 'PGREG_X2', 'PGREG_Y1', 'PGREG_Y2',
                                'PGEFCTVX', 'PGEFCTVY',
                                'PGSCALNG',
                                'PGROTANG',
                        ):
                            if (pgheader in ff_ext.header):
                                hdu.header[pgheader] = (ff_ext.header[pgheader], ff_ext.header.comments[pgheader])

                        pupilghost_center_x = ff_ext.header['PGCNTR_X'] if 'PGCNTR_X' in ff_ext.header else numpy.NaN
                        pupilghost_center_y = ff_ext.header['PGCNTR_Y'] if 'PGCNTR_Y' in ff_ext.header else numpy.NaN

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
            reduction_log.success('flat')

        #
        # Apply illumination correction, if requested
        #
        if (options['illumcorr_dir'] is None):
            reduction_log.not_selected('illumcorr')
        else:
            illumcorr_filename = podi_illumcorr.get_illumination_filename(
                options['illumcorr_dir'], filter_name, binning)

            if (not os.path.isfile(illumcorr_filename)):
                reduction_log.fail('illumcorr')
            else:
                illumcorr = pyfits.open(illumcorr_filename)
                reduction_files_used['illumination'] = illumcorr_filename

                # Search for the flatfield data for the current OTA
                for ff_ext in illumcorr:
                    if (not is_image_extension(ff_ext)):
                        continue
                    fppos_flatfield = ff_ext.header['FPPOS']
                    if (fppos_flatfield == fppos):
                        # This is the one
                        try:
                            merged /= ff_ext.data
                            logger.debug("Applying illumination correction: %s" % (illumcorr_filename))
                        except:
                            logger.warning(
                                "Unable to apply illumination correction, dimensions don't match (data: %s, i.c.: %s)" % (
                                str(merged.shape), str(ff_ext.data.shape)))
                            reduction_log.fail('illumcorr')
                            pass
                        break
                        
                illumcorr.close()
                del illumcorr
                reduction_log.success('illumcorr')

        # Finally, apply bad pixel masks 
        # Determine which region file we need
        # This function only returns valid filenames, or None otherwise
        # bpm_region_file = fpl.get_badpixel_regionfile(options['bpm_dir'], fppos)
        if (not mastercals.apply_bpm()):
            reduction_log.not_selected('badpixels')
        elif (mastercals.bpm(ota=fppos) is None):
            reduction_log.fail('badpixels')
        else:
            bpm_region_file = mastercals.bpm(ota=fppos)
        # if (bpm_region_file is None):
        #     reduction_log.not_selected('badpixels')
        # else:
            # Apply the bad pixel regions to file, marking all bad pixels as NaNs
            logger.debug("Applying BPM file: %s" % (bpm_region_file))
            mask_broken_regions(merged, bpm_region_file, reduction_log=reduction_log)
            reduction_files_used['bpm'] = bpm_region_file
            reduction_log.success('badpixels')

        #
        # Now apply the gain correction
        #

        logger.debug("GAIN setting:"+str(options['gain_correct']))
        hdu.header['GAINMTHD'] = ('none', "gain method")
        if (not options['gain_correct']):
            reduction_log.not_selected('gain')
        else:
            # Correct for the gain variations in each cell
            logger.debug("Applying gain correction (OTA %02d) - method: %s" % (ota, 
                options['gain_method'] if options['gain_method'] is not None else "default:techdata"))

            if (mastercals.apply_relative_gain()):
                relative_gains_file = mastercals.nonlinearity(mjd=mjd)
                if (relative_gains_file is None):
                    reduction_log.fail('gain')
                else:
                    reduction_files_used['gain'] = relative_gains_file
                
                    # Find the relative gain correction factor based on the non-linearity correction data
                    logger.debug("Apply gain correction from nonlinearity data")

                    for cx, cy in itertools.product(range(8), repeat=2):
                        x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy, binning)
                        merged[y1:y2, x1:x2], gain = podi_nonlinearity.apply_gain_correction(
                            merged[y1:y2, x1:x2], cx, cy, nonlin_data, return_gain=True)
                        if (gain > 0):
                            all_gains[cx,cy] /= gain
                    reduction_log.success('gain')
                    hdu.header['GAINMTHD'] = "relative"
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
                reduction_log.success('gain')
                hdu.header['GAINMTHD'] = "raw_header"

            else:
                if (techdata is not None):
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
                    reduction_log.success('gain')
                    hdu.header['GAINMTHD'] = 'techdata'
                else:
                    reduction_log.failed('gain')
                    logger.warning("GAIN correction using TECHDATA requested, but can't find a tech-data file")


        #
        # Compute the average gain, read-noise, etc and store how many cells contributed
        #
        # Set the GAIN to be used for photometry. This gain is 1.0 if we 
        # already applied the gain correction to the data
        valid_gain = (all_gains > 0) & (numpy.isfinite(all_gains))
        gain_avg = numpy.mean(all_gains[valid_gain]) if numpy.sum(valid_gain) > 0 else 1.3
        hdu.header['GAINX'] = (gain_avg, 'gain averaged over all cells')
        hdu.header['NGAIN'] = (numpy.sum(valid_gain), 'number of cells contrib. to avg. gain')

        valid_ron = (all_readnoise > 0) & (numpy.isfinite(all_readnoise))
        ron_avg = numpy.mean(all_readnoise[valid_ron]) if numpy.sum(valid_ron) > 0 else 6.0
        hdu.header['RDNOISE'] = (ron_avg, 'readnoise averaged over all cells')
        hdu.header['NRDNOISE'] = (numpy.sum(valid_ron), 'number of cells contrib. to avg. readnoise')

        valid_rone = (all_readnoise_electrons > 0) & (numpy.isfinite(all_readnoise_electrons))
        rone_avg = numpy.mean(all_readnoise_electrons[valid_rone]) if numpy.sum(valid_rone) > 0 else 8.5
        hdu.header['RDNOISEE'] = (rone_avg, 'readnoise in e- averaged over all cells')
        hdu.header['NRDNOSEE'] = (numpy.sum(valid_rone), 'number cells contrib. to avg. readnoise in e-')

        if (gain_from_flatfield is not None):
            hdu.header['GAIN'] = gain_from_flatfield
            logger.debug("Overwriting GAIN keyword: %f" % (ff_gain))
        else:
            hdu.header['GAIN'] = hdu.header['GAINX']

        # # Compute the gain across this OTA
        # print extname,"-->\n",all_gains
        # valid_gain = (all_gains > 0) & (numpy.isfinite(all_gains))
        # if (numpy.sum(valid_gain) <= 0):
        #     gain_avg = 1.3
        # else:
        #     gain_avg = numpy.mean(all_gains[valid_gain])
        # hdu.header['GAIN_RAW'] = (gain_avg, 'gain averaged over all cells')
        # hdu.header['NGAINRAW'] = (numpy.sum(valid_gain), 'number of cells contrib. to avg. gain')

        
        
                
        #
        # Optionally, normalize the frame by some header value
        #
        if (not "normalize" in options):
            reduction_log.not_selected('exptime_norm')
        else:
            norm_factor = 1
            if (options['normalize'] == "EXPTIME" or
                options['normalize'] == "EXPMEAS"):
                # get the factor from the header
                if (options['normalize'] in hdu.header):
                    norm_factor = hdu.header[options['normalize']]
                elif 'EXPTIME' in hdu.header:
                    norm_factor = hdu.header['EXPTIME']
                else:
                    norm_factor = 1.0

                if (norm_factor > 0):
                    logger.debug("normalizing data with constant %f (%s)" % (
                            norm_factor, options['normalize']))
                    # normalize the data
                    merged /= norm_factor
                    # and fix the EXPTIME/EXPMEAS header as well
                    if ('EXPTIME' in hdu.header): hdu.header['EXPTIME'] /= norm_factor
                    if ('EXPMEAS' in hdu.header): hdu.header['EXPMEAS'] /= norm_factor
            hdu.header['NORMALIZ'] = (norm_factor, "normalization constant")
            reduction_log.success('exptime_norm')
                
        #
        # If persistency correction is requested, perform it now
        #
        if (options['persistency_dir'] is None):
            reduction_log.not_selected('persistency')
            reduction_log.not_selected('saturation')
        else:
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

            if (saturated_thisframe is not None):
                logger.debug("Using the following saturation frames: %s" % (str(saturated_thisframe)))
                merged = podi_persistency.mask_saturation_defects(saturated_thisframe, ota, merged)
            hdu.header["SATMASK"] = saturated_thisframe
            reduction_log.success('saturation')

            # Also pick all files within a given MJD range, and apply the 
            # appropriate correction (which, for now, is simply masking all pixels)
            filelist = podi_persistency.select_from_saturation_tables(full_filelist, mjd, 
                                                                      [1,options['max_persistency_time']])
            if (len(filelist) <= 0):
                reduction_log.no_data('persistency')
            else:
                # Extract only the filenames from the filelist dictionary
                persistency_files_used = []
                for assoc_persistency_file, assoc_mjd in filelist.items():
                    persistency_files_used.append(assoc_persistency_file)
                logger.debug("Using persistency catalogs in %s" % (str(persistency_files_used)))
                reduction_files_used['persistency'] = persistency_files_used
                
                merged = podi_persistency.correct_persistency_effects(ota, merged, mjd, filelist)
                persistency_catalog_counter = 0
                for filename in filelist:
                    persistency_catalog_counter += 1
                    keyname = "PERSIS%02d" % (persistency_catalog_counter)
                    hdu.header[keyname] = filename
                reduction_log.success('persistency')


        #
        # If requested, determine the optimal fringe scaling
        #
        fringe_scaling = None
        fringe_scaling_median, fringe_scaling_std = -1, -1
        if (not mastercals.apply_fringe()):
            reduction_log.not_selected('fringe')
        elif (mastercals.fringe(mjd=mjd) is None or
              mastercals.fringevector(ota=ota) is None):
            reduction_log.fail('fringe')
            logger.warning("De-fringing selected, but could not find either fringe template or fringe vectors")
        else:
            logger.debug("Attempting fringe correction")
            reduction_log.attempt('fringe')
            fringe_filename = mastercals.fringe(mjd=mjd)
            fringe_vector_file = mastercals.fringevector(ota=ota)
            reduction_files_used['fringemap'] = fringe_filename
            reduction_files_used['fringevector'] = fringe_vector_file

            guide_ota = (hdu.header['CELLMODE'].find("V") >= 0) or is_guide_ota(primhdu=None, ext=hdu, data=merged)
            fringe_hdu = pyfits.open(fringe_filename)
            for ext in fringe_hdu[1:]:
                if (extname == ext.header['EXTNAME']):
                    if (options['verbose']): print("Working on fringe scaling for",extname)

                    if (guide_ota):
                        logger.info("Ignoring guide-OTA during fringe scaling")
                        reduction_log.success('fringe')
                    else:
                        fringe_scaling = podi_fringing.get_fringe_scaling(merged, ext.data, fringe_vector_file, extname=extname)

                    #  numpy.savetxt("fringe_%s" % (extname), fringe_scaling)
                    data_products['fringe-template'] = ext.data
                    if (fringe_scaling is not None):
                        # numpy.savetxt("fringescaling_%s.txt" % (extname), fringe_scaling)
                        good_scalings = three_sigma_clip(fringe_scaling[:,6], [0, 1e9])
                        fringe_scaling_median = numpy.nanmedian(good_scalings)
                        fringe_scaling_std    = numpy.nanstd(good_scalings)
                    break
            hdu.header.add_history("fringe map: %s" % fringe_filename)
            hdu.header.add_history("fringe vector: %s" % fringe_vector_file)
            fringe_hdu.close()

            # Do not determine fringe scaling if this OTA is marked as video cell
            if (guide_ota):
                pass
            elif (numpy.isfinite(fringe_scaling_median) and numpy.isfinite(fringe_scaling_std)):
                hdu.header["FRNG_SCL"] = fringe_scaling_median
                hdu.header["FRNG_STD"] = fringe_scaling_std
                hdu.header["FRNG_OK"] = (fringe_scaling is not None)
            else:
                logger.error("Unable to determine fringe scaling (encountered too many NaN values")
                reduction_log.fail('fringe')

        # Insert the DETSEC header so IRAF understands where to put the extensions
        start_x = ota_c_x * size_x #4096
        start_y = ota_c_y * size_y #4096        
        end_x = start_x + size_x #det_x2 - det_x1
        end_y = start_y + size_y #det_y2 - det_y1
        detsec_str = "[%d:%d,%d:%d]" % (start_x, end_x, start_y, end_y)
        hdu.header["DETSEC"] = (detsec_str, "position of OTA in focal plane")
                
        hdu.header['LTM1_1'] = 1.0
        hdu.header['LTM2_2'] = 1.0
        hdu.header['LTV1'] = -start_x
        hdu.header['LTV2'] = -start_y

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
        if (options['target_coords'] is not None):
            if (options['verbose']): print("Overwriting target positions",options['target_coords'])
            ra, dec = options['target_coords']
            coord_j2000 = ephem.Equatorial(ra, dec, epoch=ephem.J2000)
            if (options['verbose']): print("\n\nAdjusting ra/dec:", coord_j2000.ra, coord_j2000.dec, "\n\n")

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
            offset_total = numpy.array([0.,0.], dtype=numpy.float32)
            if (options['offset_pointing'] != None):
                offset_total += numpy.array(options['offset_pointing'])
                logger.debug("Adding offset to Ra/Dec: %s --> %s" % (str(numpy.array(options['offset_pointing'])), str(offset_total)))
            if (options['offset_dither'] != None):
                offset_total += numpy.array(options['offset_dither'])
            logger.debug("Adding user's WCS offset (ra: %f, dec: %f degrees)" % (offset_total[0], offset_total[1]))
            hdu.header['CRVAL1'] += offset_total[0] 
            hdu.header['CRVAL2'] += offset_total[1]


        # Now add the canned WCS solution
        if (not mastercals.apply_wcs()): #['wcs_distortion'] == None):
            reduction_log.not_selected('wcs_dist')
        elif (mastercals.wcs(mjd) == "plain"):
            #
            # If there's no WCS system at all, at least put in something semi-reasonable
            #
            logger.warning("Adding plain WCS as fall-back solution")
            for a,b in itertools.product(range(10), repeat=2):
                keyname = "WAT%d_%03d" % (a,b)
                if (keyname in hdu.header):
                    del hdu.header[keyname]
            for a,b in itertools.product(range(2), range(100)):
                keyname = "PV%d_%d" % (a,b)
                if (keyname in hdu.header):
                    del hdu.header[keyname]
            hdu.header['CRPIX1'] = (4-ota_c_x)*4500
            hdu.header['CRPIX2'] = (4-ota_c_y)*4500
        elif (mastercals.wcs(mjd) is None):
            reduction_log.fail("wcs_dist")
        else:
            wcs_model_fn = mastercals.wcs(mjd)
            logger.debug("Applying wcs from %s" % (str(wcs_model_fn)))
            wcsdistort = apply_wcs_distortion(wcs_model_fn, hdu, binning, reduction_log)
            reduction_files_used['wcs'] = wcs_model_fn
            if (options['simple-tan-wcs']):
                hdu.header['CTYPE1'] = "RA---TAN"
                hdu.header['CTYPE2'] = "DEC--TAN"
            
        #
        # If requested, perform cosmic ray rejection
        #
        if (options['crj'] <= 0):
            reduction_log.not_selected('crj')
        else:
            logger.debug("Starting cosmic ray removal with %d iterations" % (options['crj']))
            corrected, mask = podi_cosmicrays.remove_cosmics(data=merged, 
                                                             n_iterations=options['crj'],
                                                             gain=gain_avg,
                                                             readnoise=rone_avg,
                                                             sigclip=options['crj_sigclip'], 
                                                             sigfrac=options['crj_sigfrac'],
                                                             objlim=options['crj_objlim'],
                                                             saturation_limit=options['crj_saturation'],
                                                             binning=binning,
                                                             verbose=False,)
            merged = corrected
            reduction_log.success('crj')
            logger.debug("Done with cosmic ray removal")
        hdu.header['CRJ_ITER'] = (options['crj'], "cosmic ray removal iterations")
        hdu.header['CRJ_SIGC'] = (options['crj_sigclip'], "CR sigma clipping threshold")
        hdu.header['CRJ_SIGF'] = (options['crj_sigfrac'], "CR threshold for neighboring pixels")
        hdu.header['CRJ_OBJL'] = (options['crj_objlim'], "CR contrast to underlying object")
        hdu.header['CRJ_SATU'] = (options['crj_saturation'], "CR saturation limit")
        hdu.header['CRJ_METH'] = (options['crj_method'], "CR implementation")

        # Insert the new image data. This also makes sure that the headers
        # NAXIS, NAXIS1, NAXIS2 are set correctly
        hdu.data = merged


        if (not mastercals.apply_photflat()):
            reduction_log.not_selected('photflat')
        elif (mastercals.photflat() is None):
            reduction_log.fail('photflat')
        else:
            photflat_filename = mastercals.photflat()
            photflat = pyfits.open(photflat_filename)
            reduction_files_used['photflat'] = photflat_filename

            # ADD HERE
            # subtract the background before the photflat correction and re-add afterwards to avoid
            # intruducing too much structure into the sky background

            # Search for the flatfield data for the current OTA
            all_success = True
            try:
                photflat_ext = photflat[extname]

                try:
                    merged /= photflat_ext.data
                    logger.info("Applying photometric flatfield from %s to %s" % (photflat_filename, extname))
                except:
                    logger.warning("Unable to apply flat-field, dimensions don't match (data: %s, flat: %s)" % (
                        str(merged.shape), str(ff_ext.data.shape)))
                    reduction_log.partial_fail('photflat')
                    all_success = False
            except KeyError:
                logger.warning("Unable to find photometric flatfield data for %s in %s" % (ota_name, photflat_filename))
                reduction_log.partial_fail('photflat')
                all_success = False
                pass
            if (all_success):
                reduction_log.success('photflat')

        # #
        # # Return data as shared memory
        # #
        # # Allocate shared memory
        # shmem_buffer = multiprocessing.RawArray(ctypes.c_float, merged.size)
        # shmem_image = shmem_as_ndarray(shmem_buffer).reshape(merged.shape)
        # # insert reduced image data into shared memory buffer
        # shmem_image[:,:] = merged[:,:]

        # #
        # # Now set HDU properties and store the shared memory image
        # #
        # hdu.header['NAXIS1'] = merged.shape[1]
        # hdu.header['NAXIS2'] = merged.shape[0]
        # hdu.data = None #shmem_buffer
        # data_products['shmem_image'] = shmem_buffer

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
                logger.debug("Preparing source catalog")
                tmphdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdu.header, data=merged)])
                obsid = tmphdulist[0].header['OBSID']
                process_id = os.getpid()
                fitsfile = "%s/tmp.pid%d.%s_OTA%02d.fits" % (sitesetup.sextractor_cache_dir, process_id, obsid, ota)
                catfile = "%s/tmp.pid%d.%s_OTA%02d.cat" % (sitesetup.sextractor_cache_dir, process_id, obsid, ota)
                tmphdulist.writeto(fitsfile, overwrite=True)
                logger.debug("Wrote temp file to %s" % (fitsfile))
                sex_config_file = "%s/config/wcsfix.sex" % (sitesetup.exec_dir)
                parameters_file = "%s/config/wcsfix.sexparam" % (sitesetup.exec_dir)
                catalog_format = "-CATALOG_TYPE %s" % ("FITS_LDAC" if options['sextractor_write_fits'] else "ASCII_HEAD")
                sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s %s" % (
                    sitesetup.sextractor, sex_config_file, parameters_file,
                    catfile, catalog_format,
                    fitsfile)
                if (options['verbose']): print(sexcmd)

                start_time = time.time()
                #os.system(sexcmd)
                for sex_restarts in range(3):
                    try:
                        logger.debug("Running SourceExtractor (attempt %d)" % (sex_restarts+1))
                        ret = subprocess.Popen(sexcmd.split(),
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
                        sextractor_pid = ret.pid
                        #
                        # Wait for sextractor to finish (or die)
                        #
                        ps = psutil.Process(sextractor_pid)
                        sextractor_error = False
                        while(True):
                            try:
                                ps_status = ps.status()
                                if (ps_status in [psutil.STATUS_ZOMBIE,
                                                  psutil.STATUS_DEAD] and
                                    ret.poll() is not None):
                                    logger.critical("Sextractor died unexpectedly (%s - this is try #%d / pid=%d)" % (
                                        str(ps_status), sex_restarts+1, sextractor_pid))
                                    sextractor_error = True
                                    break
                            except psutil.NoSuchProcess:
                                pass

                            if (ret.poll() is None):
                                # sextractor completed
                                logger.debug("SourceExtractor completed successfully!")
                                break
                            time.sleep(0.1)

                        if (not sextractor_error):
                            (sex_stdout, sex_stderr) = ret.communicate()
                            if (ret.returncode != 0):
                                logger.warning("Sextractor might have a problem, check the log")
                                logger.debug("Stdout=\n"+sex_stdout)
                                logger.debug("Stderr=\n"+sex_stderr)
                            break

                    except OSError as e:
                        podi_logging.log_exception()
                        print >>sys.stderr, "Execution failed:", e
                    end_time = time.time()
                    logger.debug("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
                #
                # By now we have run Sextractor, next up is loading the catalog
                #
                try:
                    source_cat = None
                    try:
                        if (options['sextractor_write_fits']):
                            logger.debug("Reading Sextractor FITS catalog %s" % (catfile))
                            source_cat = podi_readfitscat.read_fits_catalog(catfile, 'LDAC_OBJECTS', flatten=True)
                            #logger.info("%s: Found %d sources" (extname, (-1 if source_cat is None else source_cat.shape[0]) ))
                        else:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore")
                                logger.debug("Reading Sextractor ASCII catalog (%s)" % (catfile))
                                source_cat = numpy.loadtxt(catfile)
                    except IOError:
                        logger.warning("The Sextractor catalog is empty, ignoring this OTA")
                        source_cat = None
                    else:
                        if (source_cat.shape[0] == 0 or source_cat.ndim < 2):
                            source_cat = None
                        else:
                            source_cat[:,SXcolumn['ota']] = ota

                            valid_coords = \
                                (source_cat[:, SXcolumn['x']] > 0) & \
                                (source_cat[:, SXcolumn['x']] < merged.shape[1]) & \
                                (source_cat[:, SXcolumn['y']] > 0) & \
                                (source_cat[:, SXcolumn['y']] < merged.shape[0])
                            source_cat = source_cat[valid_coords]
                            n_bad = numpy.sum(~valid_coords)
                            if (n_bad > 0):
                                logger.debug("Excluded %d sources with invalid center positions" % (n_bad))

                            # also convert all physical (pixel) FWHMs into arcseconds
                            source_cat[:, SXcolumn['fwhm_world']] = source_cat[:, SXcolumn['fwhm_image']] * 0.11
                            # print(source_cat[:, SXcolumn['fwhm_world']])

                            flags = source_cat[:,SXcolumn['flags']]
                            no_flags = (flags == 0)
                            logger.debug("Found %d sources, %d with no flags" % (source_cat.shape[0], numpy.sum(no_flags)))
                except:
                    source_cat = None
                #numpy.savetxt("ota%d" % (ota), source_cat)
                if (sitesetup.sex_delete_tmps):
                    clobberfile(fitsfile)
                    clobberfile(catfile)
                    pass


            fixwcs_data = None

        #
        # Create some PSF quality data
        #
        psf = None
        if (source_cat is not None):
            start = time.time()
            psf = psf_quality.PSFquality(
                catalog_filename=None,
                pixelscale=0.11,
                catalog=source_cat,
                image_data=merged,
                use_vignets=False,
                detector=ota,
            )
            psf.info(logger=logger)
            end = time.time()
            # psf.save2fits(fn="psf_%02d.fits" % (psf.detector))
            logger.debug("Spent %.3f seconds creating PSF model" % (end-start))

        #
        # Sample that background at random place so we can derive a median background level later on
        # This is necessary both for the pupil ghost correction AND the fringing correction
        #
        starcat = None
        if (source_cat is not None):
            ota_x, ota_y = source_cat[:,2], source_cat[:,3]
            starcat = (ota_x, ota_y)
        # Now sample the background, excluding regions close to known sources
        logger.debug("Sampling sky background")
        sky_samples = numpy.array(podi_fitskybackground.sample_background(data=merged, wcs=None, 
                                                                          starcat=starcat, 
                                                                          min_found=200, boxwidth=30,
                                                                          min_box_spacing=3))

        if (sky_samples.shape[0] <= 0): #not sky_samples.shape[1] > 4):
            logger.warning("Something went wrong with the sky-calculation: %s" % (str(sky_samples.shape)))
            # print "sky-samples =",sky_samples
            # print "sky-samples.shape =",sky_samples.shape
            sky_level_median, sky_level_mean, sky_level_std = -1, -1, -1
            sky_samples = None
        else:
            sky_level_median = numpy.median(sky_samples[:,4])
            sky_level_mean   = numpy.mean(sky_samples[:,4])
            sky_level_std    = numpy.std(sky_samples[:,4])
            hdu.header["SKY_MEDI"] = (sky_level_median, "sky-level median")
            hdu.header["SKY_MEAN"] = (sky_level_mean, "sky-level mean")
            hdu.header["SKY_STD"] = (sky_level_std, "sky-level rms")
        logger.debug("Found median sky-level: %d" % (sky_level_median))

    pupilghost_scaling = None
    pupilghost_template = None
    filter_level = get_filter_level(hdulist[0].header)
    filter_name = get_valid_filter_name(hdulist[0].header)
    if (not mastercals.apply_pupilghost()):
        logger.debug("PG not selected")
        reduction_log.not_selected('pupilghost')
    elif (False): #not 'PGAFCTD' in hdu.header or not hdu.header['PGAFCTD']):
        # This frame does not contain the keyword labeling it as affected by
        # the pupilghost. In that case we don't need to do anything
        logger.debug("This extension (%s) does not have any pupilghost problem" % (extname))
        reduction_log.not_required('pupilghost')
    elif (mastercals.pupilghost() is None):
        reduction_log.missing_data('pupilghost')
        logger.critical("Pupilghost correction requested, but no template found")
    else:
        logger.debug("Getting ready to subtract pupil ghost from science frame")
        # binning = ota_list[0].header['BINNING']
        pg_template = mastercals.pupilghost()
        logger.debug("looking for pupil ghost template %s..." % (pg_template))

        # If we have a template for this level
        logger.debug("\n   Using pupilghost template %s, filter %s ... " % (pg_template, filter_name))
        pg_hdu = pyfits.open(pg_template)

        # Getting the pupilghost scaling factors for this OTA
        reduction_files_used['pupilghost'] = pg_template
        # print "redfiles used:", reduction_files_used

        #
        # Find center position of the pupilghost. Assume that the science
        # frames have the same pupilghost centering as the flatfield used for
        # calibration, as this provides a more reliable center position.
        #
        if (pupilghost_center_x is None or pupilghost_center_y is None):
            pgcenter = 'data'
        else:

            logger.debug("Using PG center position: %f %f" % (
                pupilghost_center_x, pupilghost_center_y))

            pgcenter = numpy.array([pupilghost_center_x, pupilghost_center_y])
            if (numpy.isnan(pupilghost_center_x) or numpy.isnan(pupilghost_center_y)):
                logger.info("At least one of the center positions in NaN, setting to 'data' instead")
                pgcenter = 'data'

        # Compute the pupilghost image for this OTA at the right orientation
        logger.debug("Starting pg scaling")
        pupilghost_template = podi_matchpupilghost.compute_pupilghost_template_ota(
            hdu, pg_hdu,
            rotate=True,
            non_negative=True,
            source_center_coords=pgcenter
        )

        data_products['pupilghost-template'] = pupilghost_template
        if (pupilghost_template is not None):
            # print "PG-template:", extname, "\n", pupilghost_template
            _, pupilghost_scaling = podi_matchpupilghost.get_pupilghost_scaling_ota(
                science_hdu=hdu,
                pupilghost_frame=pg_hdu, #pupilghost_template,
                n_samples=750, boxwidth=20,
                verbose=False,
                pg_matched=True,
                return_all=True,
            )
            # print pupilghost_scaling
            data_products['pupilghost-scaling'] = pupilghost_scaling
            logger.debug("PG scaling:\n%s" % (str(pupilghost_scaling)))
            reduction_log.attempt('pupilghost')
            if (pupilghost_scaling is None):
                logger.info("No PG samples found for OTA %s" % (extname))
            else:
                logger.info("Found %d PG scaling samples" % (pupilghost_scaling.shape[0]))
        else:
            logger.debug("Could not compute PG template for OTA %s" % (extname))
            reduction_log.no_data('pupilghost')
        logger.debug("Done with pg scaling")

        # # Find the optimal scaling factor
        # logger.debug("Searching for optimal pupilghost scaling factor")
        # any_affected, scaling, scaling_std = podi_matchpupilghost.get_pupilghost_scaling(ota_list, pg_hdu)
        # logger.debug("Check if any OTA is affected: %s" % ("yes" if any_affected else "no"))
        # logger.debug("Optimal scaling factor found: %.2f +/- %.2f" % (scaling, scaling_std))

    
    # Write the reduction log
    reduction_log.write_to_header(hdu.header)

    # if one of the following headers exist, add a title before them
    for key in ['COMMENT', 'HISTORY']:
        if (key in hdu.header):
            add_fits_header_title(hdu.header, "other information", key)
            break

    data_products['hdu'] = hdu
    data_products['wcsdata'] = None #fixwcs_data
    data_products['sky-samples'] = sky_samples
    data_products['sky'] = (sky_level_median, sky_level_mean, sky_level_std)
    data_products['fringe_scaling'] = fringe_scaling
    data_products['sourcecat'] = source_cat
    data_products['tech-header'] = tech_header
    data_products['reduction_files_used'] = reduction_files_used
    data_products['psf'] = psf

    logger.debug("Done with collect_cells for this OTA")
    return data_products #hdu, fixwcs_data
    

















def apply_software_binning(hdu, softbin):

    logger = logging.getLogger("Softbin")

    logger.debug("Applying software binning of x%d to ext %s" %
                (softbin, hdu.name))

    # for software binning, bin data
    data = hdu.data
    # data_4d = data.reshape(data.shape[0]/softbin, softbin, softbin, data.shape[1]/softbin)
    # data_binned_3 = bottleneck.nanmean(data_4d, axis=-1)
    # data_binned_2 = bottleneck.nanmean(data_binned_3, axis=1)
    # hdu.data = data_binned_2
    hdu.data = rebin_image(data, softbin, operation=bottleneck.nanmean)

    # also modify all WCS relevant headers
    for hdrname in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
        hdu.header[hdrname] *= softbin

    for hdrname in [
            'CRPIX1', 'CRPIX2',
            ]:
        hdu.header[hdrname] /= softbin

    # also modify the DETSEC keyword so we can properly display the frame as 
    # IRAF mosaic
    detsec_str = hdu.header['DETSEC'][1:-1] # strip the [ and ]
    n = []
    for xy in detsec_str.split(","):
        for from_to in xy.split(":"):
            n.append(int(from_to))
    n = numpy.array(n) / softbin
    new_detsec = "[%d:%d,%d:%d]" % (n[0],n[1],n[2],n[3])
    logger.debug(new_detsec)
    hdu.header['DETSEC'] = new_detsec 

    hdu.header['SOFTBIN'] = (softbin, "software binning")

    #
    # Add here: appropriate treatment for photometric zeropoint, etc
    #

    return hdu




def parallel_collect_reduce_ota(workqueue,
                                final_results_queue,
                                intermediate_queue,
                                intermediate_results_queue,
                                filename, ota_id,
                                intermediate_ack_queue,
                                options=None,
                                shmem=None,
                                shmem_dim=None,
                                shmem_id=-1,
                                quit_signal=None,
                                complete=None,
                                ):
    """
    A minimal wrapper handling the parallel execution of collectcells.

    Parameters
    ----------
    workqueue : Queue

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

    logger = logging.getLogger("WorkManager")
    logger.debug("Process started")

    # while (True):
    #     task = workqueue.get()
    #     # if (task == None):
    #     #     logger.info("Received termination signal, shutting down")
    #     #     workqueue.task_done()
    #     #     return

    #     filename, ota_id, wrapped_pipe = task
    for _x in range(1):
        logger = logging.getLogger("WorkManager(OTA-ID:%02d)" % (ota_id))

        # Do the work
        try:
            data_products = collect_reduce_ota(filename, options=options)
        # except (KeyboardInterrupt, SystemExit):
        #     workqueue.task_done()
        #     while (not workqueue.empty()):
        #         workqueue.get()
        #         workqueue.task_done()
        #     break
        except:
            # All other problems:
            # Note them in the debug log, and raise the error for further 
            # processing
            podi_logging.log_exception("parallel_collectcells")
            raise


        return_hdu = data_products['hdu']
        # print data_products
        if (return_hdu is None):
            # workqueue.task_done()
            # continue
            if (not os.path.isfile(filename)):
                logger.critical("OTA did not return any data, missing file (%s)?" % (filename))
            else:
                logger.critical("OTA did not return any data, empty/faulty file (%s)?" % (filename))

            #
            # Send error message back via intermediate result_queue back to master 
            #

            # return_queue.put( (ota_id, data_products, shmem_id) )

            #workqueue.task_done()
            #continue
            return


        extname = return_hdu.header['FPPOS'] if (return_hdu is not None and 'FPPOS' in return_hdu.header) else "???"
        ota = return_hdu.header['OTA'] if (return_hdu is not None and 'OTA' in return_hdu.header) else -1
        logger.debug("Received OTA pre-processed data for OTA %s" % (extname))

        logger = logging.getLogger("PostProc.OTA%02d" % (ota))
        logger.debug("Trimming off pupilghost template and fringe template")

        # Trim the data section of the return data to keep transfer delays low
        pg_image = None
        try:
            pg_image = data_products['pupilghost-template']
            del data_products['pupilghost-template']
        except:
            pass

        del data_products['hdu']

        fringe_template = None
        try:
            fringe_template = data_products['fringe-template']
            del data_products['fringe-template']
        except:
            pass

        # However, we do need the headers for the intermediate processing steps
        logger.debug("Adding back header")
        data_products['header'] = return_hdu.header

        # Send the results from this OTA to the main process handler
        logger.debug("Sending results back to main process")
        try:
            intermediate_results_queue.put((ota_id, data_products, shmem_id), block=True )
        except AssertionError:
            logger.warning("assertion error / queue closed: intermediate_results_queue")
            pass

        # Now unpack the communication pipe
        logger.debug("Preparing communication pipe ...")
        #logger.debug(str(wrapped_pipe))
        #fct, params = wrapped_pipe
        #pipe = fct(*params)

        #
        # Wait to hear back with the rest of the instructions
        #
        logger.debug("Waiting to hear back with fringe/pupilghost scaling")
        while (not quit_signal.value):
            # wait for instruction for my OTA-ID
            try:
                ret = intermediate_queue.get(timeout=0.01)
            except queue.Empty:
                continue
            except IOError:
                # most likely due to closed pipe after the queue has been closed
                logger.warning("IO error / queue closed: intermediate_queue")
                break
            except:
                podi_logging.log_exception()

            # print ret, " // ", ota_id
            try:
                _ota_id, final_parameters = ret
            except:
                podi_logging.log_exception()
                print ("\n\n",ret,"\n\n")
                
            if (not _ota_id == ota_id):
                # this is not meant for me, so put it back
                try:
                    intermediate_queue.put(ret) #(_ota_id, final_parameters))
                except AssertionError:
                    # AssertionError most likely happens when the queue is 
                    # already closed --> shut-down if we encounter this problem
                    return

                # time.sleep(0.1)
            else:
                # got what I need
                break

        # acknowledge that we received the intermediate data
        logger.debug("Sending ACK message")
        try:
            intermediate_ack_queue.put(ota_id)
        except AssertionError:
            # this means the queue is closed - but since it's closed in a 
            # different process and that isn't reflected here, that likely
            # won't happen
            logger.warning("assertion error / queue closed: intermediate_ack_queue")
            pass


        # r = int(time.time()*1e6)
        # logger.info("RND: %d" % (r))
        # if (r%10 <= 2):
        #     logger.error("committing suicide by random chance (pid: %d)" % (os.getpid()))
        #     #sys.exit(1)
        #     return

        #_final_parameters = pipe.recv()
        #logger.info("OTA-ID %02d received PIPE final parameters:\n%s" % (ota_id, _final_parameters))
        
        logger.debug("OTA-ID %02d received final parameters: %s" % (ota_id, final_parameters))

        # XXX HACK
        #if (_ota_id == 22):
        #    return

        #
        # Finish work: fringe subtraction and pupilghost removal
        #

        reduction_log = data_products['reduction_log']

        #
        # XXX: UNPACK SHMEM AND REPACK WHEN DONE
        #

        logger.debug("Next up (maybe): fringe subtraction")
        try:
            # Remove fringing by subtracting the scaled fringe template
            if (not type(fringe_template) == type(None) and final_parameters['fringe-scaling-median'] > 0):
                logger.debug("Subtracting fringes (%.2f)..." % (final_parameters['fringe-scaling-median']))
                logger.debug("data: %s // fringe: %s" % (str(type(return_hdu.data)), str(type(fringe_template))))
                logger.debug("data: %s //  fringe: %s" % (str(return_hdu.data.shape), str(fringe_template.shape)))
                fringe_template *= final_parameters['fringe-scaling-median']
                logger.debug("done computing scaled fringe correction")
                return_hdu.data -= fringe_template
                logger.debug("completed fringe subtraction")
                reduction_log.success('fringe')
        except:
            reduction_log.fail('fringe')
            podi_logging.log_exception()
            pass

        logger.debug("Next up (maybe): pupilghost subtraction")
        try:
            # Also delete the pupilghost contribution
            if (pg_image is not None and final_parameters['pupilghost-scaling-median'] > 0):
                logger.debug("Subtracting pupilghost (%.2f)..." % (final_parameters['pupilghost-scaling-median']))
                return_hdu.data -= (pg_image * final_parameters['pupilghost-scaling-median'])
                reduction_log.success('pupilghost')
            else:
                reduction_log.fail('pupilghost')
        except:
            podi_logging.log_exception()
            pass

        ##################
        #
        # Now that actual reduction is done, apply software binning
        #
        ##################
        logger.debug("Next step (optional): software binning")
        if (not options['softbin'] == 0):
            # check if its a multiple of 2
            if (options['softbin'] in [2,4,8]):
                return_hdu = apply_software_binning(return_hdu, options['softbin'])
                reduction_log.success('softwarebin')
            else:
                reduction_log.fail('softwarebin')
        else:
            reduction_log.not_selected('softwarebin')

        #
        #
        # Pack the image data into shared memory to bring down communication times
        #
        #
        wx,wy = shmem_dim
        logger.debug("MPC SHMEM: %s %d,%d" % (str(shmem), wx, wy))
        shmem_image = shmem.to_ndarray() #shmem_as_ndarray(shmem).reshape((wy,wx))
        logger.debug("Packing image return into shared memory (%d x %d)" % (wx, wy))
        if (return_hdu.data.shape[1] > wx or
            return_hdu.data.shape[0] > wy):
            logger.critical("Image size exceeds allocated shared memory")
        else:
            shmem_image[:return_hdu.data.shape[0], :return_hdu.data.shape[1]] = return_hdu.data[:,:]
            return_hdu.data = numpy.array(return_hdu.data.shape)

        # Add the complete ImageHDU to the return data stream
        data_products['hdu'] = return_hdu

        # Add the results to the return_queue so the master process can assemble the result file
        logger.debug("Adding results for OTA %02d to return queue" % (ota_id))
        try:
            final_results_queue.put( (ota_id, data_products, shmem_id) )
        except AssertionError:
            logger.error("Unable to put final results into final_results_queue!")
            pass
        #time.sleep(0.1)
        # pipe.close()

        # podi_logging.print_stacktrace(logger=logger)

        #cmd_queue.task_done()
        #workqueue.task_done()
        if (complete is not None):
            complete.value = True
            pass
        logger.debug("Done with work, shutting down!")
        break

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
                              timeout=None):


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
    
    logger = logging.getLogger("CC_w_timeout")

    cc_args = (input, outputfile,
               batchmode,
               verbose,
               options,
               showsplash)

    p = multiprocessing.Process(target=collectcells, args=cc_args)
    logger.debug("Starting collectcells with timeout")
    p.start()
    logger.debug("collectcells started!")

    start_time = time.time()
    timeout = timeout if timeout > 0 else None
    p.join(timeout)
    end_time = time.time()
    logger.debug("returned from joining CC after %.3f seconds" % ((end_time-start_time)))

    if (p.is_alive()):
        logger.warning("Timeout event triggered, shutting things down ...")
        #kill_all_child_processes(process_tracker)

        logger.info("Killing collectcells after timeout...")
        # podi_logging.print_stacktrace()
        p.terminate()
        logger.info("all done after timeout problem/error!")
        return 1

    return 0




def my_system(cmd):
    logger = logging.getLogger("RunShellCmd")
    start_time = time.time()
    returncode = None
    try:
        logger.debug("Running %s" % (cmd))
        ret = subprocess.Popen(cmd.split(), 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE)
        (stdout, stderr) = ret.communicate()
        end_time = time.time()
        returncode = ret.returncode
        logger.debug("Execution done after %.3f seconds" % (end_time-start_time))
        # if (ret.returncode != 0):
        #     logger.warning("Sextractor might have a problem, check the log")
        #     logger.debug("Stdout=\n"+sex_stdout)
        #     logger.debug("Stderr=\n"+sex_stderr)
    except OSError as e:
        #podi_logging.log_exception()
        logger.debug("Error while executing: Err:%d - %s" % (e.errno, e.strerror))
    except:
        podi_logging.log_exception()
        pass
        
    return


def prestage_data(options, input):

    logger = logging.getLogger("PreStage")
    import shutil

    logger.info("Prestaging data to staging dir...")
    staged = True
    procs = []
    if (os.path.isdir(input)):
        logger.debug("Input is a directory: %s" % (input))

        # This is a directory.
        # Let's assume all files in this directory need to be prestaged
        base, dirname = os.path.split(os.path.abspath(input))
        print(base, dirname)
        tmpdir = "%s/%s" % (sitesetup.staging_dir, dirname)
        # Create the directory
        if (not os.path.isdir(tmpdir)):
            os.mkdir(tmpdir)
        # Now copy all .fits and .fits.fz files into the new directory
        for fn in os.listdir(input):

            filename = "%s/%s" % (input, fn)
            if (not os.path.isfile(filename)):
                # Not a file (probaby a directory) -> ignore it
                continue

            if (fn.endswith(".fits.fz")):
                # This is a .fits.fz file -> run funpack
                outfile = "%s/%s" % (tmpdir, fn[:-3])
                if (os.path.isfile(outfile)):
                    continue
                logger.debug("Prestaging file w/ funpack: %s --> %s" % (filename, outfile))
                # os.system("funpack -O %s %s" % (outfile, filename))
                p = multiprocessing.Process(target=my_system,
                                            kwargs={"cmd": "funpack -O %s %s" % (outfile, filename)}
                                        )
                p.start()
                procs.append(p)
            elif (fn.endswith(".fits")):
                # This a fits file -> simply copy it
                outfile = "%s/%s" % (tmpdir, fn)
                if (os.path.isfile(outfile)):
                    continue
                logger.debug("Prestaging file w/ copy: %s --> %s" % (filename, outfile))
                #shutil.copyfile(filename, outfile)
         
            
    elif (os.path.isfile(input)):
        logger.debug("input is filename")
        # its a file
        # find the basename of the file
        dirname, filename = os.path.split(os.path.abspath(input))
        if (filename.endswith(".fits.fz")):
            base = filename[:-11]
        elif (filename.endswith(".fits")):
            base = filename[:-8]
        else:
            # not sure what this input is
            return

        tmpdir = "%s/%s" % (sitesetup.staging_dir, base)
        # Create the directory
        if (not os.path.isdir(tmpdir)):
            os.mkdir(tmpdir)
        # Now copy all .fits and .fits.fz files into the new directory
        for fn in os.listdir(dirname):

            filename = "%s/%s" % (dirname, fn)
            if (not os.path.isfile(filename) or not fn.beginswith(base)):
                # Not a file (probaby a directory) -> ignore it
                continue
            
            if (fn.endswith(".fits.fz")):
                # This is a .fits.fz file -> run funpack
                outfile = "%s/%s" % (tmpdir, fn[:-3])
                print(outfile)
                if (os.path.isfile(outfile)):
                    continue
                logger.debug("Prestaging file w/ funpack: %s --> %s" % (filename, outfile))
                # os.system("funpack -O %s %s" % (outfile, filename))
                p = multiprocessing.Process(target=my_system,
                                            kwargs={"cmd": "funpack -O %s %s" % (outfile, filename)}
                                        )
                p.start()
                procs.append(p)
            elif (fn.endswith(".fits")):
                # This a fits file -> simply copy it
                outfile = "%s/%s" % (tmpdir, fn)
                if (os.path.isfile(outfile)):
                    continue
                logger.debug("Prestaging file w/ copy: %s --> %s" % (filename, outfile))
                #shutil.copyfile(filename, outfile)

    else:
        logger.error("Input (%s) is neither path nor file, can't prestage this!\n" % (input))
        staged = False

    for p in procs:
        p.join()

    logger.info("Done prestaging data !")
    
    return tmpdir, staged


def unstage_data(options, staged, input):

    logger = logging.getLogger("Unstage")

    if (options['prestage'] and staged):
        if (input.startswith(sitesetup.staging_dir)):
            files = os.listdir(input)
            for f in files:
                fn = "%s/%s" % (input, f)
                if (fn.endswith(".fits")):
                    logger.debug("deleting staged file: %s" % (fn))
                    os.remove(fn)
            try:
                os.rmdir(input)
                logger.debug("removing staging directory %s" % (input))
            except:
                logger.warning("Unable to remove staging directory %s" % (input))
                pass

    return


        










import pickle

class reduce_collect_otas (object):

    def __init__(self, options, number_cpus):

        self.options = options
        self.number_cpus = number_cpus
        self.quit = False

        self.qr_quit = multiprocessing.Value('i', False)

        self.logger = logging.getLogger("QRWorker")

        # self.logger.info("stacktrace before starting qr worker:")
        # podi_logging.print_stacktrace(logger=self.logger)

        self.info = []

        #
        # Setup all queues;
        # for each queue, start and rename the internal thread
        #
        self.queue = multiprocessing.JoinableQueue()
        self.queue._start_thread()
        self.queue._thread.name = "QueueFeederThread__JobQueue"

        self.final_results_queue = multiprocessing.Queue()
        self.final_results_queue._start_thread()
        self.final_results_queue._thread.name = "QueueFeederThread__FinalResultQueue"

        self.intermediate_queue = multiprocessing.Queue()
        self.intermediate_queue._start_thread()
        self.intermediate_queue._thread.name = "QueueFeederThread__IntermediateDataQueue"

        self.intermediate_results_queue = multiprocessing.Queue()
        self.intermediate_results_queue._start_thread()
        self.intermediate_results_queue._thread.name = "QueueFeederThread__IntermediateResultsQueue"

        self.intermediate_data_ack_queue = multiprocessing.Queue()
        self.intermediate_data_ack_queue._start_thread()
        self.intermediate_data_ack_queue._thread.name = "QueueFeederThread__AcknIntermedDataReceivedQueue"

        self.shmem_dims = (4096, 4096)


        self.intermediate_results_done = multiprocessing.Lock()
        self.intermediate_results_done.acquire()
        self.intermediate_results_queue_lock = multiprocessing.Lock()

        self.final_results_done = multiprocessing.Lock()
        self.final_results_done.acquire()
        self.final_results_queue_lock = multiprocessing.Lock()
        self.id2filename = {}
        self.filename2id = {}

        self.job_status_lock = multiprocessing.Lock()

        self.files_to_reduce = []
        
        self.active_workers = 0
        self.all_closed = False

        #
        # Start feeding the workers
        #
        self.feed_worker_thread = threading.Thread(
            target=self.feed_workers,
            name="QRThread__FeedWorkers",
            )

        #
        # Also start a worker to collect intermediate results, send 
        # intermediate data and handle data receipt acknowledgements
        #
        self.collect_intermediate_results_thread = threading.Thread(
            target=self.collect_intermediate_results,
            name="QRThread__CollectIntermediateResults",
            )
        self.collect_intermediate_data_broadcast_thread = threading.Thread(
            target=self.broadcast_intermediate_data,
            name="QRThread__BroadcastIntermediateData",
            )
        self.acknowledge_intermediate_data_thread = threading.Thread(
            target=self.acknowledge_intermediate_data_received,
            name="QRThread__AcknowledgeIntermediateDataReceived",
            )
        self.intermediate_results_complete = False
        self.intermediate_data_back_to_workers = None
        self.intermediate_results = []

        self.final_results = []
        self.collect_final_results_thread = threading.Thread(
            target=self.collect_final_results,
            name="QRThread__CollectFinalResults"
            )

        # Make all threads daemons - that way they don't keep the program from
        # not shutting down
        self.feed_worker_thread.setDaemon(True)
        self.collect_intermediate_results_thread.setDaemon(True)
        self.collect_intermediate_data_broadcast_thread.setDaemon(True)
        self.acknowledge_intermediate_data_thread.setDaemon(True)
        self.collect_final_results_thread.setDaemon(True)
        
        self.intermediate_results_complete = False

        self.shmem_list = {}

        self.intermed_results_sent = 0
        pass


    def start(self):
        self.feed_worker_thread.start()
        self.collect_intermediate_results_thread.start()
        self.collect_final_results_thread.start()
        self.collect_intermediate_data_broadcast_thread.start()
        self.acknowledge_intermediate_data_thread.start()

    def report_job_status(self):
        # report the status of all processes/jobs
        status = "\n".join(["%2d - %s (%5d, #%2d): %5s %5s %5s %10d %5s" % (
                job['ota_id'], job['filename'], 
                -1 if job['process'] is None else job['process'].pid,
                job['attempt'],
                job['intermediate_results_received'], 
                job['intermediate_data_sent'], 
                job['intermediate_data_ackd'], 
                (time.time()-job['time_of_ack']),
                bool(job['complete'].value)) for job in self.info])
        self.logger.debug("Process status:\n\n%s\n" % (status))
        # for job in self.info:
        #     self.logger.info("%2d - %s (%5d, #%2d): %5s %5s %5s %10d %5s" % (
        #         job['ota_id'], job['filename'], 
        #         -1 if job['process'] == None else job['process'].pid,
        #         job['attempt'],
        #         job['intermediate_results_received'], 
        #         job['intermediate_data_sent'], 
        #         job['intermediate_data_ackd'], 
        #         (time.time()-job['time_of_ack']),
        #         job['complete'],
        #     )
        # )

    def feed_workers(self):
        self.workers_started = 0
        should_be_working = []
        process_ids = []

        #print "Starting to feed workers"

        x = 0
        while (not self.quit):
            #
            # Make sure all workers that should be working are doing so
            #
            workers_alive = 0
            for job in self.info:
                if (job['complete'].value):
                    continue

                process = job['process']
                if (process is None):
                    continue

                pid = process.pid
                if (pid is None):
                    # not started yet
                    continue

                ps = psutil.Process(pid)
                _dead = ps.status() in [psutil.STATUS_ZOMBIE,
                                        psutil.STATUS_DEAD]
                _timeout = (job['time_of_ack'] > 0 and 
                     job['intermediate_data_sent'] >= 1 and
                     job['intermediate_data_ackd'] == True and
                    (time.time()-job['time_of_ack']) > 10 and
                            not job['complete'].value)
                if (_dead):
                    self.logger.warning("Found dead process: pid=%d ota-id:%d fn=%s" % (
                        pid, job['ota_id'], job['filename']))
                if (_timeout):
                    self.logger.warning("Found very slow process: pid=%d ota-id:%d fn=%s" % (
                        pid, job['ota_id'], job['filename']))
                if (_dead or _timeout):
                    self.logger.info("Restarting dead/slow process - ota-id:%d fn=%s" % (
                        job['ota_id'], job['filename']))
                    self.job_status_lock.acquire()
                    process.terminate()
                    process.join(timeout=0.01)
                    job['process'] = None
                    job['intermediate_data_sent'] = 0
                    #job['intermediate_results_received'] = False
                    job['intermediate_data_ackd'] = False
                    self.active_workers -= 1
                    self.report_job_status()
                    self.job_status_lock.release()
                    self.logger.debug("# active workers is now %d" % (self.active_workers))
                    continue

                workers_alive += 1
            # self.logger.debug("%d workers still alive" % (workers_alive))

            if (x==0):
                # only print this information once!
                fns = ["%s" % (f['filename']) for f in self.info]
                self.logger.debug("filenames being reduced:\n- %s" % ("\n- ".join(fns)))
                self.logger.debug("filename-list:\n- %s" % ("\n -".join(self.files_to_reduce)))

                for i, job in enumerate(self.info):
                    self.logger.debug("JOB %2d: ID=%d, FN: %s / %s" % (
                        i, job['ota_id'], job['filename'], job['args']['filename']))
                # self.report_job_status()

            x += 1

            #
            # Start new workers if we have CPUs available
            #
            if (self.active_workers < self.number_cpus and not self.quit):
                # self.logger.debug("we have some capacity to start new workers (%d < %d)" % (
                #     self.active_workers, self.number_cpus))

                #
                # Check all workers, and start one if we find one that's not alive
                #
                started_new_process = False
                # for fn in self.files_to_reduce: #self.info:
                for job in self.info:
                    if (job['process'] is None):
                        
                        job['attempt'] += 1

                        self.logger.debug("starting worker for ota %d, %s: attempt #%d (i.r.: %s)" % (
                            job['ota_id'], job['filename'], job['attempt'], str(job['intermediate_queue_msg'])))
                        #print "\n\n\nSTARTING:", id, job['filename'],"\n\n\n"

                        self.job_status_lock.acquire()
                        if (self.quit):
                            break

                        p = multiprocessing.Process(target=parallel_collect_reduce_ota, 
                                                    kwargs=job['args'])
                        p.daemon = True
                        job['process'] = p
                        job['process'].start()
                        self.job_status_lock.release()

                        job['intermediate_data_sent'] = 0
                        # if (job['intermediate_queue_msg'] is not None):
                        #     # this process is being started after intermediate data
                        #     # has already been sent
                        #     # --> re-queue one more intermediate to allow completion
                        #     self.logger.info("Resending intermediate results for ID %d, fn %s\nmsg=(%s)" % (
                        #         job['ota_id'], job['filename'], str(job['intermediate_queue_msg'])))
                        #     self.intermediate_queue.put(job['intermediate_queue_msg'], block=True)
                            
                        self.active_workers += 1
                        self.logger.debug("# active workers is now %d" % (self.active_workers))
                        #process_ids.append(job['process'].pid)
                        started_new_process = True
                        break
                if (started_new_process):
                    # self.broadcast_intermediate_data(None)
                    continue

            if (x%60 == 0):
                self.logger.debug("still feeding workers (%d < %d)" % (
                     self.active_workers, self.number_cpus))

                self.report_job_status()
                #print ",".join(["%d" % (p) for p in process_ids])
            time.sleep(0.1)
        self.logger.debug("Shutting down feed_workers")
        
    def collect_intermediate_results(self):
        self.logger.debug("Starting to collect intermediate results (%d)" % (len(self.info)))
        self.intermediate_results_collected = 0

        while (self.intermediate_results_collected < len(self.info) and
               not self.quit):

            try:
                #results = self.intermediate_results_queue.get_nowait()
                results = self.intermediate_results_queue.get(timeout=0.05)
            except queue.Empty:
                #time.sleep(0.05)
                continue

            self.intermediate_results.append(results)


            ota_id, data_products, shmem_id = results
            for job in self.info:
                if (job['ota_id'] == ota_id):
                    job['intermediate_results_received']
                    self.logger.debug("received some results from %d / %s!" % (
                        ota_id, job['filename']))
                    job['intermediate_results_received'] = True

            self.intermediate_results_collected += 1
            self.active_workers -= 1
            self.logger.debug("# active workers is now %d" % (self.active_workers))


        #print "***\n"*5,"All intermediate progress data received","\n***"*5
        if (self.quit):
            self.logger.debug("Quitting out of collect_intermediate_results")
        else:
            self.logger.debug("All intermediate progress data received")
            self.intermediate_results_done.release()
            self.intermediate_results_complete = True
        self.logger.debug("Shutting down collect_intermediate_results")

    def wait_for_intermediate_results(self):
        self.intermediate_results_done.acquire()
        self.intermediate_results_done.release()

        # Now we have all results
        return self.intermediate_results

    def get_intermediate_results(self):
        return self.intermediate_results

    def abort(self):
        self.quit = True
        self.qr_quit.value = True

        if (self.all_closed):
            return

        if (not self.intermediate_results_complete):
            self.intermediate_results_done.release()
        self.logger.debug("Terminating feeder")
        #self.feed_worker_thread.terminate()
        self.report_job_status()
        for job in self.info:
            try:
                p = job['process']
                p.join(timeout=0.1)
                self.logger.debug("terminating process for %s (alive? %s)" % (
                    job['filename'], p.is_alive()))
                p.terminate()
                p.join()
                self.logger.debug("done!")
            except:
                pass

        #
        # Join all threads to make sure everything is shut down
        #
        self.logger.debug("Joining collect_intermediate_results thread")
        self.collect_intermediate_results_thread.join()
        self.logger.debug("done joining!")

        self.logger.debug("Joining broadcast_intermediate_data thread")
        self.collect_intermediate_data_broadcast_thread.join()
        self.logger.debug("done joining!")

        self.logger.debug("Joining acknowledge_intermedidate_data thread")
        self.acknowledge_intermediate_data_thread.join()
        self.logger.debug("done joining!")

        self.logger.debug("Joining collect_final_results thread")
        self.collect_final_results_thread.join()
        self.logger.debug("done joining!")

        #
        # Close all queues
        #
        self.close_queues()

        self.all_closed = True

    def close_queues(self):
        self.logger.debug("Starting to close all queues")

        self.queue.close()
        self.queue.join_thread()
        self.logger.debug("Job queue closed")

        try:
            nq=0
            while (True):
                self.final_results_queue.get(block=True,timeout=0.001)
                nq += 1
        except queue.Empty:
            pass
        self.logger.debug("inserting sentinel")
        self.final_results_queue.put(multiprocessing.queues._sentinel)
        self.logger.debug("closing queue")
        self.final_results_queue.close()
        self.logger.debug("Joining thread")
        self.final_results_queue._thread.join(timeout=0.1) #join_thread()
        self.logger.debug("Final results queue closed (alive: %s / %d)" % (
            str(self.final_results_queue._thread.is_alive()), nq))

        try:
            nq=0
            while (True):
                self.intermediate_queue.get(block=True,timeout=0.001)
                nq += 1
        except queue.Empty:
            pass
        self.logger.debug("inserting sentinel")
        self.intermediate_queue.put(multiprocessing.queues._sentinel)
        self.logger.debug("closing queue")
        self.intermediate_queue.close()
        # self.intermediate_queue.join_thread()
        self.logger.debug("Joining thread")
        self.intermediate_queue._thread.join(timeout=0.1) #join_thread()
        self.logger.debug("intermediate data queue closed (alive: %s / %d)" % (
            str(self.intermediate_queue._thread.is_alive()), nq))
        #self.intermediate_queue._thread = None

        try:
            nq=0
            while (True):
                self.intermediate_results_queue.get(block=True,timeout=0.001)
                nq += 1
        except queue.Empty:
            pass
        self.logger.debug("inserting sentinel")
        self.intermediate_results_queue.put(multiprocessing.queues._sentinel)
        self.logger.debug("closing queue")
        self.intermediate_results_queue.close()
        self.logger.debug("Joining thread")
        self.intermediate_results_queue._thread.join(timeout=0.1) #join_thread()
        self.logger.debug("Intermediate results queue closed (alive: %s / %d)" % (
            str(self.intermediate_results_queue._thread.is_alive()), nq))
        #self.intermediate_results_queue._thread = None

        try:
            self.logger.debug("Starting to empty intermediate_data_ack_queue")
            nq=0
            while (True):
                self.intermediate_data_ack_queue.get(block=True,timeout=0.001)
                nq += 1
                self.logger.debug("Received message from intermediate_data_ack_queue: %d" % (nq))
        except queue.Empty:
            self.logger.debug("Received Queue.Empty exception while clearing intermediate_data_ack_queue")
            pass
        self.logger.debug("inserting sentinel")
        self.intermediate_data_ack_queue.put(multiprocessing.queues._sentinel)
        self.logger.debug("closing queue")
        self.intermediate_data_ack_queue.close()
        self.logger.debug("Joining thread")
        self.intermediate_data_ack_queue._thread.join(timeout=0.1) #join_thread()
        self.logger.debug("intermediate data received Ack queue closed (alive: %s / %d)" % (
            str(self.intermediate_data_ack_queue._thread.is_alive()), nq))
        #self.intermediate_data_ack_queue._thread = None

        self.logger.debug("all queues closed and threads terminated")
        
    def acknowledge_intermediate_data_received(self):
        while (not self.quit):
            try:
                ota_id = self.intermediate_data_ack_queue.get_nowait()
                for job in self.info:
                    if (job['ota_id'] == ota_id):
                        self.logger.debug("Got acknowledgement from OTA %d, fn %s" % (
                            ota_id, job['filename']))

                        self.job_status_lock.acquire()
                        if (job['intermediate_data_sent'] >= 1):
                            job['intermediate_data_ackd'] = True
                        else:
                            self.logger.warning("Received ACK before command")
                        self.report_job_status()
                        job['time_of_ack'] = time.time()
                        self.job_status_lock.release()
                        break
                continue
            except queue.Empty:
                time.sleep(0.1)
                pass
        self.logger.debug("Shutting down acknowledge_intermediate_data_received")

    def broadcast_intermediate_data(self):
        while (not self.quit):
            for job in self.info:
                if (not job['intermediate_results_received'] or
                    job['intermediate_data_ackd'] or
                    job['intermediate_data'] is None or
                    job['intermediate_queue_msg'] is None):
                    # nothing to do for this one
                    continue

                try:
                    self.intermediate_queue.put(job['intermediate_queue_msg'])
                except AssertionError:
                    pass

                self.logger.debug("Putting one set (# %d / %d) of intermediate data back in work queue" % (
                    self.intermed_results_sent, job['intermediate_data_sent']))
                job['intermediate_data_sent'] += 1
                self.intermed_results_sent += 1
            time.sleep(0.2)
        self.logger.debug("Shutting down broadcast_intermediate_data")

    def send_intermediate_data(self, ota_id, data):
        for job in self.info:
            if (job['ota_id'] == ota_id):
                self.job_status_lock.acquire()
                job['intermediate_data'] = data
                msg = (ota_id, data)
                job['intermediate_queue_msg'] = msg
                self.logger.debug("XXX: %d, %s" % (job['ota_id'], str(job['intermediate_data'])))
                job['intermediate_data_sent'] = 0
                job['intermediate_data_ackd'] = False
                self.job_status_lock.release()
                self.active_workers += 1
                break
        self.report_job_status()
        pass

    def collect_final_results(self):

        n_collected = 0
        while (not self.quit and
               n_collected < len(self.info)):

        #for i in range(len(self.info)):
            try:
                result = self.final_results_queue.get(timeout=0.05)
            except queue.Empty:
                continue


            # mark the dataset as complete
            ota_id, data_products, shmem_id = result
            self.logger.debug("received final results for ota-ID %d (%s)" % (
                ota_id, ",".join(["%d" % job['ota_id'] for job in self.info])))


            # find the right job
            found_job = False
            for job in self.info:
                if (job['ota_id'] == ota_id):
                    # fn = self.id2filename[ota_id]
                    self.job_status_lock.acquire()
                    job['complete'].value = True #<<- this is now done from within the worker process
                    self.logger.debug("File %s, ID %d marked as complete" % (job['filename'], ota_id))
                    self.active_workers -= 1
                    self.final_results.append(result)
                    self.logger.debug("# active workers is now %d" % (self.active_workers))
                    found_job = True
                    n_collected += 1
                    self.job_status_lock.release()
                    # time.sleep(1)
                    break

            if (not found_job):
                self.logger.error("Could not associate final results!")

        self.report_job_status()
        self.final_results_done.release()
        self.quit = True
        self.qr_quit.value = True
        self.logger.debug("Shutting down collect_final_results")

    def wait_for_workers_to_finish(self):
        self.final_results_done.acquire()
        self.final_results_done.release()
        self.report_job_status()

        # Now we are all done
        self.quit = True

    def get_final_results(self):
        return self.final_results

    def reduce_file(self, filename, id):

        for job in self.info:
            if (job['ota_id'] == id or job['filename'] == filename):
                self.logger.error("we are already working on %s" % (filename))
                return None

        job = {}

        self.id2filename[id] = filename
        self.filename2id[filename] = id

        #
        # Setup a new process for this file
        #
        #shmem = multiprocessing.RawArray(ctypes.c_float, 4096*4096)
        shmem = SharedMemory(ctypes.c_float, (4096,4096))
        self.shmem_list[id] = shmem
        job['complete'] = multiprocessing.Value('i', False)
        job['args'] = {
            'workqueue': self.queue,
            'intermediate_results_queue': self.intermediate_results_queue,
            'final_results_queue': self.final_results_queue,
            'intermediate_queue': self.intermediate_queue,
            'options': self.options,
            'shmem': shmem,
            'shmem_id': id,
            'shmem_dim': self.shmem_dims,
            'filename': filename,
            'ota_id': id,
            'intermediate_ack_queue': self.intermediate_data_ack_queue,
            'quit_signal': self.qr_quit,
            'complete': job['complete']
            # 'intermediate_results_queue_lock': self.intermediate_results_queue_lock,
        }
        job['intermediate_data'] = None
        job['intermediate_queue_msg'] = None
        job['filename'] = filename
        job['ota_id'] = id
        job['attempt'] = 0
        job['intermediate_data_sent'] = 0
        job['intermediate_results_received'] = False
        job['intermediate_data_ackd'] = False
        job['time_of_ack'] = -1

        self.logger.debug("Setting up reduction for %s (ID: %d) -> %s" % (filename, id, job['args']['filename']))
        job['process'] = None

        self.info.append(job)

        self.files_to_reduce.append(filename)
        pass

        #, directory, filebase, otas):

        # print "Reducing these files:"
        # for ota in otas:
        #     filename = "%s/%s.%02d.fits" % (directory, filebase, ota)
        #     print filename


    def free_shared_memory(self):

        self.logger.debug("Freeing up shared memory")
        # print self.shmem_list

        for shmem_id in self.shmem_list:
            self.shmem_list[shmem_id].free()

            








def collectcells(input, outputfile,
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

    if (options is None): options = set_default_options()

    logger = logging.getLogger('CollectCells')
    logger.debug("Starting --collectcells--")

    if (showsplash):
        splash = """\

    **********************************************************************
    * This is podi_collectcells                                          *
    * (c) 2012-2017: Ralf Kotulla, kotulla@wisc.edu                      *
    *                University of Wisconsin - Madison                   *
    *                WIYN Observatory, Inc                               *
    *                                                                    *
    * Please acknowledge the author when using any products generated    *
    * with this tool. For comments, questions or ideas for improvement   *
    * please send an email to kotulla@wisc.edu. Thank you!               *
    **********************************************************************

"""
        stdout_write(splash)

    # print "Received options:", options
    # podi_logging.print_stacktrace()
    # time.sleep(1)
    # return 1

    if (options['verbose']):
        stdout_write("\nThese are the options we are using:\n")
        for opt_name, opt_value in options.items():
            stdout_write("   %30s: %s\n" % (opt_name, opt_value))
        stdout_write("----- end options\n\n")

    # afw = podi_asyncfitswrite.async_fits_writer(1)

    staged_data = False
    if (options['prestage']):
        input, staged_data = prestage_data(options, input)
#        return

    if (os.path.isfile(input)):
        logger.debug("Input is a file: %s" % (input))
        # Assume this is one of the fits files in the right directory
        # In that case, extract the FILENAME header and convert it into 
        # the filebase we need to construct the filenames of all OTA fits files.
        hdulist = pyfits.open(input)

        filebase = hdulist[0].header['FILENAME'][:-8]
        hdulist.close()
        del hdulist

        # Split the input filename to extract the directory part
        directory, dummy = os.path.split(input)
        if (directory is None or directory == ''):
            directory = "."

    elif (os.path.isdir(input)):
        # As a safety precaution, if the first parameter is the directory containing 
        # the files, extract just the ID string to be used for this script
        if (input[-1] == "/"):
            input = input[:-1]

        basedir, filebase = os.path.split(input)
        directory = input
        logger.debug("Input is a directory: %s --> %s / %s" % (input, basedir, filebase))

    else:
        logger.error("Input (%s) is neither path nor file, aborting!\n" % (input))
        #
        # Try tracing back item by item where we stop loosing the path
        #
        fp = os.path.abspath(input)
        items = fp.split("/")
        for i in range(1, len(items)):
            part_path = "/".join(items[:i+1])
            logger.debug("BACKTRACKING: %-100s is a path/dir? %-5s / %-5s" % (
                part_path, os.path.isdir(part_path), os.path.isfile(part_path)))
            # print "BACKTRACKING: %-100s is a path/dir? %-5s / %-5s" % (
            #     part_path, os.path.isdir(part_path), os.path.isfile(part_path))
        unstage_data(options, staged_data, input)
        return

    if (outputfile is None):
        outputfile = "%s/%s.fits" % (directory, filebase)

    hdulist = None
    for ota in podi_focalplanelayout.FocalPlaneLayout().create_radially_sorted_ota_list():
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

    if (hdulist is None):
        stdout_write("Something is wrong here, can't find/open any of the files...")
        unstage_data(options, staged_data, input)
        return -1

    fpl = podi_focalplanelayout.FocalPlaneLayout(hdulist[0])
    logger.info("Focalplane Layout: %s" % (fpl.layout))

    obsid = hdulist[0].header['OBSID']

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
                objectname = header['OBJECT']
                strings_to_replace = [',', '(', ')', '/', '\\', '`', '"', '\'', '{', '}', '[', ']', '&',
                                      ' ', '*', '?']
                for _s in strings_to_replace:
                    objectname = objectname.replace(_s, '_')
                outputfile = outputfile[:start] + objectname  + outputfile[start+7:]
            elif (outputfile[start:start+6] == "%OBSID"):
                outputfile = outputfile[:start] + header['OBSID'] + outputfile[start+6:]
            elif (outputfile[start:start+7] == "%OBSSEQ"):
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
    obsid = hdulist[0].header['OBSID']

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
    if (not (
            (outputfile.endswith(".fits") and not options['compressed_hdu']) or
            (outputfile.endswith(".fits.fz") and options['compressed_hdu'])) ):
        # Output filenames have to end with .fits
        ext = ".fits.fz" if options['compressed_hdu'] else ".fits"
        logger.warning("no fits extension given to filename, adding %s" % (ext))
        outputfile += ext
        changed_outputfilename = True
    if (changed_outputfilename):
        logger.info("After revision, output filename now is %s" % (outputfile))

    if (os.path.isfile(outputfile) and not options['clobber']):
        logger.info("File %s already exists, skipping!" % (outputfile))
        print("#####################################################")
        print("#")
        print("# File %s already exists, skipping!" % (outputfile))
        print("#")
        print("#####################################################")
        print("\n")
        unstage_data(options, staged_data, input)
        return

    #
    # Start assembling the new list of HDUs
    #
    list_of_otas_to_collect = fpl.available_ota_coords
    if (options['central_only']):
        list_of_otas_to_collect = fpl.central_array_ota_coords
    if (not type(options['selectota']) == type(None)):
        list_of_otas_to_collect = options['selectota']

    logger.debug("List of otas to collect: %s" % (str(list_of_otas_to_collect)))

    filtername = hdulist[0].header['FILTER']
    if (filtername in fpl.blocked_out_otas):
        exclude_ota = fpl.blocked_out_otas[filtername]
        unblocked_otas = []
        for (ox, oy) in list_of_otas_to_collect:
            ota = ox * 10 + oy
            if (not ota in exclude_ota):
                unblocked_otas.append((ox,oy))
        list_of_otas_to_collect = unblocked_otas
        logger.debug("List of UN-BLOCKED OTAs: %s" % (str(list_of_otas_to_collect)))

    ota_list = [None] * (len(list_of_otas_to_collect)+1)
    # And add the primary HDU to make the fits file a valid one
    ota_list[0] = pyfits.PrimaryHDU()

    # Add version data
    record_pipeline_versioning(ota_list[0].header)

    # add user-defined additional keywords
    if (len(options['additional_fits_headers']) > 0):
        _firstkey = None
        for key, value in options['additional_fits_headers'].items():
            ota_list[0].header[key] = (value, "user-added keyword")
            _firstkey = key if _firstkey is None else _firstkey
        add_fits_header_title(ota_list[0].header, "User-added keywords", _firstkey)

    ota_list[0].header['BINNING'] = (binning, "binning factor")

    # Creates the persistency index.cat file
    if (options['persistency_dir'] is not None):
        logger.info("Gathering persistency data...")
        podi_persistency.get_list_of_saturation_tables(options['persistency_dir'])
        logger.debug("Done gathering persistency data!")

    #
    # Set up the parallel processing environment
    #
    #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
    # queue = multiprocessing.JoinableQueue()
    # return_queue = multiprocessing.Queue()
    # intermediate_queue = multiprocessing.Queue()

    processes = []

    #
    # For easier user-understanding, add special keywords for each reduction step
    #

    # worker_args = (queue, return_queue, options)
    shmem_dims = (4096, 4096)
    # kw_worker_args= {
    #     'queue': queue,
    #     'return_queue': return_queue,
    #     'intermediate_queue': intermediate_queue,
    #     'options': options,
    #     'shmem': None,
    #     'shmem_dim': shmem_dims,
    # }
    shmem_list = []
    number_extensions = 0

    list_of_otas_being_reduced = []
    ota_ids_being_reduced = []
    intermediate_results = []

    # Set up all the communication pipes to communicate data back to the process
    communication = {}
    ota_from_otaid = {}
    mp_params = {}

    logger.debug("Setting up parallel workers to collect & reduce OTAs")
    worker = reduce_collect_otas(options, number_cpus)
    for ota_id in range(len(list_of_otas_to_collect)):
        ota_c_x, ota_c_y = list_of_otas_to_collect[ota_id]        
        ota = ota_c_x * 10 + ota_c_y

        ota_from_otaid[ota_id+1] = ota

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

        # # Setup some way of communicating with the workers
        # msg_pipe = multiprocessing.Pipe(duplex=False)
        # pipe_recv, pipe_send = msg_pipe
        # # Make sure to pass the receiver pipe to the worker so it can listen 
        # # for instructions. Do some python wrapping to be able to transport a 
        # # pipe object though a pipe/queue
        # wrapped_pipe = multiprocessing.reduction.reduce_connection(pipe_recv)
        
        # queue.put( (filename, ota_id+1, wrapped_pipe) )
        # mp_params[ota] = (filename, ota_id+1, wrapped_pipe)
        # # del wrapped_pipe
        ota_ids_being_reduced.append(ota_id+1)

        # # save some data we need later on for the intermediate results
        # intres = {'ota-id': ota_id+1,
        #           'sent': False,
        #           'queued': True,
        #           'pipe-recv': pipe_recv,
        #           'pipe-send': pipe_send,
        # }
        # intermediate_results.append(intres)

        logger.debug("Queuing up reduction for OTA %s" % (filename))
        worker.reduce_file(filename, ota_id+1)

        
    logger.debug("list_of_otas_being_reduced=\n%s" % (str(list_of_otas_being_reduced)))

    logger.info("Performing instrumental detrending")
    podi_logging.ppa_update_progress(0, "Starting work")
    worker.start()

    logger.debug("Waiting for de-trending to proceed before continuing!")
    worker.wait_for_intermediate_results()

    # Create all processes to handle the actual reduction and combination
    #print "Creating",number_cpus,"worker processes"
    # if ('profile' in options or number_cpus == 0):
    #     # 
    #     # If profiling is activated, run one only one processor and in non-multiprocessing mode
    #     #
    #     # Tell all workers to shut down when no more data is left to work on
    #     #for i in range(len(processes)):
    #     if (verbose): stdout_write("Sending quit command!\n")
    #     queue.put((True,None,None))

    #     while (True):
    #         print "Doing work single-processed"
    #         cmd_quit, filename, ota_id = queue.get()
    #         if (cmd_quit):
    #             queue.task_done()
    #             break

    #         # Do the work
    #         data_products = collect_reduce_ota(filename, options=options)

    #         # Add the results to the return_queue so the master process can assemble the result file
    #         # print "Adding results for OTA",ota_id,"to return queue"
    #         # return_queue.put( (hdu, ota_id, wcsfix_data) )
    #         return_queue.put( (ota_id, data_products) )
    #         queue.task_done()
    # else:
    #     for i in range(number_cpus):
    #         logger.debug("Starting a new process...")
    #         #p = multiprocessing.Process(target=parallel_collect_reduce_ota, args=worker_args)

    #         # Allocate shared memory for data return
    #         new_shmem = multiprocessing.RawArray(ctypes.c_float, 4096*4096)
    #         kw_worker_args['shmem'] = new_shmem
    #         kw_worker_args['shmem_id'] = len(shmem_list)
    #         shmem_list.append(new_shmem)

    #         p = multiprocessing.Process(target=parallel_collect_reduce_ota, kwargs=kw_worker_args)
    #         p.start()
    #         processes.append(p)
    #         if (not process_tracker == None):
    #             if (verbose): print "Adding current slave process to process-tracker...", i
    #             process_tracker.put(p.pid)
    #             if (verbose): print "done adding to process-tracker"

    #         # Tell all workers to shut down when no more data is left to work on
    #         #queue.put(None)
    #         time.sleep(0.01)



    
        

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
    expmeas = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header else hdulist[0].header['EXPTIME']
    date_obs = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], date_format)
    date_mid = date_obs + datetime.timedelta(seconds=(0.5*expmeas))
    date_end = date_obs + datetime.timedelta(seconds=expmeas)
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
    ota_list[0].header['MJD-MID'] = (mjd + (0.5 * expmeas)/86400., 
                                     "MJD at exposure mid-point")
    ota_list[0].header['MJD-END'] = (mjd + (expmeas/86400.), 
                                     "MJD at exposure end-point")
    add_fits_header_title(ota_list[0].header, "Additional time stamps", 'DATE-MID')

    #
    # Gather some information from the guider, compute statistics,
    # and create the guide star diagnostic plots
    #
    guide_files = podi_guidestars.get_guidephotom_filelist(directory, obsid)
    logger.debug("Found the following guide photometry files:\n-- %s" % (
        "\n-- ".join(guide_files)))
    guide_plot_filename = outputfile[:-5]+".guide.png"
    guide_plot_title = "%(OBSID)s: %(OBJECT)s (%(FILTER)s, %(EXPTIME)d s)" % hdulist[0].header
    guidestats = podi_guidestars.draw_guidestarplot(
        guide_files, 
        title=guide_plot_title,
        plot_filename=guide_plot_filename)
    #
    # Now add some of the guide information to the FITS header
    #
    ota_list[0].header['N_GUIDES'] = (len(guidestats), 
                                      'number of guide stars')
    total_flux_max, total_flux_min, total_guide_samples = 0, 0, 0
    for idx, starfile in enumerate(guidestats):
        star = guidestats[starfile]
        total_flux_max += star['flux_max']
        total_flux_min += star['flux_min']
        total_guide_samples += star['n_guide_samples']
        ota_list[0].header['SGMA1F_%d' % (idx+1)] = star['flux_1sigma']
        ota_list[0].header['SGMA3F_%d' % (idx+1)] = star['flux_3sigma']

    # print "phot:", total_flux_max, total_flux_min
    photometricity = (total_flux_min / total_flux_max) \
                     if ((total_flux_max > 0) and (total_flux_min > 0)) else -1. 
    ota_list[0].header['PHOTQUAL'] = (photometricity,
                                      'degree of photometric stability (1=photometric)')
    ota_list[0].header['GUIDSMPL'] = (0 if total_guide_samples <= 0 
                                      else total_guide_samples/len(guidestats),
                                      'avg number of guide frames per star')
    add_fits_header_title(ota_list[0].header, "Guider statistics", "N_GUIDES")

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
    ota_headers = {}
    intermediate_results_sent = {}
    pupilghost_scaling = None
    ota_missing_empty = []

    global_reduction_log = ReductionLog()

    # otas_checked_in = []
    # otas_not_checked_in = [x*10+y for (x,y) in list_of_otas_being_reduced]
    # timeouts_left = 3
    # while (len(otas_checked_in) < len(list_of_otas_being_reduced)):
    #     need_another_process = False
    #     received_results = False
    #     try:
    #         ota_id, data_products, shmem_id = return_queue.get(timeout=sitesetup.per_ota_timeout)
    #         received_results = True
    #     except Queue.Empty:
    #         if (len(processes) < len(list_of_otas_being_reduced)):
    #             #
    #             # This means we haven't sent off all OTAs for reduction
    #             #
    #             pass
    #             # Just wait a little longer
    #             continue
    #         else:
    #             #
    #             # All OTAs have an assigned process, so we are 
    #             # only waiting for them to finish
    #             #

    #             timeouts_left -= 1
    #             logger.warning("Received timeout (%d s), %d left before re-submission!" % (sitesetup.per_ota_timeout, timeouts_left))
    #             if (timeouts_left <= 0):
    #                 logger.error("No timeouts remaining, program stalled, adding another process to compensate!")

    #                 #
    #                 # Add another instance
    #                 #
    #                 need_another_process = True

    #                 #
    #                 # pick one of the OTAs that's missing/not yet returned
    #                 #
    #                 _ota_id = otas_not_checked_in[0]
    #                 logger.error("OTAs yet to report intermed. results back (%2d/%2d): %s" % (
    #                     len(otas_not_checked_in), len(list_of_otas_being_reduced), 
    #                     ",".join(["%02d" % x for x in otas_not_checked_in])))
    #                 logger.info("Re-queing work for OTA-ID %d" % (_ota_id))

    #                 q_params = mp_params[_ota_id]
    #                 queue.put( q_params )
    #                 timeouts_left = 3

    #             # unstage_data(options, staged_data, input)
    #             # return None
    #             # i -= 1
    #         pass
    #     except (KeyboardInterrupt, SystemExit):
    #         while (not return_queue.empty()):
    #             return_queue.get()
    #         raise
    #         unstage_data(options, staged_data, input)
    #         return

    #     #
    #     # We received one entry. Check if we need to start another process
    #     # If all processes are running we should have as many processes as
    #     # we have OTAs to be reduced
    #     #
    #     # Doing this in here ensures we only have a limited number of processes 
    #     # doing active work, hence keeping the machine from overloading.
    #     # 
    #     if (len(processes) < len(list_of_otas_being_reduced) or need_another_process):
    #         # We don't have enough processes yet, start another one

    #         # Allocate shared memory for data return
    #         new_shmem = multiprocessing.RawArray(ctypes.c_float, 4096*4096)
    #         kw_worker_args['shmem'] = new_shmem
    #         kw_worker_args['shmem_id'] = len(shmem_list)
    #         shmem_list.append(new_shmem)

    #         p = multiprocessing.Process(target=parallel_collect_reduce_ota, kwargs=kw_worker_args)
    #         p.start()
    #         processes.append(p)
    #         logger.debug("Starting another process for another OTA")
    #         if (need_another_process):
    #             logger.info("Adding another process to make up for the dead one!")
    #         if (not process_tracker == None):
    #             if (verbose): print "Adding current slave process to process-tracker...", i
    #             process_tracker.put(p.pid)
    #             if (verbose): print "done adding to process-tracker"
    #         # Also send another quit command for this process
    #         #queue.put(None)

    #     if (not received_results):
    #         continue

    #     if (timeouts_left < 3):
    #         logger.info("Received valid data, resetting timeout counter!")
    #         timeouts_left = 3

    for i, intermed_results in enumerate(worker.get_intermediate_results()):

        if (intermed_results is None):
            logger.error("Illegal results")
            
        ota_id, data_products, shmem_id = intermed_results
        # print(data_products['psf'])

        logger.debug("Received intermediate results from OTA-ID %02d" % (ota_id))
        podi_logging.ppa_update_progress(int(50.*(i+1)/len(list_of_otas_being_reduced)), "Reducing")

        # _ota = ota_from_otaid[ota_id]
        # otas_checked_in.append(_ota)
        # idx = otas_not_checked_in.index(_ota)
        # try:
        #     del otas_not_checked_in[idx]
        # except:
        #     podi_logging.log_exception()
        # logger.debug("Received intermed. results from OTAs (%2d/%2d): %s" % (
        #     len(otas_checked_in), len(list_of_otas_being_reduced), 
        #     ",".join(["%02d" % x for x in otas_checked_in])))
        # logger.debug("OTAs yet to report intermed. results back (%2d/%2d): %s" % (
        #     len(otas_not_checked_in), len(list_of_otas_being_reduced), 
        #     ",".join(["%02d" % x for x in otas_not_checked_in])))


        # # Mark this ota as not fully complete. This is important later on when
        # # we report intermediate results back for completion
        # intermediate_results_sent[ota_id] = False

        header = data_products['header']
        
        if (header is None):
            logger.warning("OTA %d reporting empty header, ignoring" %(ota_id))
            ota_missing_empty.append(ota_id)
            continue

        extname = header['EXTNAME']
        ota_headers[extname] = header

        sky_samples[header['EXTNAME']] = data_products['sky-samples']

        if (not type(data_products['fringe_scaling']) == type(None)):
            fringe_scaling = data_products['fringe_scaling'] if type(fringe_scaling) == type(None) else \
                numpy.append(fringe_scaling, data_products['fringe_scaling'], axis=0)
            logger.debug("XXX fringe scaling:\n%s" % (str(fringe_scaling)))

        # Add something about pupilghost scaling XXXXXXXXXXXXX
        if (not type(data_products['pupilghost-scaling']) == type(None)):
            pupilghost_scaling = data_products['pupilghost-scaling'] if (type(pupilghost_scaling) == type(None)) \
                else numpy.append(pupilghost_scaling, data_products['pupilghost-scaling'], axis=0)
            
        # print "\n\n\n\n",pupilghost_scaling,"\n\n\n\n"

    logger.info("Received all intermediate data")

    ############################################################################
    #
    # To do now: 
    # 1) Compute the global values for pupil-ghost correction and fringe scaling
    # 2) Send the results back to all worker threads so they can finish up their 
    #    work
    #
    ############################################################################
        
    # 
    # Now combine all sky-samples to compute the global background level.
    # Take care to exclude OTAs marked as guide-OTAs and those not covered 
    # by the narrow-band filters.
    # 
    logger.debug("Combining sky-samples from all OTAs into global sky value")
    sky_samples_global = None #numpy.empty(0)
    # logger.info("Filtername: %s (%s)" % (
    #     ota_list[0].header['FILTER'] if 'FILTER' in ota_list[0].header else "<NO FILTER KEYWORD>",
    #     ota_list[0].header['FILTERID'] if 'FILTERID' in ota_list[0].header else "<NO FILTERID KEYWORD>"
    #     ))
    valid_ext = fpl.otas_for_photometry[get_valid_filter_name(ota_list[0].header)]
    sky_global_median = -1.
    for ext in sky_samples:
        # print ext, valid_ext, int(ext[3:5])
        ota_number = int(ext[3:5])
        ota_name = "OTA%02d.SCI" % (ota_number)
        not_a_guiding_ota = False

        if (not ota_name in ota_headers):
            logger.debug("Found OTA that I don't understand: %s" % (ota_name))
            # We don't know about this OTA, skip it
            continue

        if (sky_samples[ext] is None):
            logger.debug("Found empty sky-sample list: %s" % (ota_name))
            continue

        sky_plus_ota = numpy.empty((sky_samples[ext].shape[0], sky_samples[ext].shape[1]+1))
        sky_plus_ota[:, :sky_samples[ext].shape[1]] = sky_samples[ext][:,:]
        sky_plus_ota[:, -1] = ota_number

        not_a_guiding_ota = ota_headers[ota_name]['CELLMODE'].find("V") < 0
        if (ota_number in valid_ext and not_a_guiding_ota and sky_samples[ext] is not None):
            if (sky_samples_global is None):
                sky_samples_global = sky_plus_ota 
            else:
                sky_samples_global = numpy.append(sky_samples_global, sky_plus_ota, axis=0)
        
        if (sky_samples_global is None):
            continue
        logger.debug("Entire list of sky-samples now contains % 4d entries" % (sky_samples_global.shape[0]))

    logger.debug("Done with collecting sky sample data")
    # numpy.savetxt("skyglobal", sky_samples_global)

    try:
        sky_samples_clipped = three_sigma_clip(sky_samples_global[:,4], nsigma=3)
        sky_median_clipped = numpy.median(sky_samples_clipped)
        # print "clipped sky median", numpy.median(sky_samples_clipped)
        # print "minimum sky value", numpy.min(sky_samples_global[:,4])
    except:
        sky_median_clipped = invalid_sky_level_value
        
    try:
            
        sky_global_median = numpy.median(sky_samples_global[:,4])
        sky_global_min = numpy.min(sky_samples_global[:,4])
        sky_global_std = numpy.std(sky_samples_global[:,4])
        number_sky_samples = sky_samples_global.shape[0]

        percentiles = list(sigma_to_percentile(numpy.array([-1,1,-2,2])))
        sky_sigmas = scipy.stats.scoreatpercentile(sky_samples_global[:,4], percentiles)
        sky_bottom5percent = scipy.stats.scoreatpercentile(sky_samples_global[:,4], 5)
        sky_bottom5boxes = numpy.average(numpy.sort(sky_samples_global[:,4])[:5])
        # print sky_1sigmas
    except:
        logger.warning("Problem determining the global sky level")
        sky_global_median = invalid_sky_level_value
        sky_global_min = invalid_sky_level_value
        sky_global_std = invalid_sky_level_value
        number_sky_samples = 0
        sky_sigmas = [invalid_sky_level_value, invalid_sky_level_value, invalid_sky_level_value, invalid_sky_level_value]
        sky_bottom5percent = invalid_sky_level_value
        sky_bottom5boxes = invalid_sky_level_value

    logger.debug("Found global median sky-value = %.1f" % (sky_global_median))
    ota_list[0].header["SKYLEVEL"] = (sky_global_median, "median global sky level")
    ota_list[0].header["SKYBG"] = (sky_global_median, "median global sky background")
    ota_list[0].header['SKYBGCLP'] = (sky_median_clipped, "3-sigma clipped median sky level")
    ota_list[0].header['SKYBGMIN'] = (sky_global_min, "miminum sky level")
    ota_list[0].header['SKYBGSTD'] = (sky_global_std, "std.dev. of sky level")
    ota_list[0].header['SKYSMPLS'] = (number_sky_samples, "number of sky level samples")
    ota_list[0].header['SKYL1SIG'] = (sky_sigmas[0], "sky level lower 1st quantile")
    ota_list[0].header['SKYU1SIG'] = (sky_sigmas[1], "sky level upper 1st quantile")
    ota_list[0].header['SKYL2SIG'] = (sky_sigmas[2], "sky level lower 2nd quantile")
    ota_list[0].header['SKYU2SIG'] = (sky_sigmas[3], "sky level upper 2nd quantile")
    ota_list[0].header['SKY_LO5P'] = (sky_bottom5percent, "sky level, 5 percent")
    ota_list[0].header['SKY_LO5S'] = (sky_bottom5boxes, "sky level, lowest 5 samples")
    add_fits_header_title(ota_list[0].header, "Derived global data", 'SKYLEVEL')
    logger.debug("All sky-related FITS header entries written")

    ota_list[0].header['CRJ'] = (options['crj'] > 0, "cosmic ray removal done")
    ota_list[0].header['CRJ_ITER'] = (options['crj'], "cosmic ray removal iterations")
    ota_list[0].header['CRJ_SIGC'] = (options['crj_sigclip'], "CR sigma clipping threshold")
    ota_list[0].header['CRJ_SIGF'] = (options['crj_sigfrac'], "CR threshold for neighboring pixels")
    ota_list[0].header['CRJ_OBJL'] = (options['crj_objlim'], "CR contrast to underlying object")
    ota_list[0].header['CRJ_SATU'] = (options['crj_saturation'], "CR saturation limit")
    ota_list[0].header['CRJ_METH'] = (options['crj_method'], "CR implementation")
    add_fits_header_title(ota_list[0].header, "Cosmic-ray rejection setting", 'CRJ')
        
    #
    # Compute the global fringe scaling 
    # 
    fringe_scaling_median = fringe_scaling_std = 0
    if (options['fringe_dir'] is not None and fringe_scaling is not None): #.shape[0] > 0):
        logger.debug("Computing average fringe scaling amplitude")
        # Determine the global scaling factor
        good_scalings = three_sigma_clip(fringe_scaling[:,6], [0, 1e9])
        fringe_scaling_median = numpy.median(good_scalings)
        fringe_scaling_std    = numpy.std(good_scalings)
    else:
        fringe_scaling_median, fringe_scaling_std = -1., -1.
        
    # and log the values in the primary header
    ota_list[0].header["FRNG_SCL"] = (fringe_scaling_median, "median fringe scaling")
    ota_list[0].header["FRNG_STD"] = (fringe_scaling_std, "fringe scaling uncertainty")
    logger.debug("Fringe scaling: %.3f +/- %.3f" % (fringe_scaling_median, fringe_scaling_std))
    add_fits_header_title(ota_list[0].header, "Metadata combined from individual OTAs", 'FRNG_SCL')

    #
    # Compute the global median pupil-ghost contribution
    #
    pupilghost_scaling_median = pupilghost_scaling_std = 0.
    if (options['pupilghost_dir'] != None and not type(pupilghost_scaling) == type(None)):
        logger.debug("Computing global pupilghost scaling")

        pupilghost_scaling_median, pg_background = podi_matchpupilghost.iterate_reject_scaling_factors(
            pupilghost_scaling, iterations=3,
            significant_only=False)

        logger.info("Found PG scaling: %.3f (BG: %.2f)" % (pupilghost_scaling_median, pg_background))

    else:
        pupilghost_scaling_median = -1.
    ota_list[0].header["PUPLGFAC"] = (pupilghost_scaling_median, "pupilghost scaling")
    logger.debug("Pupilghost scaling: %.3f +/- %.3f" % (pupilghost_scaling_median, pupilghost_scaling_std))

    ############################################################################
    #
    # Send back the intermediate global numbers and finish work
    #
    # Again make sure we do not overload the machine with too many concurrent,
    # active processes.
    # Combine this with reading all the remaining values returned from the 
    # worker processes
    # 
    ############################################################################

    logger.info("Computed all intermediate data parameters")
    intermed_results = {
        "pupilghost-scaling-median": pupilghost_scaling_median,
        "pupilghost-scaling-std": pupilghost_scaling_std,
        "fringe-scaling-median": fringe_scaling_median,
        "fringe-scaling-std": fringe_scaling_std,
    }
    
    # Send off the initial bunch of results to the worker threads
    # logger.debug("Intermediate results:\n%s" % (str(intermediate_results_sent)))
    logger.debug("Received %d intermediate results!" % (len(intermediate_results)))

    for ota_id in range(len(list_of_otas_to_collect)):
        worker.send_intermediate_data(ota_id+1, intermed_results)


    logger.debug("waiting for workers to finish")
    worker.wait_for_workers_to_finish()
    logger.info("All OTAs have been de-trended")

    # podi_logging.print_stacktrace()

    logger.debug("Getting final results")
    results = worker.get_final_results()
    logger.debug("Received %d final results" % (len(results)))

    logger.debug("finishing up processing")
    worker.abort()

    # return

    # n_intermed_results_sent = 0
    # otas_to_be_finalized = []
    # finalization_message = {}
    # for i in range(number_cpus):
        
    #     #if (i < len(intermediate_results)):
    #     for j in range(len(intermediate_results)):

    #         if (intermediate_results[j]['sent']):
    #             continue

    #         target_ota_id = intermediate_results[j]['ota-id']
    #         if (target_ota_id in ota_missing_empty):
    #             logger.debug("OTA %d is listed as missing" % (target_ota_id))
    #             continue

    #         pipe_send = intermediate_results[j]['pipe-send']

    #         # Sent the intermediate results
    #         logger.debug("Sending finalization data back to ota-id %02d" % (target_ota_id))
    #         #pipe_send.send(intermed_results)
    #         intermediate_queue.put((target_ota_id, intermed_results))
    #         otas_to_be_finalized.append(target_ota_id)
    #         finalization_message[target_ota_id] = (target_ota_id, intermed_results)

    #         intermediate_results[j]['sent'] = True
    #         n_intermed_results_sent += 1
    #         break

    #     else:
    #         break
    # logger.debug("Sent off %d intermediate results back to workers!" % (n_intermed_results_sent))


    #
    # Finish off work
    #
    # logger.debug("Original OTA list: %s" % (" - ".join(["%d%d"% (x,y) for (x,y) in list_of_otas_being_reduced])))
    # logger.debug("Known faulty/missing: %s" % (" - ".join(ota_missing_empty)))
    # print "---------"
    # print "expected:",list_of_otas_being_reduced
    # print "missing:", ota_missing_empty
    # print "missing_x",
    # for x in ota_missing_empty:
    #     try:
    #         print list_of_otas_being_reduced[x],
    #     except:
    #         print "XXX",
    #         pass
    # print
    # # print "missing_x: ", [list_of_otas_being_reduced[x] for x in ota_missing_empty]
    # print "---------"

    # recv_start = time.time()
    # timeouts_left = 3
    # n_expected_results = len(list_of_otas_being_reduced) - len(ota_missing_empty)
    # otas_checked_in_final = []
    # #    for i in range(len(list_of_otas_being_reduced) - len(ota_missing_empty)):
    # retry_count = {}
    # while (len(otas_checked_in_final) < len(otas_to_be_finalized)):
    #     try:
    #         ota_id, data_products, shmem_id = return_queue.get(timeout=sitesetup.per_ota_timeout)
    #         received_results = True
    #     except Queue.Empty:
    #         timeouts_left -= 1
    #         if (timeouts_left <= 0):
    #             # We have a timeout - meaning that one of the worker process
    #             # either didn't get the final message, or has died since
    #             #
    #             # find an OTA that has'nt reported back yet
    #             for _ota in otas_to_be_finalized:
    #                 if (not _ota in otas_checked_in_final):
    #                     # this OTA hasn't report back yet
    #                     intermediate_queue.put(finalization_message[_ota])
    #                     if (not _ota in retry_count):
    #                         retry_count[_ota] = 0
    #                     else:
    #                         retry_count[_ota] +=1 
    #                     timeouts_left = 3
    #                     if (retry_count[_ota] > 2):
    #                         break
    #         continue
    #     except (KeyboardInterrupt, SystemExit):
    #         while (not return_queue.empty()):
    #             return_queue.get()
    #         raise
    #         unstage_data(options, staged_data, input)
    #         return

    #     if (received_results):
    #         otas_checked_in_final.append(ota_id)

    #     logger.info("OTAS final: %s" % (",".join(["%02d" % (i) for i in otas_checked_in_final])))
    #     logger.info("OTAS req'd: %s" % (",".join(["%02d" % (i) for i in otas_to_be_finalized])))
    #         #otas_to_be_finalized.append(target_ota_id)

    #     timeouts_left = 3
    #     # We received a final answer, so if necessary send off another intermediate results
    #     try:
    #         logger.debug("received final answer from OTA-ID %02d [c=%d, FP=%s], expecting %d [%d] more" % (
    #             ota_id, i, 
    #             str(list_of_otas_to_collect[ota_id-1]), 
    #             n_expected_results-(i+1), n_expected_results))
    #     except:
    #         podi_logging.log_exception()
    #         print "Error with logger.info statement"
    #         pass

    #     for j in range(len(intermediate_results)):
    #         if (intermediate_results[j]['sent']):
    #             continue
    #         target_ota_id = intermediate_results[j]['ota-id']
    #         if (target_ota_id in ota_missing_empty):
    #             #logger.info("OTA %d is listed as missing" % (target_ota_id))
    #             continue
    #         pipe_send = intermediate_results[j]['pipe-send']
    #         # Sent the intermediate results
    #         logger.debug("Sending finalization data back to ota-id %02d" % (target_ota_id))

    #         intermediate_queue.put((target_ota_id, intermed_results))
    #         otas_to_be_finalized.append(target_ota_id)
    #         #pipe_send.send(intermed_results)

    #         intermediate_results[j]['sent'] = True
    #         n_intermed_results_sent += 1
    #         break

    #     # for j in range(len(intermediate_results)):
    #     #     if (not intermediate_results[j]['sent']):
    #     #         target_ota_id = intermediate_results[j]['ota-id']
    #     #         pipe_send = intermediate_results[j]['pipe-send']
    #     #         pipe_send.send(intermed_results)
    #     #         intermediate_results[j]['sent'] = True
    #     #         logger.info("sending new instructions:\nIntermed results sent:\n%s" % (str(intermediate_results_sent)))
    #     #         break

    #     for j in range(len(intermediate_results)): 
    #         if (intermediate_results[j]['ota-id'] == ota_id):
    #             # Close the communication pipes
    #             logger.debug("Closing intermediate results pipe for ota %d" % (ota_id))
    #             intermediate_results[j]['pipe-send'].close()
    #             intermediate_results[j]['pipe-recv'].close()

    worker.abort()

    psf_quality_data = {}

    for final_result in worker.get_final_results():

        ota_id, data_products, shmem_id = final_result
        logger.debug("Working on final results from %d" % (ota_id))

        hdu = data_products['hdu']
        if (hdu is None):
            logger.warning("Empty HDU received!")
            continue

        #
        # Unpack shared memory buffer using the information from the fits header
        #
        wx, wy = hdu.data[0], hdu.data[1]
        logger.debug("Unpacking shared memory: (ID: %d, %d x %d)" % (shmem_id, wx,wy))
        shmem_buffer = worker.shmem_list[ota_id]
        # logger.debug("shmem: %s" % (str(shmem_buffer)))
        shmem_image = shmem_buffer.to_ndarray() #).reshape(shmem_dims)
        # shmem_image = shmem_as_ndarray(shmem_buffer).reshape(shmem_dims)
        hdu.data = numpy.copy(shmem_image[:wy,:wx])

        ota_list[ota_id] = hdu

        extname = hdu.header['EXTNAME']

        wcsfix_data = data_products['wcsdata']

        if ('reduction_files_used' in data_products):
            files_this_frame = data_products['reduction_files_used']
            # print "\n\n\n\n\nfiles_this_frame=\n\n",files_this_frame,"\n\n\n"
            podi_associations.collect_reduction_files_used(master_reduction_files_used, files_this_frame)

        global_gain_sum += (hdu.header['GAIN'] * hdu.header['NGAIN'])
        global_gain_count += hdu.header['NGAIN']

        if (data_products['sourcecat'] is not None):
            global_source_cat = data_products['sourcecat'] if (global_source_cat is None) \
                else numpy.append(global_source_cat, data_products['sourcecat'], axis=0)

        if ('tech-header' in data_products):
            all_tech_headers.append(data_products['tech-header'])
            # print "techdata for ota",ota_id,"\n",data_products['tech-header']


        #
        # Collect all information for the global reduction log
        #
        ota_reduction_log = data_products['reduction_log']
        global_reduction_log.combine(ota_reduction_log)

        #
        # Collect the data for the PSF quality plot
        #
        psf = data_products['psf']
        if (psf is not None):
            logger.debug("ADDING PSF for OTA %s" % (psf.detector))
            psf_quality_data[psf.detector] = psf

    worker.free_shared_memory()

    # recv_end = time.time()
    # logger.debug("RECEIVING: %f" % ((recv_end - recv_start)))
    
    podi_logging.ppa_update_progress(60, "Detrending done, starting calibration")

    #
    # Update the global gain variables
    #
    ota_list[0].header['GAIN'] = (-1, "global average gain [e-/ADU]")
    ota_list[0].header['NGAIN'] = (0, "number of cells contribution to gain")

    # XXX ADD HERE
    # COPY GAIN METHOD FROM ANY CHILD HEADER TO PRIMARY HEADER

    if (global_gain_count > 0):
        ota_list[0].header['GAIN'] = global_gain_sum / global_gain_count
        ota_list[0].header['NGAIN'] = global_gain_count

    # Compute the noise of the sky-level based on gain and readnoise XXXXXXX
    ota_list[0].header["SKYNOISE"] = (invalid_sky_level_value, "noise level of sky background")
    if (sky_global_median > 0 and global_gain_count > 0 and ota_list[0].header['GAIN'] > 0):
        ota_list[0].header["SKYNOISE"] = math.sqrt( 8.**2 + sky_global_median*ota_list[0].header['GAIN'])

    logger.debug("all data received from worker processes!")
    logger.info("Starting post-processing")
    additional_reduction_files = {}

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
    first_inherited_header = None
    for extension in range(1, len(ota_list)):
        ota = ota_list[extension]

        if (ota is None):
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
                first_inherited_header = keyword if first_inherited_header is None \
                                         else first_inherited_header

            # By now the value should exist in the primary header, 
            # so delete it from each of the extensions
            del ota.header[header]
                
        for header in headers_to_delete_from_otas:
            # As above, make sure header exists
            if (not header in ota.header):
                continue
            del ota.header[header]

    add_fits_header_title(ota_list[0].header, "Exposure- and instrument-specific data", first_inherited_header)

    #
    # Correct all pixel coordinates by the softbin factor to keep things consistent
    #
    if (not options['softbin'] == 0):
        # check if its a multiple of 2
        sb = options['softbin']
        if (sb in [2,4,8]):
            # We did some software binning
            
            #
            # Correct all source positions given in pixels
            #
            if (global_source_cat is not None):
                headers = ['x', 'y', 'fwhm_image']
                for h in headers:
                    global_source_cat[:,SXcolumn[h]] /= sb

            # Correct the sky sampling positions and box sizes
                       
    # print(psf_quality_data)
    if (options['create_qaplots'] and
            (psf_quality_data is not None) and
            (global_source_cat is not None)):
        plotfilename = create_qa_filename(outputfile, "psf", options)
        try:
            psf_quality.make_psf_plot(ota_listing=psf_quality_data,
                                  title="%(OBSID)s (%(OBJECT)s - %(FILTER)s - %(EXPTIME)ds)" % ota_list[0].header,
                                  output_filename=plotfilename,
                                  plotformat=options['plotformat'])
        except:
            pass

    #
    # Fix the WCS if requested
    # New WCS matching using CCMatch
    #
    ota_list[0].header['WCSFIXED'] = (False, "Was WCS calibration performed")
    ota_list[0].header['WCSCAL'] = (False, "Was astrometric solution found")
    ota_list[0].header['WCSXRMS'] = (-1., "RMS in x-dir of astrometric solution")
    ota_list[0].header['WCSYRMS'] = (-1., "RMS in y-dir of astrometric solution")
    ota_list[0].header['ASTRMCAT'] = ("", "astrometric reference catalog")
    add_fits_header_title(ota_list[0].header, "World Coordinate System", 'WCSFIXED')

    enough_stars_for_fixwcs = not type(global_source_cat) == type(None) \
                              and global_source_cat.shape[0]>3
    if (options['fixwcs']):
        if (enough_stars_for_fixwcs):
            logger.debug("Found enough stars for astrometric calibration")
            global_reduction_log.attempt('wcscal')
        else:
            logger.info("Couldn't find enough stars for astrometric calibration!")
            global_reduction_log.fail('wcscal')
    else:
        global_reduction_log.not_selected('wcscal')

    logger.debug("Next up: fixwcs")
    if (options['fixwcs'] 
        and enough_stars_for_fixwcs):

        logger.info("Performing astrometric calibration")
        # The entire source catalog is located in: --> global_source_cat
        # Current HDUList is in: --> ota_list

        import dev_ccmatch
        if (dev_ccmatch.create_debug_files): numpy.savetxt("debug.wcs_raw", global_source_cat)
        #
        # Create a OTA coordinate grid to speed up matching sources from ODI to sources 
        # in reference catalog
        #
        ccmatched = dev_ccmatch.ccmatch(source_catalog=global_source_cat,
                                        reference_catalog=None, # meaning ccmtch will obtain it
                                        input_hdu=ota_list, 
                                        mode=sitesetup.fixwcs_mode,
                                        max_pointing_error=sitesetup.max_pointing_error,
                                        max_rotator_error=sitesetup.max_rotator_error,
                                        min_contrast=sitesetup.min_wcs_quality,
                                        catalog_order=sitesetup.wcscalib_order)

        ota_list[0].header['WCSFIXED'] = ccmatched['success']

        if (ccmatched['success']):
            # Use the fixed HDUList
            ota_list = ccmatched['hdulist']

        ota_list[0].header['ASTRMCAT'] = ccmatched['astrmcat'] #"2MASS"
        ota_list[0].header['WCSMXPOS'] = (ccmatched['max_pointing_error_searched'],
                                          "maximum pointing offset searched for success")
        ota_list[0].header['WCSEXPOS'] = (ccmatched['max_pointing_error'],
                                          "maximum pointing offset allowed by config")
        ota_list[0].header['WCSMXROT'] = (str(sitesetup.max_rotator_error).replace(' ',''),
                                          "maximum pointing offset compensated")
        ota_list[0].header['WCSPLIST'] = (str(sitesetup.max_pointing_error).replace(' ',''),
                                          "maximum pointing error allowed")
        ota_list[0].header['WCSQUALS'] = ("["+','.join(["%.3f" % a for a in ccmatched['contrasts']])+"]",
                                          "WCS quality")
        ota_list[0].header['WCSMINQ']  = (sitesetup.min_wcs_quality,
                                          "Minimum WCS quality for successful calibration")
        ota_list[0].header['WCSCAL'] = ccmatched['valid_wcs_solution']

        master_reduction_files_used = podi_associations.collect_reduction_files_used(
            master_reduction_files_used,
            {'wcs-reference': ccmatched['catalog_filenames']})

        #
        # Now add some headers to visualize the WCS quality (i.e. number of 
        # matched sources for each OTA)
        #
        ota_list[0].header['WCSDETAL'] = "     0     1     2     3     4     5     6     7"
        matched_cat = ccmatched['matched_src+2mass']
        for otay in range(8)[::-1]:
            wcs_qual_string = ""
            for otax in range(8):

                _ota = otax * 10 + otay

                try:
                    n_src_odi = numpy.sum((global_source_cat[:,SXcolumn['ota']] == _ota))
                except:
                    n_src_odi = 0

                try:
                    n_matched = numpy.sum((matched_cat[:, SXcolumn['ota']] == _ota))
                except:
                    n_matched = 0

                if (n_src_odi <= 0):
                    wcs_qual_string += "     ."
                elif (n_matched >= 1e5):
                    wcs_qual_string += " +++++"
                else:
                    wcs_qual_string += " %5d" % (n_matched)
            
            logger.debug("WCSDTL_%d --> %s" % (otay, wcs_qual_string))
            ota_list[0].header["WCSDTL_%d" % otay] = wcs_qual_string
                
        if (not ccmatched['valid_wcs_solution']):

            # This disabled the photometric calibration afterwards
            enough_stars_for_fixwcs = False
            global_reduction_log.fail('wcscal')

        else:
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
                                     array=odi_2mass_cat[:,  0], disp='F12.8'),
                       pyfits.Column(name='ODI_DEC', format='D', unit='degrees', 
                                     array=odi_2mass_cat[:,  1], disp='F12.8'),
                       pyfits.Column(name='REF_RA', format='D', unit='degrees',
                                     array=odi_2mass_cat[:, -2], disp='F12.8'),
                       pyfits.Column(name='REF_DEC', format='D', unit='degrees',
                                     array=odi_2mass_cat[:, -1], disp='F12.8'),
                       pyfits.Column(name='OTA', format='I', unit='',
                                     array=odi_2mass_cat[:, 8].astype(numpy.int), disp='I2.2'),
            ]
            coldefs = pyfits.ColDefs(columns)
            matchedhdu = pyfits.BinTableHDU.from_columns(coldefs)
            matchedhdu.name = "WCSCAL.CAT"
            matchedhdu.header['MATCHRAD'] = (2., "matching radius in arcsec")
            ota_list.append(matchedhdu)

            # Compute the WCS quality statistics
            # This also writes to the Primary and OTA-level headers
            wcs_quality = dev_ccmatch.global_wcs_quality(odi_2mass_cat, ota_list)
            # print "WCS quality =",wcs_quality
            ota_list[0].header['WCSXRMS'] = \
                wcs_quality['full']['RMS-RA'] if numpy.isfinite(wcs_quality['full']['RMS-RA']) else -1.
            ota_list[0].header['WCSYRMS'] = \
                wcs_quality['full']['RMS-DEC'] if numpy.isfinite(wcs_quality['full']['RMS-DEC']) else -1.

            global_reduction_log.success('wcscal')


        # Compute the image quality using all detected sources
        # Make sure to only include round source (elongation < 1.3) and decent 
        # photometry (i.e. high enoug S/N)
        good_fwhm_values = (global_source_cat[:, SXcolumn['flags']] == 0) & \
                           (global_source_cat[:, SXcolumn['elongation']] < 1.3) & \
                           (global_source_cat[:, SXcolumn['mag_err_3.0']] <= 0.2)
        seeing = global_source_cat[:, SXcolumn['fwhm_world']] #* 3600. # convert to arcsec
        seeing_all = numpy.median(seeing)
        seeing_clipped = three_sigma_clip(seeing)
        seeing_filtered = numpy.median(seeing_clipped)
        ota_list[0].header['FWHM_ALL'] = (seeing_all, "median FWHM of all sources")
        ota_list[0].header['FWHM_FLT'] = (seeing_filtered, "median FWHM, 3sigma clipped")
        ota_list[0].header['FWHMNFLT'] = (seeing_clipped.shape[0], "number of src in FWHM comp")

        ota_list[0].header['FWHMSTAR'] = (seeing_all, "median FWHM of SDSS-matched stars")
        ota_list[0].header['SEEING'] = (seeing_all, "Seeing [arcsec]")
        ota_list[0].header['SEEING_N'] = (seeing_clipped.shape[0], "number of stars in seeing comp")

    logger.debug("Next up: fixwcs & qaplots")
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
                                             matched_odierror=odi_2mass_matched[:, SXcolumn['mag_err_auto']],
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

        plotfilename = create_qa_filename(outputfile, "psfshape", options)
        podi_diagnosticplots.diagplot_psfshape_map(ra=global_source_cat[:, SXcolumn['ra']],
                                                   dec=global_source_cat[:, SXcolumn['dec']],
                                                   elongation=global_source_cat[:, SXcolumn['elongation']],
                                                   angle=global_source_cat[:, SXcolumn['position_angle']],
                                                   fwhm=global_source_cat[:, SXcolumn['fwhm_world']],
                                                   ota=global_source_cat[:, SXcolumn['ota']],
                                                   output_filename=plotfilename,
                                                   title=diagnostic_plot_title,
                                                   ota_outlines=ota_outlines,
                                                   options=options,
                                                   also_plot_singleOTAs=options['otalevelplots'],
                                                   title_info=title_info)

    podi_logging.ppa_update_progress(70, "Astrometric calibration complete")

    #
    # Add some default photometric calibration keywords
    #
    ota_list[0].header['PHOTMCAT'] = (None, "catalog used for photometric calibration")
    add_fits_header_title(ota_list[0].header, "Photometric calibration results", 'PHOTMCAT')
    ota_list[0].header['PHOTFILT'] = (None, "filter from reference catalog")
    
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

    logger.debug("Next up: photcalib")

    if (not options['photcalib']):
        # update global reduction log
        global_reduction_log.not_selected('photcal')

    elif ((options['photcalib'] and options['fixwcs'] and not enough_stars_for_fixwcs) or 
          (options['photcalib'] and not options['fixwcs'])):
        # update reduction log
        global_reduction_log.skip_after_fail('photcal')

    elif (options['photcalib'] 
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

        if (options['auto_photflat']):
            global_reduction_log.attempt('autophotflat')

            logger.info("Starting initial photometric calibration")
            photcalib_details = {}
            zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = \
                podi_photcalib.photcalib(global_source_cat, outputfile,
                                         filter_name,
                                         exptime=exptime,
                                         diagplots=options['create_qaplots'],
                                         plottitle=titlestring,
                                         otalist=ota_list,
                                         options=options,
                                         detailed_return=photcalib_details,
                                         plot_suffix="init")
            podi_photcalib.write_photcalib_headers(
                hdr=ota_list[0].header,
                zeropoint_median=zeropoint_median,
                zeropoint_std=zeropoint_std,
                zeropoint_exptime=zeropoint_exptime,
                filter_name=filter_name,
                photcalib_details=photcalib_details,
            )

            #
            # Calculate a photometric flat-field from the initial photometric
            # calibration data
            #
            logger.info("Calculating (auto-)photometric flatfield from data")
            raw_tbhdu = create_odi_sdss_matched_tablehdu(
                photcalib_details['odi_sdss_matched_raw'],
                photcalib_details=photcalib_details,
            )
            ota_list.append(raw_tbhdu)

            pf_hdu, pf_interpol = podi_photflat.create_photometric_flatfield(
                input_hdus=[pyfits.HDUList(ota_list)],
                strict_ota=False,
                return_interpolator=True,
            )

            # create the photflat diagnostic plot
            (fluxcorr, all_extnames) = pf_interpol
            photflat_diagplot_fn = create_qa_filename(outputfile, "autophotflat", options)
            photflat_diagplot_title = \
                "%(OBJECT)s: auto phot.flat.\n" \
                "%(OBSID)s -- %(FILTER)s -- %(EXPMEAS).1f" % ota_list[0].header
            # self-generated photometric flatfield
            podi_diagnosticplots.diagplot_photflat(
                extnames=all_extnames,
                data=fluxcorr,
                title=photflat_diagplot_title,
                output_filename=photflat_diagplot_fn,
                options=options,
                n_sigma=3, force_symmetric=False,
            )

            # change the extension name of the initial CAT.PHOTCALIB table to
            # avoid duplicates with the final CAT.PHOTCALIB table created below.
            raw_tbhdu.name = "CAT.PHOTCALIB.RAW"
            logger.info("#sources in catalog: %d" % (photcalib_details['odi_sdss_matched_raw'].shape[0]))

            photflat_allsuccess = True
            for ext in ota_list:
                if (not is_image_extension(ext)):
                    continue
                try:
                    ff = pf_hdu[ext.name].data
                    logger.debug("Applying photometric flat-field to OTA %s [sky: %f]" % (ext.name, sky_global_median))

                    if (sky_global_median is not None):
                        ext.data = (ext.data - sky_global_median) / ff + sky_global_median
                    else:
                        ext.data /= ff
                except KeyError:
                    logger.warning("No auto-photflat data for OTA %s" % (ext.name))
                    photflat_allsuccess = False
            if (photflat_allsuccess):
                global_reduction_log.success('autophotflat')
            else:
                global_reduction_log.partial_fail('autophotflat')

            #
            # Correct the photometry catalog using the information from the
            # photometric flatfield
            #
            logger.info("Correcting existing photometry for photometric flatfielding")
            # raw_cat = photcalib_details['odi_sdss_matched_raw']
            photflat_correction = podi_photflat.lookup_corrections(
                ota=global_source_cat[:, SXcolumn['ota']],
                x=global_source_cat[:, SXcolumn['x']],
                y=global_source_cat[:, SXcolumn['y']],
                photflat=pf_hdu,
            )
            # Now correct all photometry in the ODI source catalog
            # the next run of photcalib will then create the new matched
            # ODI+REF catalog using the newly corrected photometry
            for col in ['mag_aper_2.0','mag_aper_3.0','mag_aper_4.0',
                'mag_aper_5.0','mag_aper_6.0','mag_aper_8.0',
                'mag_aper_10.0','mag_aper_12.0','mag_auto']:
                global_source_cat[:,SXcolumn[col]] -= photflat_correction
            # print "Phtoflat corrections:\n", photflat_correction.shape, "\n", photflat_correction
            logger.debug("computed %d photometry corrections from photometric flatfield" % (
                photflat_correction.shape[0]))

            #
            # Finally, run the photometric calibration again, using the
            # corrected, photometrically flat-fielded data
            #
            logger.info("Starting final photometric calibration")


        if (not options['auto_photflat']):
            logger.info("Starting photometric calibration")

        photcalib_details = {}
        zeropoint_median, zeropoint_std, odi_sdss_matched, zeropoint_exptime = \
            podi_photcalib.photcalib(global_source_cat, outputfile, filter_name, 
                                     exptime=exptime,
                                     diagplots=options['create_qaplots'],
                                     plottitle=titlestring,
                                     otalist=ota_list,
                                     options=options,
                                     detailed_return=photcalib_details)

        podi_photcalib.write_photcalib_headers(
            hdr=ota_list[0].header,
            zeropoint_median=zeropoint_median,
            zeropoint_std=zeropoint_std,
            zeropoint_exptime=zeropoint_exptime,
            filter_name=filter_name,
            photcalib_details=photcalib_details,
        )

        # ota_list[0].header['PHOTMCAT'] = (photcalib_details['catalog'])
        # ota_list[0].header['PHOTFILT'] = (photcalib_details['reference_filter'])
        #
        # ota_list[0].header["PHOTZP"] = (zeropoint_median, "phot. zeropoint corr for exptime")
        # ota_list[0].header["PHOTZPSD"] = (zeropoint_std, "zeropoint std.dev.")
        # ota_list[0].header["PHOTZP_X"] = (zeropoint_exptime, "phot zeropoint for this frame")
        # ota_list[0].header["PHOTZPSP"] = (photcalib_details['zp_upper1sigma'], "phot ZP upper 1sigma limit")
        # ota_list[0].header["PHOTZPSM"] = (photcalib_details['zp_lower1sigma'], "phot ZP lower 1sigma limit")
        # ota_list[0].header["PHOTZPER"] = (photcalib_details['stderrofmean'], "phot ZP std.err of the mean")
        # ota_list[0].header["PHOTZP_N"] = (photcalib_details['n_clipped'], "number stars in clipped distrib.")
        # ota_list[0].header["PHOTZPN0"] = (photcalib_details['n_raw'], "total number of matched ref stars")
        #
        # ota_list[0].header["MAGZERO"] = (photcalib_details['median'], "phot. zeropoint corr for exptime")
        # ota_list[0].header["MAGZSIG"] = (photcalib_details['std'], "phot ZP dispersion")
        # ota_list[0].header["MAGZERR"] = (photcalib_details['stderrofmean'], "phot ZP uncertainty")

        if (photcalib_details['median'] > 0 and
            photcalib_details['median'] < 50):
            global_reduction_log.success('photcal')
        else:
            global_reduction_log.fail('photcal')

        flux_scaling = 1.0
        if (photcalib_details['n_raw'] > 0 and
            photcalib_details['median'] > 0 and
            photcalib_details['median'] < 50):
            flux_scaling = math.pow(10., -0.4*(zeropoint_exptime - 25.0))
        ota_list[0].header["FLXSCALE"] = flux_scaling
        # For swarp to work properly, also copy the FLXSCALE keyword into every image extension
        for i in range(len(ota_list)):
            if (is_image_extension(ota_list[i])):
                ota_list[i].header["FLXSCALE"] = flux_scaling

        # # Add some information on what apertures were used for the photometric calibration
        # ota_list[0].header['MAG0MODE'] = (photcalib_details['aperture_mode'], "how was aperture determined")
        # ota_list[0].header['MAG0SIZE'] = (photcalib_details['aperture_size'], "what aperture size was used")
        # ota_list[0].header['MAG0_MAG'] = (photcalib_details['aperture_mag'], "id string for magnitude")
        # ota_list[0].header['MAG0_ERR'] = (photcalib_details['aperture_magerr'], "is string for mag error")
        #
        # if (not photcalib_details['radialZPfit'] == None):
        #     ota_list[0].header['RADZPFIT'] = True
        #     ota_list[0].header['RADZP_P0'] = photcalib_details['radialZPfit'][0]
        #     ota_list[0].header['RADZP_P1'] = photcalib_details['radialZPfit'][1]
        #     ota_list[0].header['RADZP_E0'] = photcalib_details['radialZPfit_error'][0]
        #     ota_list[0].header['RADZP_E1'] = photcalib_details['radialZPfit_error'][1]
        #
        # if (not photcalib_details['zp_restricted'] == None):
        #     (sel_median, sel_std, sel_psigma, sel_msigma, sel_n, sel_medodimag, sel_maxodimag, sel_minodimag) = photcalib_details['zp_restricted']
        #     ota_list[0].header['ZPRESMED'] = sel_median
        #     ota_list[0].header['ZPRESSTD'] = sel_std
        #     ota_list[0].header['ZPRES_SP'] = sel_psigma
        #     ota_list[0].header['ZPRES_SM'] = sel_msigma
        #     ota_list[0].header['ZPRES__N'] = sel_n
        #     ota_list[0].header['ZPRES_MD'] = sel_medodimag
        #     ota_list[0].header['ZPRES_MX'] = sel_maxodimag
        #     ota_list[0].header['ZPRES_MN'] = sel_minodimag
        #
        # if (not photcalib_details['zp_magnitude_slope'] == None):
        #     fit, uncert = photcalib_details['zp_magnitude_slope']
        #     ota_list[0].header['ZPSLP_P0'] = fit[0]
        #     ota_list[0].header['ZPSLP_P1'] = fit[1]
        #     ota_list[0].header['ZPSLP_E0'] = uncert[0]
        #     ota_list[0].header['ZPSLP_E1'] = uncert[1]
        #
        # ref_ZP = -99. if not filter_name in reference_zeropoint else reference_zeropoint[filter_name][0]
        # ota_list[0].header['MAGZREF'] = (ref_ZP, "reference photometric zeropoint")
        #
        # # Also compute the zeropoint after correction for airmass
        # zp_airmass1 = -99.
        # if (filter_name in atm_extinction):
        #     zp_airmass1 = zeropoint_median + (ota_list[0].header['AIRMASS']-1) * atm_extinction[filter_name]
        # ota_list[0].header['MAGZ_AM1'] = (zp_airmass1, "phot Zeropoint corrected for airmass")
        #
        # # Add some information whether or not we performed a color-term correction
        # colorterm_correction = (not photcalib_details['colorterm'] == None)
        # ota_list[0].header['MAGZ_CT'] = colorterm_correction
        # ota_list[0].header['MAGZ_COL'] = photcalib_details['colorcorrection'] if colorterm_correction else ""
        # ota_list[0].header['MAGZ_CTC'] = photcalib_details['colorterm'] if colorterm_correction else 0.0

        # Compute the sky-brightness 
        sky_arcsec = sky_global_median / (0.11**2) # convert counts/pixel to counts/arcsec*2
        sky_mag = -99.
        if (sky_arcsec > 0 and zeropoint_exptime < 99):
            sky_mag = -2.5 * math.log10(sky_arcsec) + zeropoint_exptime
        ota_list[0].header['SKYMAG'] = sky_mag

        master_reduction_files_used = podi_associations.collect_reduction_files_used(
            master_reduction_files_used, {'photcalib-reference': photcalib_details['reference_catalog_files']}
        )

        # Now convert the matched source catalog into a valid FITS extension 
        # and add it to the output.
        if (odi_sdss_matched is not None and odi_sdss_matched.shape[0] > 0):
            logger.debug("Adding matched SDSS=ODI source catalog to output as FITS extension")

            ref_tbhdu = create_odi_sdss_matched_tablehdu(
                photcalib_details['odi_sdss_matched_raw'],
                photcalib_details=photcalib_details,
                extname="CAT.PHOTREF",
            )
            logger.info("#sources in catalog: %d" % (
                photcalib_details['odi_sdss_matched_raw'].shape[0]))
            ota_list.append(ref_tbhdu)

            match_tablehdu = create_odi_sdss_matched_tablehdu(odi_sdss_matched, photcalib_details)
            logger.debug("final photcalib cat: %d" % (odi_sdss_matched.shape[0]))
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
                            'PHOTMCAT', 'PHOTFILT',
                            'MAG0MODE', 'MAG0SIZE', 'MAG0_MAG', 'MAG0_ERR',
            ]:
                if (hdrname in ota_list[0].header):
                    match_tablehdu.header[hdrname] = ota_list[0].header[hdrname]
            ota_list.append(match_tablehdu)

            # Also use the matched catalog to determine the seeing of only stars
            star_seeing = odi_sdss_matched[:, SXcolumn['fwhm_world']+2]
            print(star_seeing)
            cleaned = three_sigma_clip(star_seeing)
            seeing = numpy.median(cleaned)
            logger.debug("Seeing is %.2f arcsec" % (seeing))
            ota_list[0].header['FWHMSTAR'] = (seeing, "median FWHM of SDSS-matched stars")
            ota_list[0].header['SEEING'] = (seeing, "Seeing [arcsec]")
            ota_list[0].header['SEEING_N'] = (cleaned.shape[0], "number of stars in seeing comp")

            # Compute a approximate detection limit
            # This assumes an aperture with diameter of 2x the seeing 
            # (i.e. radius = seeing).
            aperture_area = (seeing / 0.11)**2 * math.pi
            readnoise = 8.
            if (sky_global_median >= 0):
                bgcounts = aperture_area * (sky_global_median + readnoise**2)
                counts_sn1 = math.sqrt(bgcounts)
                limiting_mag = -2.5*math.log10(counts_sn1) + zeropoint_exptime
            else:
                limiting_mag = -99.
            ota_list[0].header['PHOTDPTH'] = limiting_mag

    # Create the image quality plot
    # Also create a diagnostic plot for the Seeing.
    # choose the sdssm-matched catalog if available, otherwise use the raw source catalog.
    logger.debug("Next up: photcalib & seeingplots")
    logger.debug("fixwcs: %s, create_qaplots: %s, photcalib: %s, enough_stars_for_fixwcs: %s" % (
        str(options['fixwcs']), str(options['create_qaplots']), str(options['photcalib']), str(enough_stars_for_fixwcs)
    ))
    if (options['fixwcs'] and 
        options['create_qaplots'] and 
        options['photcalib'] and
        enough_stars_for_fixwcs):

        if (options['photcalib'] and odi_sdss_matched is not None and odi_sdss_matched.shape[0] > 0):
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



    #
    # Convert the sky samples into a FITS table extension and append it to the output
    # Also make sure to compute real Ra/Dec values from OTA, X, and Y
    #
    logger.debug("Creating the SKYLEVEL table extension")
    if (sky_samples_global is not None):
        sky_samples_final = numpy.empty((0,sky_samples_global.shape[1]))
        sky_otas = set(sky_samples_global[:,-1])

        # print sky_otas
        logger.debug("SKY=OTAs: %s" % (str(sky_otas)))
        for sky_ota in sky_otas:
            sky_extname = "OTA%02d.SCI" % (sky_ota)
            # find the right HDU extension to compute the WCS solution
            for cur_ext in ota_list:
                if ('EXTNAME' in cur_ext.header and
                    cur_ext.header['EXTNAME'] == sky_extname):
                    logger.debug("Computing Ra/Dec for sky-samples in %s" % (sky_extname))
                    # Found the right OTA
                    in_this_ota = (sky_samples_global[:, -1] == sky_ota)
                    if (numpy.sum(in_this_ota) <= 0):
                        continue
                    sky_ota = sky_samples_global[in_this_ota]
                    wcs = astWCS.WCS(cur_ext.header, mode='pyfits')
                    ota_radec = numpy.array(wcs.pix2wcs(sky_ota[:,2]-1.0, sky_ota[:,3]-1.0))
                    # print ota_radec
                    sky_ota[:, 0:2] = ota_radec
                    sky_samples_final = numpy.append(sky_samples_final, sky_ota, axis=0)
                    break
        sky_columns = [
            # RA of box center
            pyfits.Column(name='RA', format='D', unit='degrees',
                          array=sky_samples_final[:,0], disp='F10.6'),
            # DEC of box center
            pyfits.Column(name='DEC', format='D', unit='degrees',
                          array=sky_samples_final[:,1], disp='F10.6'),
            # X position of box center
            pyfits.Column(name='X', format='D', unit='pixel',
                          array=sky_samples_final[:,2], disp='F10.4'),
            # Y position of box center
            pyfits.Column(name='Y', format='D', unit='pixel',
                          array=sky_samples_final[:,3], disp='F10.4'),
            # median counts in box
            pyfits.Column(name='INTENSITY', format='D', unit='counts',
                          array=sky_samples_final[:,4], disp='E12.5E2'),
            # distance to closest source
            pyfits.Column(name='MIN_D', format='D', unit='pixel',
                          array=sky_samples_final[:,2], disp='F8.3'),
            # source OTA
            pyfits.Column(name='OTA', format='I2', unit='',
                          array=sky_samples_final[:,-1], disp='I2.2'),
        ]
        sky_coldefs = pyfits.ColDefs(sky_columns)
        sky_tbhdu = pyfits.BinTableHDU.from_columns(sky_coldefs)
        sky_tbhdu.name = "SKYLEVEL"
        # Copy a bunch of headers from the primary HDUu to the SKYLEVEL hdu
        for key in [
                "SKYLEVEL", "SKYBG", 'SKYBGCLP', 'SKYBGMIN', 'SKYBGSTD', 'SKYSMPLS',
                'SKYL1SIG', 'SKYU1SIG', 'SKYL2SIG', 'SKYU2SIG', 'SKY_LO5P', 'SKY_LO5S',
                ]:
            sky_tbhdu.header[key] = ota_list[0].header[key]
        ota_list.append(sky_tbhdu)
        # numpy.savetxt("skyfinal", sky_samples_final)

    logger.debug("Saving reduction log for non-sidereal corrections")
    if (options['nonsidereal'] is None):
        global_reduction_log.not_selected('nonsidereal')
    else:
        logger.info("Starting non-sidereal WCS modification")
        apply_nonsidereal_correction(ota_list, options, logger, 
                                     reduction_log=global_reduction_log)

        if ('ref' in options['nonsidereal'] and
            os.path.isfile(options['nonsidereal']['ref'])):
            master_reduction_files_used = podi_associations.collect_reduction_files_used(
                master_reduction_files_used, 
                {"nonsidereal-reference": options['nonsidereal']['ref']})



    #
    # Add some information about the filter bandpasses to the output file
    #
    logger.debug("Adding filter bandpass information")
    filtername = ota_list[0].header['FILTER']
    bandpass = ("???", 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,) if (not filtername in filter_bandpass) else filter_bandpass[filtername]
    (filter_def, filter_mean_pos, filter_center_pos, 
     filter_fwhm, filter_fwhm_left, filter_fwhm_right, 
     filter_max, filter_mean, filter_area, filter_left5, filter_right5) = bandpass
    # Keywords for compatibility with AuCaP
    ota_list[0].header['PHOTCLAM'] = (filter_center_pos/10., "central wavelength of filter [nm]")
    ota_list[0].header['PHOTBW']   = (filter_fwhm/10.,       "RMS width of filter [nm]")
    ota_list[0].header['PHOTFWHM'] = (filter_fwhm/10.,       "FWHM of filter [nm]")
    # Extra keywords with more specifics
    ota_list[0].header['FILTFILE'] = (filter_def,            "filename of filter definition")
    ota_list[0].header['FILTAVGP'] = (filter_mean_pos/10.,   "weighted center position [nm]")
    ota_list[0].header['FILTCTRP'] = (filter_center_pos/10., "filter center [nm]")
    ota_list[0].header['FILTFWHM'] = (filter_fwhm/10.,       "filter FWHM [nm]")
    ota_list[0].header['FILTBLUE'] = (filter_fwhm_left/10.,  "blue filter edge at half max [nm]")
    ota_list[0].header['FILTRED']  = (filter_fwhm_right/10., "red filter edge at half max [nm]")
    ota_list[0].header['FILTMAXT'] = (filter_max,            "filter maximum transmission")
    ota_list[0].header['FILTAVGT'] = (filter_mean,           "filter average transmission")
    ota_list[0].header['FILTEQWD'] = (filter_area/10.,       "filter equivalent width [nm]")
    ota_list[0].header['FILTBLU5'] = (filter_left5/10.,      "blue filter edge at 5% max [nm]")
    ota_list[0].header['FILTRED5'] = (filter_right5/10.,     "red filter edge at 5% max [nm]")
    add_fits_header_title(ota_list[0].header, "Filter bandpass definitions", 'PHOTCLAM')


    # Now all processes have returned their results, terminate them 
    # and delete all instances to free up memory
    logger.debug("Joining all processes to make sure they terminate and free memory")
    for cur_process in processes:
        logger.debug("Joining process %s, PID %d" % (cur_process, cur_process.pid))
        cur_process.join(timeout=0.1)
        if (cur_process.is_alive()):
            logger.debug("Process is still alive, terminate it and re-join")
            cur_process.terminate()
            cur_process.join(timeout=0.1)


    #
    # Create an association table from the master reduction files used.
    #
    logger.debug("Creating master association table and adding it to FITS")
    master_reduction_files_used = podi_associations.collect_reduction_files_used(master_reduction_files_used, 
                                                               additional_reduction_files)
    assoc_table = podi_associations.create_association_table(master_reduction_files_used)
    ota_list.append(assoc_table)


    #
    # Add some almanac data to the current frame
    #
    logger.debug("Adding almanac data to primary FITS header")
    podi_almanac.add_ephem_data_to_header(ota_list[0].header, None)

    #
    # Check if we find a shift-history, and if so, create the shift history plot
    #
    logger.debug("Create shift-history plot")
    used_ot_shifting = False
    for _ox, _oy in fpl.available_ota_coords:
        _ota = _ox * 10 + _oy
        shift_history_file = "%s/%s.%02d_shift.fits" % (directory, filebase, _ota)
        if (os.path.isfile(shift_history_file)):
            used_ot_shifting = True
            break
        else:
            shift_history_file += ".fz"
            if (os.path.isfile(shift_history_file)):
                used_ot_shifting = True
                break
    if (used_ot_shifting):
        logger.info("Creating shift-history plot from %s" % (shift_history_file))
        plottitle = "Shift-history for %(OBSID)s (%(OBJECT)s, %(FILTER)s, %(EXPTIME)ds, #shifts=%%(nshifts)d)" \
                    % ota_list[0].header
#            obsid, object, , exptime)
        podi_shifthistory.view_shift_history(
            filename=shift_history_file,
            plot_filename=outputfile[:-5]+".otshift",
            extension_list=options['plotformat'],
            title=plottitle,)

    logger.debug("Writing global reduction log to primary FITS header")
    global_reduction_log.write_to_header(ota_list[0].header)

    #print "Waiting for a bit"
    #afw.wait()
    #print "done waiting, writing output file"
    #print ota_list
    if (options['compressed_hdu']):
        for iext, ext in enumerate(ota_list):
            if (type(ext) == pyfits.hdu.image.ImageHDU):
                logger.info("Compressing %s" % (ext.name))
                ch = pyfits.CompImageHDU(data=ext.data,
                                         header=ext.header,
                                         name=ext.name,
                                         compression_type='RICE_1',
                                         quantize_level=4,)
                ota_list[iext] = ch

    hdulist = pyfits.HDUList(ota_list)
    for i in range(1, len(hdulist)):
        if 'SIMPLE' in hdulist[i].header:
            del hdulist[i].header['SIMPLE']
    hdulist.verify()

    #print "hdulist=",hdulist

    podi_logging.ppa_update_progress(80, "Reduction and calibration complete")

    if (outputfile is not None):
        logger.debug("Complete, writing output file %s" % (outputfile))
        start_time = time.time()
        clobberfile(outputfile)
        hdulist.writeto(outputfile, overwrite=True, output_verify='ignore')
        end_time = time.time()
        logger.debug("All work completed successfully, output written to %s (took %.3f seconds)" % (
            outputfile, (end_time-start_time)))

    logger.debug("Unstaging data!")
    unstage_data(options, staged_data, input)
    logger.debug("Done unstaging data!")

    # podi_logging.print_stacktrace()

    if (batchmode):
        logger.info("All work completed successfully, parsing output for further processing")
        return hdulist

    #     clobberfile(outputfile)
    #     hdulist.writeto(outputfile, overwrite=True)
    #     # afw.write(hdulist, outputfile)

    # else:
            

    # # afw.finish(userinfo=True)

    logger.debug("Work done, returning output filename: %s" % (outputfile))
    return outputfile


def apply_nonsidereal_correction(ota_list, options, logger=None, reduction_log=None):
    # This WCS in this frame should be corrected for non-sidereal tracking
    # Tracking rates are given in arcseconds per hour
    # Note that d_ra is given as dRA*cos(dec)

    if (logger is None):
        logger = logging.getLogger("NonsiderealCorr")

    if (options['nonsidereal'] is None):
        logger.debug("No nonsidereal option set, skipping non-sidereal correction")
        return
    if (not ('dra' in options['nonsidereal'] and
             'ddec' in options['nonsidereal'] and
             'ref_mjd' in options['nonsidereal'])):
        logger.debug("One of the nonsidereal option missing, skipping non-sidereal correction")
        logger.debug("Available options are:"+str(options['nonsidereal']))
        if (reduction_log is not None):
            reduction_log.fail('nonsidereal')
        return
             
    mjd = ota_list[0].header['MJD-OBS']
    delta_t_hours = (mjd - options['nonsidereal']['ref_mjd']) * 24.
    dra_t = options['nonsidereal']['dra'] * delta_t_hours / 3600.
    ddec_t = options['nonsidereal']['ddec'] * delta_t_hours / 3600.

    ota_list[0].header['NSIDPMRA'] = (options['nonsidereal']['dra'], "proper motion dRa*cos(dec) [arcsec/hr]")
    ota_list[0].header['NSIDPMDE'] = (options['nonsidereal']['ddec'], "proper motion dDec [arcsec/hr]")
    ota_list[0].header['NSIDBASE'] = (options['nonsidereal']['ref_mjd'], "MJD of reference frame")
    ota_list[0].header['NSIDDMJD'] = (mjd - options['nonsidereal']['ref_mjd'], "time diff to ref. frame [days]")
    ota_list[0].header['NSID_DHR'] = (delta_t_hours, "time diff to ref. frame [hours]")
    add_fits_header_title(ota_list[0].header, "Non-sidereal correction", 'NSIDPMRA')
    logger.debug("Tracking rates are dRA=%(dra)f dDEC=%(ddec)f arcsec/hour" % options['nonsidereal'])
    logger.debug("Time-offset to reference frame: %f hours" % (delta_t_hours))

    for ext in range(len(ota_list)):
        if (not is_image_extension(ota_list[ext])):
            continue
        declination = ota_list[ext].header['CRVAL2']
        dra_corrected = dra_t / math.cos(math.radians(declination))

        ota_list[ext].header['CRVAL1'] -= dra_corrected
        ota_list[ext].header['CRVAL2'] -= ddec_t

        if ('NSIDDRA'in ota_list[ext].header or
            'NSIDDDEC' in ota_list[ext].header):
            logger.error("This OTA already has a non-sidereal correction applied!")
            continue

        ota_list[ext].header['NSIDDRA'] = dra_corrected
        ota_list[ext].header['NSIDDDEC'] = ddec_t
        logger.debug("Adding offset of %f %f arcsec to OTA %s" % (
            dra_corrected*3600., ddec_t*3600, ota_list[ext].header['EXTNAME'])
        )

    if (reduction_log is not None):
        reduction_log.success('nonsidereal')
    return [dra_t, ddec_t]


    
def create_odi_sdss_matched_tablehdu(odi_sdss_matched, photcalib_details=None,
                                     extname=None):

    # Note:
    # The +2 in the collumn indices accounts for the fact that the Ra/Dec 
    # from the SDSS catalog are inserted after the Ra/Dec of the ODI 
    # catalog, shifting all other columns to the right by two positions.
    columns = [\
        #
        # Ra/Dec from ODI catalog
        #
        pyfits.Column(name='ODI_RA', disp='F12.8',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 0]),
        pyfits.Column(name='ODI_DEC', disp='F12.8',
                      format='D', unit='degrees', 
                      array=odi_sdss_matched[:, 1]),

        ]

    if (photcalib_details['catalog'] == "SDSS"):
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

    
        #
        # Ra/Dec from SDSS
        #
        columns.append(pyfits.Column(name='SDSS_RA', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 2]))
        columns.append(pyfits.Column(name='SDSS_DEC', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 3]))

        #
        # Now add the SDSS magnitudes
        #
        columns.append(pyfits.Column(name='SDSS_MAG_U', disp='F7.4',
                                     format='E', unit='mag', 
                                     array=odi_sdss_matched[:, SDSScolumn['u']]))
        columns.append(pyfits.Column(name='SDSS_ERR_U', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['u_err']]))

        columns.append(pyfits.Column(name='SDSS_MAG_G', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['g']]))
        columns.append(pyfits.Column(name='SDSS_ERR_G', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['g_err']]))

        columns.append(pyfits.Column(name='SDSS_MAG_R', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['r']]))
        columns.append(pyfits.Column(name='SDSS_ERR_R', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['r_err']]))

        columns.append(pyfits.Column(name='SDSS_MAG_I', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['i']]))
        columns.append(pyfits.Column(name='SDSS_ERR_I', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['i_err']]))

        columns.append(pyfits.Column(name='SDSS_MAG_Z', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['z']]))
        columns.append(pyfits.Column(name='SDSS_ERR_Z', disp='F7.4',
                      format='E', unit='mag', 
                      array=odi_sdss_matched[:, SDSScolumn['z_err']]))
        # end SDSS

    elif (photcalib_details['catalog'] == "UCAC4"):
        #
        # Ra/Dec from SDSS
        #
        columns.append(pyfits.Column(name='UCAC_RA', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 2]))
        columns.append(pyfits.Column(name='UCAC_DEC', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 3]))

        catalog_columns = ['RA', 'DEC',
                           'MAG_UCAC', 'ERR_UCAC',
                           'MAG_B', 'ERR_B',
                           'MAG_V', 'ERR_V',
                           'MAG_G', 'ERR_G',
                           'MAG_R', 'ERR_R',
                           'MAG_I', 'ERR_I',
                       ]
        UCACcolumn = {}
        for name in catalog_columns:
            UCACcolumn[name] = len(UCACcolumn) + len(SXcolumn)

        ucac_columns = [ ('UCAC_MAG', "UCAC photometry magnitude", 'MAG_UCAC'), 
                         ('UCAC_ERR', "UCAC photometry magnitude error", 'ERR_UCAC'), 
                         ('APASS_MAG_B', "UCAC/APASS magnitude B", 'MAG_B'), 
                         ('APASS_ERR_B', "UCAC/APASS mag error B", 'ERR_B'), 
                         ('APASS_MAG_V', "UCAC/APASS magnitude V", 'MAG_V'), 
                         ('APASS_ERR_V', "UCAC/APASS mag error V", 'ERR_V'), 
                         ('APASS_MAG_g', "UCAC/APASS magnitude g", 'MAG_G'), 
                         ('APASS_ERR_g', "UCAC/APASS mag error g", 'ERR_G'), 
                         ('APASS_MAG_r', "UCAC/APASS magnitude r", 'MAG_R'), 
                         ('APASS_ERR_r', "UCAC/APASS mag error r", 'ERR_R'), 
                         ('APASS_MAG_i', "UCAC/APASS magnitude i", 'MAG_I'), 
                         ('APASS_ERR_i', "UCAC/APASS mag error i", 'ERR_I'), 
                    ]
        for (name, disp, col) in ucac_columns:
            columns.append(pyfits.Column(name=name, disp='F7.4',
                                         format='E', unit='mag', 
                                         array=odi_sdss_matched[:, UCACcolumn[col]]))

        # end UCAC
    elif (photcalib_details['catalog'] == "IPPRef"):
        #
        # Ra/Dec from IPPRef
        #
        columns.append(pyfits.Column(name='IPP_RA', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 2]))
        columns.append(pyfits.Column(name='IPP_DEC', disp='F12.8',
                                     format='D', unit='degrees', 
                                     array=odi_sdss_matched[:, 3]))

        catalog_columns = ['RA', 'DEC',
                           'MAG_G', 'ERR_G',
                           'MAG_R', 'ERR_R',
                           'MAG_I', 'ERR_I',
                           'MAG_Z', 'ERR_Z',
                       ]
        IPPcolumn = {}
        for name in catalog_columns:
            IPPcolumn[name] = len(IPPcolumn) + len(SXcolumn)

        ipp_columns = [ ('IPP_MAG_G', "synth. IPP magnitude g", 'MAG_G'), 
                        ('IPP_ERR_G', "synth. IPP mag error g", 'ERR_G'), 
                        ('IPP_MAG_R', "synth. IPP magnitude r", 'MAG_R'), 
                        ('IPP_ERR_R', "synth. IPP mag error r", 'ERR_R'), 
                        ('IPP_MAG_I', "synth. IPP magnitude i", 'MAG_I'), 
                        ('IPP_ERR_I', "synth. IPP mag error i", 'ERR_I'), 
                        ('IPP_MAG_Z', "synth. IPP magnitude z", 'MAG_Z'), 
                        ('IPP_ERR_Z', "synth. IPP mag error z", 'ERR_Z'), 
                    ]
        for (name, disp, col) in ipp_columns:
            columns.append(pyfits.Column(name=name, disp='F7.4',
                                         format='E', unit='mag', 
                                         array=odi_sdss_matched[:, IPPcolumn[col]]))

        # end UCAC

    elif (photcalib_details['catalog'] in sitesetup.catalog_directory and
          photcalib_details['catalog'] in sitesetup.catalog_mags):

        # print "SXcolumns:", len(SXcolumn)
        # print "cat cols:", len(sitesetup.catalog_mags[photcalib_details['catalog']])
        # print odi_sdss_matched.shape
        # numpy.savetxt("photcal.cat", odi_sdss_matched)
        for i_key, key in enumerate(sitesetup.catalog_mags[photcalib_details['catalog']]):
            if (key == ''):
                continue
            table_key = ('ref_%s' % (key)).upper()

            if (key in ['ra', 'dec']):
                col = i_key + 2
                unit = 'degrees'
            else:
                col = i_key + len(SXcolumn)
                unit = 'mag'

            columns.append(pyfits.Column(name=table_key, format='D', unit=unit,
                                     array=odi_sdss_matched[:, col]))

    else:
        logger.warning("Catalog not properly configured for output as TABLE extension")
        pass


    columns.append(pyfits.Column(name='ODI_FWHM', disp='F10.8',
                                 format='D', unit='degrees', 
                                 array=odi_sdss_matched[:, SXcolumn['fwhm_world']+2]))

    columns.append(pyfits.Column(name='ODI_MAG_AUTO', format='E', unit='mag',
                                 array=odi_sdss_matched[:,SXcolumn['mag_auto']+2],
                                 disp='F7.4'))
    columns.append(pyfits.Column(name='ODI_ERR_AUTO', format='E', unit='mag',
                                 array=odi_sdss_matched[:,SXcolumn['mag_err_auto']+2],
                                 disp='F6.4'))
        


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
                                     disp='F8.4'
                                     )
        )
        columns.append(pyfits.Column(name='ODI_ERR_D%02d' % (int(apsize*10.)),
                                     format='E', unit='mag',
                                     array=odi_sdss_matched[:,SXcolumn[col_magerr]+2],
                                     disp='F8.4'
                                 )
        )


    #
    # ODI source position in detector coordinates
    #
    columns.append(pyfits.Column(name='ODI_X',
                                     format='E', unit='pixel',
                                     array=odi_sdss_matched[:,SXcolumn['x']+2],
                                     disp='F7.2'
                                 )
    )
    columns.append(pyfits.Column(name='ODI_Y',
                                     format='E', unit='pixel',
                                     array=odi_sdss_matched[:,SXcolumn['y']+2],
                                     disp='F7.2'
                                 )
    )
    columns.append(pyfits.Column(name='ODI_OTA',
                                     format='I', unit='',
                                     array=odi_sdss_matched[:,SXcolumn['ota']+2].astype(numpy.int),
                                     disp='I2.2'
                                 )
    )
    if (photcalib_details['use_for_calibration_mask'] is not None):
        columns.append(pyfits.Column(name='USED4CALIB',
                                         format='L', unit='',
                                         array=photcalib_details['use_for_calibration_mask'],
                                     )
        )

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)

    if (extname is None):
        tbhdu.name = "CAT.PHOTCALIB"
    else:
        tbhdu.name = extname

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
        pyfits.Column(name='RA',  format='D', disp='F12.8', array=catalog[:,0]),
        pyfits.Column(name='DEC', format='D', disp='F12.8', array=catalog[:,1]),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)

    tbhdu.name = "CAT.WCSREF"
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
                             disp='F12.8'),
               pyfits.Column(name='DEC',            format='D', unit='degrees',
                             array=source_cat[:,SXcolumn['dec']], 
                             disp='F12.8'),
               pyfits.Column(name='X',              format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['x']],
                             disp='F8.3'),
               pyfits.Column(name='Y',              format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['y']], 
                             disp='F8.3'),
               pyfits.Column(name='FWHM_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['fwhm_image']],
                             disp='F6.2'),
               pyfits.Column(name='FWHM_WORLD',     format='E', unit='arcsec',
                             array=source_cat[:,SXcolumn['fwhm_world']],
                             disp='F10.6'),
               pyfits.Column(name='BACKGROUND',     format='E', unit='counts',
                             array=source_cat[:,SXcolumn['background']],
                             disp='F8.1'),
               pyfits.Column(name='FLAGS',          format='I', unit='',
                             array=source_cat[:,SXcolumn['flags']],
                             disp='B8.8'),
               pyfits.Column(name='OTA',            format='I', unit='',
                             array=source_cat[:,SXcolumn['ota']],
                             disp='I2.2'),

               pyfits.Column(name='MAG_AUTO',        format='E', unit='mag',
                             array=source_cat[:,SXcolumn['mag_auto']],
                             disp='F8.4'),
               pyfits.Column(name='MAGERR_AUTO',     format='E', unit='mag',
                             array=source_cat[:,SXcolumn['mag_err_auto']],
                             disp='F8.4'),
        
               pyfits.Column(name='FLUX_MAX',       format='E', unit='counts',
                             array=source_cat[:,SXcolumn['flux_max']],
                             disp='F8.1'),
               pyfits.Column(name='AWIN_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['major_axis']],
                             disp='F8.3'),
               pyfits.Column(name='BWIN_IMAGE',     format='E', unit='pixel',
                             array=source_cat[:,SXcolumn['minor_axis']],
                             disp='F8.3'),
               pyfits.Column(name='THETAWIN_IMAGE', format='E', unit='degrees',
                             array=source_cat[:,SXcolumn['position_angle']],
                             disp='F5.1'),
               pyfits.Column(name='ELONGATION',     format='E', unit='',
                             array=source_cat[:,SXcolumn['elongation']],
                             disp='F8.4'),
               pyfits.Column(name='ELLIPTICITY',    format='E', unit='',
                             array=source_cat[:,SXcolumn['ellipticity']],
                             disp='F7.5'),
           ]

    # Add all aperture photometry
    for i in range(len(SXapertures)):
        apsize = SXapertures[i]
        col_mag = "mag_aper_%0.1f" % (apsize)
        col_magerr = "mag_err_%0.1f" % (apsize)
        columns.append(pyfits.Column(name='MAG_D%02d' % (int(apsize*10.)), 
                                     format='E', unit='mag',
                                     array=source_cat[:,SXcolumn[col_mag]],
                                     disp='F8.4'
                                 )
        )
        columns.append(pyfits.Column(name='MAGERR_D%02d' % (int(apsize*10.)),
                                     format='E', unit='mag',
                                     array=source_cat[:,SXcolumn[col_magerr]],
                                     disp='F8.4'
                                 )
        )

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)

    tbhdu.name = "CAT.ODI"
    return tbhdu

















returnvalue_meaning = {
    0: "OK",
    1: "timeout",
    2: "exception",
    3: "unknown",
    4: "profiler",
}


if __name__ == "__main__":

    if (len(sys.argv) <= 1 or sys.argv[1] == "-help"):
        #print help('podi_matchpupilghost')
        import podi_collectcells as me
        print(me.__doc__)
        sys.exit(0)

    # m = multiprocessing.Manager()
    # process_tracker = m.Queue()

    # Read the input directory that contains the individual OTA files
    input = get_clean_cmdline()[1]

    # Assign a fallback output filename if none is given 
    if (len(get_clean_cmdline())>2):
        outputfile = get_clean_cmdline()[2]
    else:
        print("No output filename has been given, setting to default mergedcells.fits")
        outputfile = "mergedcells"
    print("Writing results into",outputfile)

    # Set the options for collectcells to some reasonable start values
    options = set_default_options()

    # Setup everything we need for logging
    podi_logging.setup_logging(options)

    # Then read the actual given parameters from the command line
    options = read_options_from_commandline(options)

    #
    # Initialize the log for easier debugging later on (if needed)
    #
    logger = logging.getLogger("CCRoot")
    logger.debug("COLLECTCELLS STARTUP")
    logger.debug("Execute command:\n%s" % (" ".join(sys.argv)))
    logger.debug("Current working directory: %s" % (os.getcwd()))
    logger.debug("Executable: %s" % (os.path.abspath(__file__)))
    logger.debug("COLLECTCELLS INPUT: %s (%s)" % (input, os.path.abspath(input)))


    #
    # Collect all cells, perform reduction and write result file
    #
    retvalue = 3
    try:
        if (cmdline_arg_isset('-profile')):
            options['profile'] = True
            import cProfile, pstats
            cProfile.run("""collectcells(input, outputfile,
                     options=options)""", "profiler")
            p = pstats.Stats("profiler")
            p.strip_dirs().sort_stats('time').print_stats()
            p.sort_stats('time').print_stats()
            retvalue = 4
        else:
            # print "collectcells with timeout:", cmdline_arg_isset("-timeout")
            #time.sleep(5)

            if (cmdline_arg_isset("-timeout")):
                
                timeout = float(cmdline_arg_set_or_default("-timeout", 900))
                print("Setting timeout to",timeout,"seconds")
                retvalue = collectcells_with_timeout(input, outputfile, options=options,
                                                     timeout=timeout,)
            else:
                start_time = time.time()
                collectcells(input, outputfile, 
                             #process_tracker=process_tracker, 
                             options=options)
                end_time = time.time()
                logger.debug("collectcells returned to __main__ after %.3f seconds" % ((end_time-start_time)))
                retvalue = 0
    except:
        print("Cleaning up left over child processes")
        podi_logging.log_exception()
        #kill_all_child_processes(process_tracker)
        retvalue = 2

    finally:
        #
        # Send some debug.log closing statement and shutdown all logging
        #
        logger.debug("COLLECTCELLS RETURN-VALUE: %d (%s)" % (retvalue, returnvalue_meaning[retvalue]))
        logger.debug("COLLECTCELLS SHUTDOWN")
        podi_logging.shutdown_logging(options)

    #
    # Adding some final information before shutting down
    # This should help find the problem inside PPA
    #
    # time.sleep(1)
    # print "All threads should be closed now!"
    # podi_logging.print_stacktrace(stdout=True)

    #
    # return the return value as determined above to let the calling program 
    # know if all execution was completed successfully or if there were some 
    # errors/problems during execution
    #
    sys.exit(retvalue)
