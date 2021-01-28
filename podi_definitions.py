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

This module contains a number of global constants required during reduction, as
well as functions that are widely used throughout the other modules.  Some
general definitions useful in a number of podi scripts

"""

import sys
import os
import numpy
import ctypes
import math
import numpy
import astropy.io.fits as pyfits
import subprocess
import scipy
import scipy.ndimage
import scipy.special
import itertools
import logging
import bottleneck
from bottleneck import nanmean, nanmedian

from podi_commandline import *
import podi_sitesetup as sitesetup
from wiyn_filters import *

#
#
#


############
#                  | | |
#  Experimental    | | |
#                  V V V
############

#pyfits.USE_MEMMAP = False

############
#                  A A A
#  Experimental    | | |
#                  | | |
############






sdss_photometric_column = {"u":  2,
                           "g":  4,
                           "r":  6,
                           "i":  8,
                           "z": 10,
                           }
ucac_photometric_column = {"UCAC_Red": 2,
                           "B":  4,
                           "V":  6,
                           "g":  8,
                           "r": 10,
                           "i": 12,
                           }

cellmode_ids = {
    "S": 0,
    "V": 1,
    "D": 1,
    # Add other cellmodes and some numerical representation here
}




#
# Source Extractor flags
# 
# 0b00000001 =   1: near neighbor
# 0b00000010 =   2: deblended source
# 0b00000100 =   4: >1 pixel saturated
# 0b00001000 =   8: truncated / image boundary
# 0b00010000 =  16: aperture data incomplete/corrupted
# 0b00100000 =  32: isophotal data incomplete/corrupted
# 0b01000000 =  64: memory overflow during deblending
# 0b10000000 = 128: memory overflow during extraction
#
# The following flags define the FLAGS that exclude a source
#
# WCS can handle saturated/deblended/crowded sources
sexflag_wcs  = 0b11111000
#
# Photometric calibration needs perfect stars
sexflag_phot = 0b11111111



#
# Define some names for the columns in the source extractor catalog
#
SXapertures = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0] #arcsec
SXcolumn_names = [
    'ra', 'dec',
    'x', 'y', 
    'fwhm_image', 'fwhm_world',
    'background',
    'flags',
    'ota',
    'mag_aper_2.0',
    'mag_aper_3.0',
    'mag_aper_4.0',
    'mag_aper_5.0',
    'mag_aper_6.0',
    'mag_aper_8.0',
    'mag_aper_10.0',
    'mag_aper_12.0',
    'mag_err_2.0',
    'mag_err_3.0',
    'mag_err_4.0',
    'mag_err_5.0',
    'mag_err_6.0',
    'mag_err_8.0',
    'mag_err_10.0',
    'mag_err_12.0',
    'mag_auto',
    'mag_err_auto',
    'flux_max',
    'major_axis',
    'minor_axis',
    'position_angle',
    'elongation',
    'ellipticity',
    'source_id'
]
# Convert columns into dictionary to make look-up easier
SXcolumn = {}
for name in SXcolumn_names:
    SXcolumn[name] = len(SXcolumn)

SXcolumn_descriptions = [
    "ALPHAWIN_J2000      Right ascension of barycenter (J2000)                      [deg]",
    "DELTAWIN_J2000      Declination of barycenter (J2000)                          [deg]",
    "XWIN_IMAGE          Object position along x (double precision)                 [pixel]",
    "YWIN_IMAGE          Object position along y (double precision)                 [pixel]",
    "FWHM_IMAGE          FWHM assuming a gaussian core                              [pixel]",
    "FWHM_WORLD          FWHM assuming a gaussian core                              [pixel]",
    "BACKGROUND          Background at centroid position [count]",
    "FLAGS               Extraction flags",
    "EXT_NUMBER          FITS extension number",
    "MAG_APER_2.0        magnitude in aperture with radius 2.0''                    [mag]",
    "MAG_APER_3.0        magnitude in aperture with radius 3.0''                    [mag]",
    "MAG_APER_4.0        magnitude in aperture with radius 4.0''                    [mag]",
    "MAG_APER_5.0        magnitude in aperture with radius 5.0''                    [mag]",
    "MAG_APER_6.0        magnitude in aperture with radius 6.0''                    [mag]",
    "MAG_APER_8.0        magnitude in aperture with radius 8.0''                    [mag]",
    "MAG_APER_10.0       magnitude in aperture with radius 10.0''                   [mag]",
    "MAG_APER_12.0       magnitude in aperture with radius 12.0''                   [mag]",
    "MAG_ERR_2.0         magnitude uncertainty in aperture with radius 2.0''        [mag]",
    "MAG_ERR_3.0         magnitude uncertainty in aperture with radius 3.0''        [mag]",
    "MAG_ERR_4.0         magnitude uncertainty in aperture with radius 4.0''        [mag]",
    "MAG_ERR_5.0         magnitude uncertainty in aperture with radius 5.0''        [mag]",
    "MAG_ERR_6.0         magnitude uncertainty in aperture with radius 6.0''        [mag]",
    "MAG_ERR_8.0         magnitude uncertainty in aperture with radius 8.0''        [mag]",
    "MAG_ERR_10.0        magnitude uncertainty in aperture with radius 10.0''       [mag]",
    "MAG_ERR_12.0        magnitude uncertainty in aperture with radius 12.0''       [mag]",
    "MAG_AUTO            Kron-like elliptical aperture magnitude                    [mag]",
    "MAGERR_AUTO         RMS error for AUTO magnitude                               [mag]",
    "FLUX_MAX            Peak flux above background                                 [count]",
    "AWIN_IMAGE          Profile RMS along major axis                               [pixel]",
    "BWIN_IMAGE          Profile RMS along minor axis                               [pixel]",
    "THETA_IMAGE         Position angle (CCW/x)                                     [deg]",
    "ELONGATION          A_IMAGE/B_IMAGE",
    "ELLIPTICITY         1 - B_IMAGE/A_IMAGE",
]

   
reference_zeropoint = {
    "odi_g": [26.0, 26.2],
    "odi_r": [26.1, 26.3],
    "odi_i": [25.6, 25.8],
    "odi_z": [24.6, 24.8],
}
atm_extinction = {
    "odi_g": 0.2,
    "odi_r": 0.12,
    "odi_i": 0.058,
    "odi_z": 0.04,
}
photzp_colorterms = {
    'sdss_dr13': {
        "odi_g": [ [0.0000,  0.1600], 'g', 'r'],
        "odi_z": [ [0.0000, -0.1277], 'i', 'z'],
    },
    'sdss': {
        "odi_g": [ [0.0000,  0.1600], 'g', 'r'],
        "odi_z": [ [0.0000, -0.1277], 'i', 'z'],
    },
    'panstarrs': {
        # "odi_u": not available in PanSTARRS
        "odi_g": [ [0.0336, -0.0359], 'g', 'r'],
        # "odi_r":
        # "odi_i": No color term (|slope| < 0.01 ([ 0.00435847 -0.00343126])
        "odi_z": [ [0.0113, -0.0424], 'i', 'z'],
    }
}

#   "odi_r": [ 0.0047, 'g', 'i'],
#   "odi_r": [ 0.0074, 'g', 'r'],
#   "odi_r": [-0.0001, 'r', 'i'],
#   "odi_i": [-0.0027, 'r', 'i'],
#   "odi_i": [ 0.0022, 'i', 'z'],
#   "odi_i": [-0.0006, 'r', 'z'],


#
# Header keyword names in the TECHDATA extension
# to get the actual keywords, add the cell-id string, 
# i.e. "%02d%d%d" % (ota, cellx, celly)
#
techdata_keywords = [
    "GN__", "GN_E",
    "RN__", "RN_E",
    "RNE_", "RNEE", 
]
backup_gain = 1.3
backup_readnoise = 6.5
backup_readnoise_electrons = 9.0

invalid_sky_level_value = -99999.99


def get_valid_filter_name(hdr):
    """

    Convert the header FILTERID entry into a valid filtername. This is necessary
    as some of the filter names change over time, but the filter-id should
    remain unchanged.

    """

    try:
        filter_id = hdr['FILTERID'].strip()
        if (filter_id in list_of_valid_filter_ids):
            filter_name = list_of_valid_filter_ids[filter_id]
            return filter_name
    except:
        try: 
            filter = hdr['FILTER']
            return filter
        except:
            return "unknown"
    return 'unknown'





def get_cellmode(primhdr, cellhdr, focalplanelayout):
    """

    Check if the specified cell, identified by OTA and CELL-ID, is either broken
    or marked as a guide/video cell.

    """

    ota = int(primhdr['FPPOS'][2:4])

    ota_name = "OTA%02d" % ota
    extname = "OTA%02d.SCI" % ota

    cell = cellhdr['EXTNAME']

    # Check if this is one of the permanently broken cells
    wm_cellx, wm_celly = cellhdr['WN_CELLX'], cellhdr['WN_CELLY']
    broken = focalplanelayout.is_cell_broken(primhdr['OTA_ID'], wm_cellx, wm_celly)
    if (broken):
        #stdout_write ("Rejecting borken cell %s %s %s %s\n" %  (primhdr['OTA_ID'], ota, wm_cellx, wm_celly))
        return -1

    # It's not one of the broken cells, but it might still be a guide/video cell
    idx = wm_cellx + 8 * wm_celly
    cellmode = primhdr['CELLMODE']
    this_cellmode = cellmode[idx]
    cell_id = cellmode_ids[this_cellmode]

    return cell_id


def stdout_write(str):
    """

    Write a given text to stdout and flush the terminal. Pretty much what print
    does, but flushing the output afterwards.

    """

    sys.stdout.write(str)
    sys.stdout.flush()
    return

def clobberfile(filename):
    """

    Delete a file if it already exists, otherwise do nothing.

    """

    if (os.path.isfile(filename)):
        os.remove(filename)
    return



from types import *   
def shmem_as_ndarray( raw_array ):
    """

    Helper function needed to allocate shared memory numpy-arrays.

    """
    _ctypes_to_numpy = {
        ctypes.c_char : numpy.int8,
        ctypes.c_wchar : numpy.int16,
        ctypes.c_byte : numpy.int8,
        ctypes.c_ubyte : numpy.uint8,
        ctypes.c_short : numpy.int16,
        ctypes.c_ushort : numpy.uint16,
        ctypes.c_int : numpy.int32,
        ctypes.c_uint : numpy.int32,
        ctypes.c_long : numpy.int32,
        ctypes.c_ulong : numpy.int32,
        ctypes.c_float : numpy.float32,
        ctypes.c_double : numpy.float64
    }
    address = raw_array._wrapper.get_address()
    size = raw_array._wrapper.get_size()
    dtype = _ctypes_to_numpy[raw_array._type_]
    class Dummy(object): pass
    d = Dummy()
    d.__array_interface__ = {
         'data' : (address, False),
         'typestr' : ">f4", #FloatType, #"uint8", #numpy.uint8.str,
         'descr' : "", #"UINT8", #numpy.uint8.descr,
         'shape' : (size/4,),
         'strides' : None,
         'version' : 3
    }
    return numpy.asarray(d)#.view( dtype=numpy.float32 )




def sexa2deg(sexa):
    """

    Convert a sexa-decimal coordinate string (e.g. 12:30:00.0") into a float
    (12.5 in the example).

    """

    components = sexa.split(":")
    
    if (len(components) != 3):
        return 1e99
    
    deg = float(components[0])
    min = float(components[1])
    sec = float(components[2])
    
    return math.copysign(math.fabs(deg) + math.fabs(min/60.0) + math.fabs(sec/3600.0), deg)

def deg2sexa(deg, signed=False):
    """

    Convert a float coordinate into the more user-friednly sexa-decimal format

    """

    unsigned = math.fabs(deg)

    degrees = math.floor(unsigned)
    rest = (unsigned - degrees) * 60.0

    minutes = math.floor(rest)
    rest = (rest - minutes) * 60.0

    seconds = math.floor(rest)

    num = [math.copysign(degrees, deg), minutes, seconds]

    if (signed):
        text = "%+03d:%02d:%04.1f" % (int(math.copysign(degrees, deg)), int(minutes), seconds)
    else:
        text = "%02d:%02d:%04.1f" % (int(math.copysign(degrees, deg)), int(minutes), seconds)

    return text, num



headers_to_inherit = [
    'RA', 'DEC', 'TARGRA', 'TARGDEC', 'TELRAOFF', 'TELDECOF', 
    'FILTER', 'FILTERID', 'FILTDSCR', 'EXPTIME',
    'OBSID', 'OBJECT', 'OBSTYPE',
    'WCSASTRM',
    'EXPMEAS',
    
    'ORIGIN', 'INSTRUME',
    'FILENAME', 
    'OBSLOGIN',
    'RADESYS', 'TIMESYS', 'LSTHDR',
    'OBSERVAT', 'TELESCOP',
    'OBSERVER', 'PROPOSER', 'PROPID', 'PROGID', 'TACID', 'PROPPERD',
    'DATE-OBS', 'TIME-OBS', 'MJD-OBS', 'DATE',
    'ZD', 'AIRMASS',
    'TELFOCUS',
    'TRACK',
    'ELMAP', 'AZMAP', 'ROTPORT',
    'FOLDPOS','OBSBLOCK',
    'ADCMODE', 'ADCANG1', 'ADCANG2', 'ADCJD',
    'ROTSTART', 'ROTEND', 'ROTOFF', 
    
    'TEMPSTAT', 'DEWAR', 'COOLHEAD', 'COLPLATE', 'FOCPLATE', 'DEWPRESS',
    'FLTARM1A', 'FLTARM1B', 'FLTARM1C',
    'FLTARM2A', 'FLTARM2B', 'FLTARM2C',
    'FLTARM3A', 'FLTARM3B', 'FLTARM3C',
    'SHUTDIR', 'SHUTOPEN', 'SHUTCLOS',
    'CONTROLR',
    'IMAGESWV',
    ]

headers_to_delete_from_otas = [
    'CELLGAP1', 'CELLGAP2',
    'CNAXIS1', 'CNAXIS2',
    'NAMPS', 'NEXTEND', 
    'PRESCAN1', 'PRESCAN2',
    'OVRSCAN1', 'OVRSCAN2',
    'IMNAXIS1', 'IMNAXIS2',
    'EXTEND'
    ]


    


def inherit_headers(header, primary_header):
    """Copy all headers from the primary header into the current header

    """

    for header in headers_to_inherit:
        if (not header in primary_header):
            print("Problem with header ",header)
            continue

        card = primary_header.ascardlist()[header]
        header[card.key] = (card.value, card.comment)

        

def rebin_image(data, binfac, operation=numpy.mean):
    """

    Apply a binning factor to a data array.

    Parameters
    ----------
        
    data : ndarray

        Input data array. Only tested to work on two-dimensional data.

    binfac : int

        binning factor, e.g. 2 for 2xw binning. Only identical binning in both
        dimensions is supported at the present.

    operation : function (default: numpy.mean)

        What operation to use when combining the pixels. All functions operating
        on ndarrays are supported, but typical cases are likely one of the
        following:

        * numpy.mean
        * numpy.average
        * numpy.sum
        * numpy.median

        Or from the bottleneck package:

        * bottleneck.nanmean
        * bottleneck.nanmedian

    """
    if (binfac < 1):
        stdout_write("Rebinning at the moment only supports binning to larger pixels with binfac>1\n")
        return None
    elif (binfac == 1):
        return data

    out_size_x, out_size_y = int(math.ceil(data.shape[0]*1.0/binfac)), int(math.ceil(data.shape[1]*1.0/binfac))

    if (out_size_x*binfac != data.shape[0] or out_size_y*binfac != data.shape[1]):
        # The input array size is not a multiple of the new binning
        # Create a slightly larger array to hold the data to be rebinned
        container = numpy.zeros(shape=(out_size_x*binfac, out_size_y*binfac))

        # And insert the original data
        container[0:data.shape[0], 0:data.shape[1]] = data[:,:]
    else:
        container = data 
        
#    rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)

    rb1 = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
    rb2 = operation(rb1, axis=-1)
    rebinned = operation(rb2, axis=1)

#    rb1 = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
#    rb2 = nanmedian(rb1, axis=-1)
#    rebinned = nanmedian(rb2, axis=1)

#    #.nanmean(axis=-1).nanmean(axis=1)

    return rebinned



def center_coords(hdr):
    """

    Return the center coordinates of a given HDU, based on the WCS information
    (does not include distortion).

    """
    
    try:
        centerx, centery = hdr['NAXIS1']/2, hdr['NAXIS2']/2
    except:
        centerx, centery = 2048., 2048.

    center_ra  = (centerx-hdr['CRPIX1'])*hdr['CD1_1'] + (centery-hdr['CRPIX2'])*hdr['CD1_2'] + hdr['CRVAL1']
    center_dec = (centerx-hdr['CRPIX1'])*hdr['CD2_1'] + (centery-hdr['CRPIX2'])*hdr['CD2_2'] + hdr['CRVAL2']

    return center_ra, center_dec

    

def break_region_string(str_region):
    """

    Break down a IRAF-like string (e.g. [0:100,0:200]) into its four components
    and return them separately.

    """
    reg = str_region[1:-1]
    x,dummy,y = reg.partition(",")
    x1,dummy,x2 = x.partition(":")
    y1,dummy,y2 = y.partition(":")
    return int(x1)-1, int(x2)-1, int(y1)-1, int(y2)-1

def extract_region(data, str_region):
    """Extract a region based on the a IRAF-like region definition.

    See also
    --------
    `break_region_string`

    """
    x1,x2,y1,y2 = break_region_string(str_region)
    return data[y1:y2+1, x1:x2+1]


def insert_into_array(data, from_region, target, target_region):
    """

    Copy data from one array into another array, with source and target
    coordinates specified in the IRAF-style format.

    """
    fx1, fx2, fy1, fy2 = break_region_string(from_region)
    tx1, tx2, ty1, ty2 = break_region_string(target_region)

    if (fx2-fx1 != tx2-tx1 or fy2-fy1 != ty2-ty1):
        print("Dimensions do not match, doing nothing")
    else:
        target[ty1:ty2+1, tx1:tx2+1] = data[fy1:fy2+1, fx1:fx2+1]

    return 0

def mask_broken_regions(datablock, regionfile, verbose=False, reduction_log=None):
    """

    Mask out broken regions in a data array. Regions are defined via a ds9
    region file to allow for ay creation by the user.

    """

    logger = logging.getLogger("MaskBrokenRegion")

    if (not os.path.isfile(regionfile)):
        if (not reduction_log == None):
            reduction_log.failed('badpixels')
        return datablock

    counter = 0
    file = open(regionfile)
    for line in file:
        if (line[0:3] == "box"):
            coords = line[4:-2]
            coord_list = coords.split(",")
                        
            if (not type(datablock) == type(None)):
                x, y = int(float(coord_list[0])), int(float(coord_list[1]))
                dx, dy = int(0.5*float(coord_list[2])), int(0.5*float(coord_list[3]))
                #mask[y-dy:y+dy,x-dx:x+dx] = 1

                x1 = numpy.max([0, x-dx])
                x2 = numpy.min([datablock.shape[1], x+dx])
                y1 = numpy.max([0, y-dy])
                y2 = numpy.min([datablock.shape[0], y+dy])
                datablock[y1:y2, x1:x2] = numpy.NaN
                logger.debug("Masking block X=%d-%d, y=%d-%d" % (x1,x2,y1,y2))

                # print x,x+dx,y,y+dy
            counter += 1

        if (line[0:8] == "# vector"):
            coords = line[9:]
            end = coords.find(")")
            coords = coords[:end]
            coord_list = coords.split(",")

            if (not type(datablock) == type(None)):
                x, y = int(float(coord_list[0])), int(float(coord_list[1]))
                angle = int(float(coord_list[3]))

                for cx,cy in itertools.product(range(8), repeat=2):
                    x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy) 
                    if (x1 <= x <= x2) and (y1 <= y <=y2):

                        horiz = False
                        vert = False
                   
                        if (357 <= angle or angle <= 3) or (177 <= angle <= 183):
                            horiz = True

                        if (87 <= angle <= 93) or (267 <= angle <= 273):
                            vert = True

                        if (horiz): 
                            datablock[y, x1:x2] = numpy.NaN
                            logger.debug("Masking line x=%d-%d, y=%d" % (x1,x2,y))

                        if (vert): 
                            datablock[y1:y2, x] = numpy.NaN
                            logger.debug("Masking column x=%d, y=%d-%d" % (x,y1,y2))
                            #print line,"vertical",x,y
                           
                        if (not horiz and not vert):
                            logger.warning("vector line/column region ambiguous in %s\n>>> %s" % (
                                regionfile, line))
                            #print "vector line/column region ambiguous"

                        break
     
                counter += 1

    file.close()

    logger.debug("Marked %d bad pixel regions" % (counter))

    if (not reduction_log == None):
        reduction_log.success('badpixels')

    return datablock


def is_image_extension(hdu):
    """

    Check if a given HDU is a Image extension

    """
    if (hdu == None):
        return False

    if (type(hdu) == pyfits.hdu.image.ImageHDU or
        type(hdu) == pyfits.hdu.compressed.CompImageHDU):

        try:
            if (type(hdu.data) == type(None)):
                return False
        except:
            # Can't access the .data block, hence this can't be an image extension
            return False

        return True

    elif (type(hdu) == pyfits.hdu.image.PrimaryHDU):
        # This might or might not be an image
        if ('NAXIS' in hdu.header and
            'NAXIS1' in hdu.header and
            'NAXIS2' in hdu.header and
            hdu.header['NAXIS'] == 2 and
            hdu.header['NAXIS1'] > 0 and
            hdu.header['NAXIS2'] > 0):
            # This looks like it might be an image
            return True
            
    return False



def get_svn_version():
    """

    Return current SVN version as given by the `svnversion` command.

    """

    try:
        p = subprocess.Popen('svnversion -n', shell=True, stdout=subprocess.PIPE)
        svn_version, err = p.communicate()
        ret = p.wait()
        if (ret != 0):
            svn_version = "problem_with_svn"
    except:
        svn_version="no_svnversion_found"

    return svn_version

def log_svn_version(hdr):
    """

    Add SVN version number to FITS header.

    """
    svn = get_svn_version()
    hdr["QPIPESVN"] = (svn, "QuickReduce Revision")
    return



def rotate_around_center(data, angle, mask_limit = 0.1, verbose=True, safety=1, mask_nans=True, spline_order=3):
    """

    Rotate a given data array. Rotation center is the center of the data frame.

    """
    logger = logging.getLogger("RotateFrame")
    logger.debug("Rotating data block by %.1f deg ..." % (angle))

    # Prepare mask so we can mask out undefined regions
    mask = numpy.zeros(shape=data.shape)
    mask[numpy.isnan(data)] = 1.

    # Replace NaN's with some numer to make interpolation work
    data[numpy.isnan(data)] = 0.


    # Now rotate both the image and its mask
    rotated = scipy.ndimage.interpolation.rotate(input=data, angle=angle, axes=(1, 0), 
                                                 reshape=False, order=spline_order,
                                                 mode='constant', cval=0, )

    if (mask_nans):
        logger.debug("masking NaNs")
        rotated_mask = scipy.ndimage.interpolation.rotate(input=mask, angle=angle, axes=(1, 0), 
                                                          reshape=False, order=1,
                                                          mode='constant', cval=1, )

        # Blur out the NaN mask to make sure we get clean edges. This approach
        # is on the conservative side, rather clipping some pixels too many then to 
        # add artificial edges that later are hard to remove and/or deal with
        filter_gaussed = scipy.ndimage.filters.gaussian_filter(input=rotated_mask, order=0, sigma=safety)

        # And finally apply the mask
        # rotated[rotated_mask > mask_limit] = numpy.NaN
        rotated[filter_gaussed > mask_limit] = numpy.NaN
        
    # and return the results
    logger.debug("done!")
    return rotated



def get_filter_level(header):
    """

    Return the level of the installed filter, based on the information in the
    FLTARM header keywords. If more than one filter is active, this function
    returns the level of the lowest populated filter arm.

    """
    filter_level = 0
    filter_count = 0
 
    for lvl in range(1,4):
        for arm in "ABC":
            keyword = "FLTARM%d%s" % (lvl, arm)
            # print keyword

            if (header[keyword].strip() == "IN"):
                filter_count += 1
                if (filter_count == 1):
                    filter_level = lvl

    if (filter_count > 1):
        return -1

    return filter_level
            




def cell2ota__extract_data_from_cell(data_in=None):
    """
    Don't use anymore!
    """
    if (data_in == None):
        return numpy.zeros(shape=(494,480))

    return data_in[0:494, 0:480]

def cell2ota__get_target_region(x, y, binning=1, trimcell=None):
    """

    Get the location of a given cell in the monolithic OTA array, accounting for
    binning (only 1x1 and 2x2 supported).

    """
    #taken from ODI OTA Technical Datasheet (det area 480x494, streets 11/28 px)

    # Y-coordinates of cell numbers are counted top down,
    # but pixel coordinates are counted bottom up
    _y = 7 - y

    if (trimcell == None):
        trimcell = 0

    if (binning == 1):
        y1 = (505*_y)  #was 503 
        y2 = y1 + 494
        x1 = 508 * x
        x2 = x1 + 480
    elif (binning == 2):
        y1 = int(math.ceil(252.5 * _y))
        y2 = y1 + 247
        x1 = 254 * x
        x2 = x1 + 240
        
    return x1+trimcell, x2-trimcell, y1+trimcell, y2-trimcell




def three_sigma_clip(input, ranges=[-1e9,1e9], nrep=3, return_mask=False, nsigma=3):
    """

    Perfom an iterative 3-sigma clipping on the passed data array. 

    """

    old_valid = numpy.isfinite(input)
    valid = (input > ranges[0]) & (input < ranges[1])
                
    for rep in range(nrep):
        if (numpy.sum(valid) < 1):
            valid = old_valid
            break

        lsig = scipy.stats.scoreatpercentile(input[valid], 16)
        hsig = scipy.stats.scoreatpercentile(input[valid], 84)
        median = numpy.median(input[valid])
        sigma = 0.5 * (hsig - lsig)

        mingood = numpy.max([median - nsigma*sigma, ranges[0]])
        maxgood = numpy.min([median + nsigma*sigma, ranges[1]])

            #print median, sigma
        old_valid = valid
        valid = (input > mingood) & (input < maxgood)
        
    if (return_mask):
        return input[valid], valid

    return input[valid]


def derive_ota_outlines(otalist):
    """

    For each OTA (extension) in the pased list, derive the sky-position of all 4
    corners.

    """
    from astLib import astWCS

    all_corners = []
    for ext in range(len(otalist)):
        if (not is_image_extension(otalist[ext])):
            continue
#type(otalist[ext]) == pyfits.hdu.image.ImageHDU):
            
        wcs = astWCS.WCS(otalist[ext].header, mode='pyfits')
            
        corner_coords = []
        corner_coords.append(wcs.pix2wcs(                            0,                             0))
        corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'],                             0))
        corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'], otalist[ext].header['NAXIS2']))
        corner_coords.append(wcs.pix2wcs(                            0, otalist[ext].header['NAXIS2']))

        all_corners.append(corner_coords)

    return all_corners


def create_qa_filename(outputfile, plotname, options):
    """

    Return the filename for a given diagnostic plot, accounting for
    user-specified preferences.

    """
    if (options['structure_qa_subdirs']):
        dirname, basename = os.path.split(outputfile)
        if (dirname == None or dirname == ''): 
            dirname = "."
        qa_plotdir = "%s/%s/" % (dirname, options['structure_qa_subdir_name'])
        if (not os.path.isdir(qa_plotdir)):
            os.mkdir(qa_plotdir)
        qa_plotfile = "%s/%s" % (qa_plotdir, plotname)
    else:
        qa_plotfile = "%s.%s" % (outputfile[:-5], plotname)

    return qa_plotfile

def create_qa_otaplot_filename(plotname, ota, structure_qa_subdirs):
    """

    Return the filename for a given OTA-level diagnostic plot, accounting for
    user-specified preferences.

    """

    if (structure_qa_subdirs):
        # in this case, plotname is to be taken as directory
        if (not os.path.isdir(plotname)):
            os.mkdir(plotname)
        # The actual filenames are the directory and the OTA
        qa_plotfile = "%s/OTA%02d" % (plotname, ota)
    else:
        qa_plotfile = "%s_OTA%02d" % (plotname, ota)

    return qa_plotfile








#
# Additional functions to handle binned data more elegantly
#
def get_binning(hdr):
    """

    Get the binning factor of a given frame/cell based on its header.

    """
    # print "determining binning factor"
    
    # The CCDBIN1 keyword should be foung in primary header, so try this one first
    try:
        binfac = hdr['CCDBIN1']
        # print binfac
        return int(binfac)
    except:
        pass

    # Couldn't find the CCDBIN1 header, so this is likely not a primary HDU
    # Try the CCDSUM keyword next
    try:
        ccdsum = hdr['CCDSUM']
        # print ccdsum
        items = ccdsum.split()
        return int(items[0])
    except:
        pass

    #
    # If we still don't find a valid header, look at the size of the data section
    #

    # Still no luck finding a header? Assume it's 1 and continue
    return 1


def get_collected_image_dimensions(binning):
    """

    Return the dimension of the monolithic OTA frame, accounting for binning

    """
    if (binning == 1):
        sizex, sizey = 4096, 4096
    elif (binning == 2):
        sizex, sizey = 2048, 2048

    return sizex, sizey


def extract_datasec_from_cell(data, binning, trimcell=None):
    """

    Return the science region of a given cell, accounting for binning

    """

    if (trimcell == None):
        trimcell=0

    if (binning == 1):
        dx, dy = 480, 494
    elif (binning == 2):
        dx, dy = 240, 247

    # print "extracting datasec", dx, dy
    return data[0+trimcell:dy-trimcell, 0+trimcell:dx-trimcell]


def extract_biassec_from_cell(data, binning):
    """

    Return the overscan region of a given cell, accounting for binning

    """

    if (binning == 1):
        dx1, dx2, dy1, dy2 = 500, 550, 0, 494
    elif (binning == 2):
        dx1, dx2, dy1, dy2 = 260, 277, 0, 246

    # print dx1, dx2, dy1, dy2
    return data[dy1:dy2, dx1:dx2]


def match_catalog_areas(src, to_match, radius):
    """

    Match the area coverage of two catalogs, allowing for some extra coverage
    around the edges.

    """

    #
    # Account for the cos(dec) effect
    #
    max_dec = numpy.max(numpy.fabs(src[:,1]))
    if (max_dec > 85): max_dec = 85
    cos_dec = numpy.cos(numpy.radians(max_dec))

    #
    # Apply cos(dec) correction to the src catalog 
    # (make sure to copy the catalogs to prevent altering the 
    # actual input catalogs) ...
    #
    src_radec = src[:,0:2].copy()
    src_radec[:,0] *= cos_dec
    #src_radec[:,0] *= numpy.cos(numpy.radians(src_radec[:,1]))
    src_tree = scipy.spatial.cKDTree(src_radec)

    #
    # ... and to the match catalog
    #
    match_radec = to_match[:,0:2].copy()
    match_radec[:,0] *= cos_dec

    #
    # Now determine, for each match source (in 2MASS), the distance 
    # to the 1 (k=1) closest source in the src (ODI) catalog. If the distance is
    # less than 'radius', the source is close enough and is considered matched
    #
    d,l = src_tree.query(match_radec, k=1, eps=0.05, p=2, distance_upper_bound=radius)
    matched = to_match[d<radius]

    return matched






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
                print("% 15s : %s" % (key, filename))

    # print reduction_step
    # print short_filename

    columns = [\
        pyfits.Column(name='correction',    format='A25',  array=reduction_step),
        pyfits.Column(name='filename_full', format='A375', array=full_filename),
        pyfits.Column(name='filename',      format='A100', array=short_filename),
        ]

    coldefs = pyfits.ColDefs(columns)
    try:
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    except:
        tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.name = "ASSOCIATIONS"

    return tbhdu



def add_fits_header_title(header, title, before_keyword):
    header.add_blank(before=before_keyword)
    header.add_blank(title, before=before_keyword)
    header.add_blank(before=before_keyword)
    return

def sigma_to_percentile(s):
    return 100 * 0.5*(1 + scipy.special.erf(s/math.sqrt(2)))




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
            objectname = header['OBJECT']
            strings_to_replace = [',', '(', ')', '/', '\\', '`', '"', '\'', '{', '}', '[', ']', '&',
                                      ' ', '*', '?']
            for _s in strings_to_replace:
                objectname = objectname.replace(_s, '_')
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



def read_wipecells_list():

    logger = logging.getLogger("ReadWipeCells")

    if (not cmdline_arg_isset("-wipecells")):
        return None

    wipecells = {}
    wc = get_cmdline_arg("-wipecells")
    logger.debug("wipecells: %s" % (wc))
    for _wc in wc.split(","):
        # print _wc
        _ota,_cell = _wc.split(".")
        ota = int(_ota)
        ota_name = "OTA%02d.SCI" % (ota)
        cellx,celly = int(_cell[0]), int(_cell[1])
        if (not ota_name in wipecells):
            wipecells[ota_name] = []
        wipecells[ota_name].append((cellx, celly))
    logger.debug("wipe-cells final: %s" % (str(wipecells)))

    return wipecells

def wipecells(ext, wipecells_list, binning=1, fillvalue=numpy.NaN):

    logger = logging.getLogger("WipeCells")
    
    if (not is_image_extension(ext) or
        not ext.name in wipecells_list):
        return

    # We have some cells to wipe
    for (cell_x, cell_y) in wipecells_list[ext.name]:
        # Get coordinates for this cell, for the given binning
        logger.debug("Wiping out cell %d,%d in OTA %s" % (cell_x, cell_y, ext.name))
        x1, x2, y1, y2 = cell2ota__get_target_region(cell_x, cell_y, binning=binning)
        #_was = bottleneck.nanmedian(ext.data[y1:y2, x1:x2].astype(numpy.float32))
        ext.data[y1:y2, x1:x2] = fillvalue
        #logger.info("%d:%d, %d:%d --> %f --> %f"  %(
        #    x1,x2,y1,y2, _was, bottleneck.nanmedian(ext.data[y1:y2, x1:x2].astype(numpy.float32))))

    return

def is_guide_ota(primhdu, ext, w=20):

    logger = logging.getLogger("IsGuideOTA")

    binning = primhdu.header['BINNING']
    skylevel = primhdu.header['SKYLEVEL']
    gain = primhdu.header['GAIN']
    skynoise = primhdu.header['SKYNOISE']

    logger.debug("Checking OTA %s (bin=%d, sky=%.1f, skynoise=%.2f)" % (
        ext.name, binning, skylevel, skynoise))

    if (not is_image_extension(ext)):
        logger.debug("extension is not a valid image extension")
        return False

    excesses = numpy.empty((8,8))
    excesses[:,:] = numpy.NaN

    for cx, cy in itertools.product(range(8), repeat=2):

        #
        # Get pixel coord for this cell
        #
        
        x1,x2,y1,y2 = cell2ota__get_target_region(cx, cy, binning=binning, trimcell=0)
        x21 = (x2-x1)//2

        # extract the mean value in the bottom corner
        corner = bottleneck.nanmean(ext.data[y1:y1+w, x1:x1+w].astype(numpy.float32))

        # also get the value in the bottom center
        center = bottleneck.nanmean(ext.data[y1:y1+w, x1+x21-w//2:x1+x21+w//2].astype(numpy.float32))

        excess = corner - center
        #print ext.name, cx, cy, corner, center, excess
            
        excesses[cx,cy] = excess

    _mean = bottleneck.nanmean(excesses)
    _median = bottleneck.nanmedian(excesses)

    is_guideota = (_median > 10*skynoise)
    logger.debug("Found corner excess mean=%.1f, median=%.1f --> guide-OTA: %s" % (
        _mean, _median, "YES" if is_guideota else "NO"))

    return is_guideota
