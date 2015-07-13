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

import os
import sys
import logging
import numpy
import pyfits

def cmdline_arg_isset(arg):
    """

    Check if the given argument was given on the command line

    """
    # Go through all command line arguments and check
    # if the requested argument is one of them
    for cur in sys.argv[1:]:
        name,sep,value = cur.partition("=")
        if (name == arg):
            return True
    return False


def get_cmdline_arg(arg):
    """

    Retreive the value of a command-line argument

    """

    # Check all arguments if 
    for cur in sys.argv[1:]:
        name,sep,value = cur.partition("=")
        if (name == arg):
            if (sep == ""):
                return None
            return value
    return None


def get_clean_cmdline():
    """
    
    Return all values entered on the command line with all command-line flags
    removed.

    Example:
        
        The user called a program ``./podi_something param1 -flag1 -flag2 param2``

        This function would then return a list containing
        ['./podi_something', 'param1', 'param2']

    """
    
    list = []
    for cur in sys.argv:
        if (cur[0] != "-" or cur[1].isdigit()):
            list.append(cur)
    return list


def cmdline_arg_set_or_default(name, defvalue):
    """

    Return the value of a command line argument. If no argument was passed (for
    example -flag= ), assign the specified default value instead.

    """

    if (cmdline_arg_isset(name)):
        val = get_cmdline_arg(name)
        if (val == None):
            return defvalue
        return val
    return defvalue


import podi_sitesetup as sitesetup
from podi_definitions import *


def read_comma_separated_list(inp):

    print "INPINPINP",inp
    ret_list = []
    if (type(inp) == list):
        return inp
    elif (not inp == None):
        for _in in inp.split(","):
            if (os.path.isdir(_in) or os.path.isfile(_in)):
                ret_list.append(_in)
    print ret_list
    return ret_list


    
def read_options_from_commandline(options=None, ignore_errors=False):
    """
    Read all command line options and store them in the options dictionary.

    """

    logger = logging.getLogger("ReadOptions")

    if (options == None):
        options = set_default_options()

    options['verbose'] = cmdline_arg_isset("-verbose")

    # Handle all reduction flags from command line
    if (cmdline_arg_isset("-cals")):
        cals_dir = read_comma_separated_list(get_cmdline_arg("-cals"))
        if (cals_dir == []):
            logger.critical("The specified cals-directory (%s) does not exist!!!" % (
                get_cmdline_arg("-cals")))
            if (not ignore_errors):
                sys.exit(0)

        options['bias_dir'] = cals_dir
        options['dark_dir'] = cals_dir
        options['flat_dir'] = cals_dir
        options['techdata'] = cals_dir
        logger.debug("tech-data: %s" % (options['techdata']))

        options['illumcorr_dir'] = cals_dir

    options['bias_dir'] = read_comma_separated_list(
        cmdline_arg_set_or_default("-bias", options['bias_dir']))
    options['dark_dir'] = read_comma_separated_list(
        cmdline_arg_set_or_default("-dark", options['dark_dir']))
    options['flat_dir'] = read_comma_separated_list(
        cmdline_arg_set_or_default("-flat", options['flat_dir']))

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

    options["update_persistency_only"] = cmdline_arg_isset("-update_persistency_only")

    options['fringe_dir'] = cmdline_arg_set_or_default('-fringe', "auto")
    options['fringe_vectors'] = cmdline_arg_set_or_default("-fringevectors", options['fringe_vectors'])

    options['pupilghost_dir'] = cmdline_arg_set_or_default('-pupilghost', None)

    options['fixwcs'] = cmdline_arg_isset("-fixwcs")

    # For now assume that the WCS template file is located in the same directory as the executable
    options['wcs_distortion'] = sitesetup.exec_dir + "/"
    options['wcs_distortion'] = cmdline_arg_set_or_default("-wcs", options['wcs_distortion'])
    if (not os.path.isfile(options['wcs_distortion']) and 
        not os.path.isdir(options['wcs_distortion'])):
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
        logger.info("Applying a user-defined WCS offset of %.3f, %.3f arcseconds\n" % (options['offset_pointing'][0], options['offset_pointing'][1]))

    #
    # Read all offsets from command line
    # For convenience, there are two sets of offset parameters, that internally simply 
    # get added up. The reason for this is to make specifying them on the command line 
    # easier, since the pointing offsets stay constant across a dither pattern, while 
    # the dither offsets change.
    #
    options['target_coords'] = None
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
        logger.info("writing plots as %s" % (str(options['plotformat'])))
        
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
            comma_pos = key_value.find(",")
            if (comma_pos < 0):
                # no fitsheader given, ignore this one
                continue
            elif (comma_pos > 8):
                logger.warning("FITS keyword given in -addfitskey option (%s) exceeds 8 character limit" % (
                    key_value[:comma_pos]))
                continue
            else:
                key = key_value[:comma_pos]
                value = key_value[comma_pos+1:]
            logger.info("Adding fits keyword %s = %s" % (key, value))

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
                ns['ref_obsid'] = None
                try:
                    ns['ref_mjd'] = float(ns['ref'])
                    logger.debug("Found reference MJD (%f) on command line" % (ns['ref_mjd']))
                except:
                    if (os.path.isfile(ns['ref'])):
                        hdulist = pyfits.open(ns['ref'])
                        if ("MJD-OBS" in hdulist[0].header):
                            ns['ref_mjd'] = hdulist[0].header['MJD-OBS'] * 1.0
                            ns['ref_obsid'] = hdulist[0].header['OBSID']
                            logger.debug("Found reference MJD (%f) in file %s" % (ns['ref_mjd'], ns['ref']))
                        hdulist.close()
                        del hdulist
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

    if (not cmdline_arg_isset('-illumcorr')):
        options['illumcorr_dir'] = None
    else:
        options['illumcorr_dir'] = cmdline_arg_set_or_default("-illumcorr", options['illumcorr_dir'])
        logger.debug("Using illumination correction from %s" % (options['illumcorr_dir']))

    if (cmdline_arg_isset("-keepsex")):
        sitesetup.sex_delete_tmps = False

    options['prestage'] = cmdline_arg_isset("-prestage")

    options['softbin'] = int(cmdline_arg_set_or_default("-softbin", 0))

    if (cmdline_arg_isset("-selectota")):
        str_otas = cmdline_arg_set_or_default("-selectota", "33")
        otas = str_otas.split(",")
        options['selectota'] = [None] * len(otas)
        for idx, ota in enumerate(otas):
            x = int(ota[0])
            y = int(ota[1])
            options['selectota'][idx] = (x,y)
        
#    options['selectota'] = int(cmdline_arg_set_or_default("-selectota", None))

    return options









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
    options['exec_dir'] = sitesetup.exec_dir

    options['update_persistency_only'] = False
    options['persistency_dir'] = None
    options["persistency_map"] = None
    options['max_persistency_time'] = sitesetup.persistency_duration

    options['fringe_dir'] = None
    options['fringe_vectors'] = None #"%s/.fringevectors/" % (options['exec_dir'])

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

    options['illumcorr_dir'] = None

    options['prestage'] = False

    options['softbin'] = 0
    options['selectota'] = None

    return options


