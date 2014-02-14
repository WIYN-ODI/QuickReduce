#!/usr/bin/env python
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
into the filnal step. While this is overkill for stellar fields, this approach
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


import os
import sys
import pyfits
import subprocess

from podi_collectcells import *


swarp_exec = "swarp"
tmp_dir = "/scratch/"
single_dir = "."

import podi_logging
import logging



def swarpstack():

    logger = logging.getLogger("SwarpStack - Prepare")

    # Figure out the config path
    abspath = os.path.abspath(sys.argv[0])
    dirname, filename = os.path.split(abspath)
    swarp_default = "%s/.config/swarp.default" % (dirname)

    stacked_output = get_clean_cmdline()[1]

    inputfiles = get_clean_cmdline()[2:]

    pixelscale = float(cmdline_arg_set_or_default("-pixelscale", 0))
    subtract_back = cmdline_arg_isset("-bgsub")

    reuse_singles = cmdline_arg_isset("-reusesingles")

    use_nonsidereal = cmdline_arg_isset("-nonsidereal")
    options = read_options_from_commandline()
    if (use_nonsidereal):
        logger.info("Applying the non-sidereal motion correction")

        # First apply the non-sidereal correction to all input frames
        # Keep frames that are already modified
        for i in range(len(inputfiles)):
            hdulist = pyfits.open(inputfiles[i])
            # Assemble the temporary filename for the corrected frame
            
            corrected_filename = "%(single_dir)s/%(obsid)s.nonsidereal.fits" % {
                "single_dir": single_dir,
                "obsid": hdulist[0].header['OBSID'],
            }

            # Check if the corrected file already exists - if not create it
            #if (not os.path.isfile(corrected_filename)):
            logger.info("Correcting input frame %s --> %s" % (
                inputfiles[i], corrected_filename))
            # Apply the non-sidereal option
            apply_nonsidereal_correction(hdulist, options, logger)
            hdulist.writeto(corrected_filename, clobber=True)
            # else:
            #     logger.info("Input-frame %s with correction already exists (%s)" % (
            #         inputfiles[i], corrected_filename))

            # Now change the filename of the input list to reflect 
            # the corrected file
            inputfiles[i] = corrected_filename
        #
        # By now all frames have the non-sidereal correction applied,
        # so we can go ahead and stack them as usual
        #

    target_dimension = float(cmdline_arg_set_or_default('-dimension', -1))

    add_only = cmdline_arg_isset("-add") and os.path.isfile(stacked_output)
    if (add_only):
        logger.info("Activating ADD mode")

    if (stacked_output.endswith(".fits")):
        stacked_output = stacked_output[:-5]

    header_only_file = "%s/preswarp.fits" % (tmp_dir)

    reference_file = cmdline_arg_set_or_default("-reference", None)
    # Make sure the reference file is a valid file
    if (not reference_file == None and os.path.isfile(reference_file)):
        logger.info("Using %s as reference file" % (reference_file))
    else:
        reference_file = None

    if (add_only or not reference_file == None):
        #
        # This is the simpler add-only mode
        #

        if (not reference_file == None):
            output_info = pyfits.open(reference_file)
        else:
            # Open the existing output header and get data from there
            output_info = pyfits.open(stacked_output+".fits")

        logger.info("Stack information...")
        logger.info("   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header))

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']

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
                  'weightout': "%s/preswarp.weight.fits" % (tmp_dir),
                  'combine_type': 'AVERAGE',
              }

        if (pixelscale > 0):
            swarp_opts += " -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.4f " % (pixelscale)

        swarp_opts += " -SUBTRACT_BACK %s " % ("Y" if subtract_back else "N")

        logger.debug("SWARP options for pre-stack:\n"+" ".join(swarp_opts.split()))

        # 
        # First create only the output header so we can pass some information 
        # to the user
        #
        swarp_cmd = "%(swarp)s %(opts)s -HEADER_ONLY Y %(files)s" % {
            'swarp': swarp_exec,
            'opts': swarp_opts,
            'files': " ".join(inputfiles),
        }
        logger.debug("swarp_cmd=\n"+swarp_cmd)


        try:
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
            logger.debug("swarp stdout:\n"+swarp_stdout)
            if (len(swarp_stderr) > 0 and ret.returncode != 0):
                logger.warning("swarp stderr:\n"+swarp_stderr)
            else:
                logger.debug("swarp stderr:\n"+swarp_stderr)
        except OSError as e:
            podi_logging.log_exception()
            print >>sys.stderr, "Execution failed:", e

        #
        # some information about the resulting stack is in the output-file
        #

        output_info = pyfits.open(header_only_file)
        print "Stack information..."
        print "   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header)

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']
        

    
    #
    # Prepare the individual frames, rectified and re-projected 
    # to the final grid
    #
    logger = logging.getLogger("SwarpStack - Singles")
    single_prepared_files = []
    for singlefile in inputfiles:
        hdulist = pyfits.open(singlefile)
        obsid = hdulist[0].header['OBSID']

        # assemble all swarp options for that run
        dic = {'singledir': single_dir,
               'obsid': obsid,
               'pixelscale': pixelscale,
               'pixelscale_type': "MANUAL" if pixelscale > 0 else "MEDIAN",
               'center_ra': out_crval1,
               'center_dec': out_crval2,
               'imgsizex': out_naxis1,
               'imgsizey': out_naxis2,
               'resample_dir': tmp_dir,
               'inputfile': singlefile,
               'swarp_default': swarp_default,
           }

        swarp_opts = """\
                 -c %(swarp_default)s \
                 -IMAGEOUT_NAME %(singledir)s/%(obsid)s.fits \
                 -WEIGHTOUT_NAME %(singledir)s/%(obsid)s.weight.fits \
                 -PIXEL_SCALE %(pixelscale)f \
                 -PIXELSCALE_TYPE %(pixelscale_type)s \
                 -COMBINE Y \
                 -COMBINE_TYPE AVERAGE \
                 -CENTER_TYPE MANUAL \
                 -CENTER %(center_ra)f,%(center_dec)f \
                 -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
                 -RESAMPLE_DIR %(resample_dir)s \
                 -SUBTRACT_BACK N \
                 %(inputfile)s \
                 """ % dic

        single_file = "%(singledir)s/%(obsid)s.fits" % dic

        # print swarp_opts
        swarp_cmd = "%s %s" % (swarp_exec, swarp_opts)

        if (add_only and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, skipping it" % (single_file))
        elif (reuse_singles and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, re-using it" % (single_file))
            single_prepared_files.append(single_file)
        else:
            logger.info("Preparing file %s, please wait ..." % (singlefile))
            logger.debug(" ".join(swarp_cmd.split()))
            try:
                ret = subprocess.Popen(swarp_cmd.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
                (swarp_stdout, swarp_stderr) = ret.communicate()
                logger.debug("swarp stdout:\n"+swarp_stdout)
                if (len(swarp_stderr) > 0 and ret.returncode != 0):
                    logger.warning("swarp stderr:\n"+swarp_stderr)
                else:
                    logger.debug("swarp stderr:\n"+swarp_stderr)
                
                #print "\n".join(swarp_stderr)
                single_prepared_files.append(single_file)
            except OSError as e:
                podi_logging.log_exception()
                print >>sys.stderr, "Execution failed:", e


    #
    # If in "add" mode, rename the previous output file and add it to the list of input files
    #
    if (add_only):

        if (len(single_prepared_files) < 1):
            logger.info("No new files were added, so there's nothing to do.")
            return

        prev = 1
        while (True):
            filename = "%s.prev%02d.fits" % (stacked_output, prev)
            if (not os.path.isfile(filename)):
                break
            prev += 1
            continue
                
        # Rename the current output file and its weights
        old_stacked = "%s.prev%02d.fits" % (stacked_output, prev)
        old_weight = "%s.prev%02d.weight.fits" % (stacked_output, prev)

        os.rename(stacked_output+".fits", old_stacked)
        logger.debug("renamed old stack %s -> %s" % (stacked_output+".fits", old_stacked))

        os.rename(stacked_output+".weight.fits", old_weight)
        logger.debug("renamed old stack weight %s -> %s" % (stacked_output+".weight.fits", old_weight))

        # Also add the new re-named old stacked file to list of input files
        single_prepared_files.append(old_stacked)
        logger.debug("Adding now old stack file to input list")
    #
    # Done re-naming the old file
    #


    #
    # Now all single files are prepared, go ahead and produce the actual stack
    #
    dic['combine_type'] = "AVERAGE"
    dic['imageout'] = stacked_output+".fits"
    dic['weightout'] = stacked_output+".weight.fits"
    dic['prepared_files'] = " ".join(single_prepared_files)
    dic['bgsub'] = "Y" if subtract_back else "N"

    swarp_opts = """\
                 -c %(swarp_default)s \
                 -IMAGEOUT_NAME %(imageout)s \
                 -WEIGHTOUT_NAME %(weightout)s \
                 -COMBINE_TYPE %(combine_type)s \
                 -PIXEL_SCALE %(pixelscale)f \
                 -PIXELSCALE_TYPE %(pixelscale_type)s \
                 -COMBINE Y \
                 -COMBINE_TYPE %(combine_type)s \
                 -CENTER_TYPE MANUAL \
                 -CENTER %(center_ra)f,%(center_dec)f \
                 -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
                 -RESAMPLE_DIR %(singledir)s \
                 -SUBTRACT_BACK %(bgsub)s \
                 -WEIGHT_TYPE MAP_WEIGHT \
                 """ % dic

#                 -RESAMPLE N

    #
    # Now use some brains to figure out the best way of setting the background 
    # subtraction to get s nice smooth background that does not over-subtract 
    # the target.
    #
    if (target_dimension > 0 and subtract_back):
        dic['BACK_TYPE'] = "AUTO"
        dic['BACK_SIZE'] = 128
        dic['BACK_FILTERSIZE'] = 3

        # Rule of thum: larger objects: make both filtersize and back_size 
        # larger
        # first compute a reference size for the default settings
        ref_size = dic['BACK_SIZE'] * dic['BACK_FILTERSIZE'] \
                   * pixelscale / 60. \
                   * 0.1  # the last factor is a fudge-factor
        logger.debug("Reference size: %f" % (ref_size))
        # Now scale up the filtersize, making sure it stays between 3 and 7
        filtersize = int(math.floor(math.sqrt(target_dimension / ref_size) * dic['BACK_FILTERSIZE']))
        logger.debug("Simple filter size: %d" % (filtersize))
        if (filtersize < 3): filtersize = 3
        if (filtersize > 7): filtersize = 7

        # in a next step, modify the backsize parameter. Make sure it does not
        # become too large or too small
        backsize = (target_dimension * 60. / pixelscale) / filtersize

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
        logger.info("Adding background parameters:\n\n"+bg_opts+"\n\n")
        swarp_opts += bg_opts


    logger = logging.getLogger("SwarpStack - FinalStack")
    logger.info("Starting final stacking...")
    # print swarp_opts

    swarp_cmd = "%s %s %s" % (swarp_exec, swarp_opts, " ".join(single_prepared_files))
    logger.debug(" ".join(swarp_cmd.split()))
    try:
        ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        (swarp_stdout, swarp_stderr) = ret.communicate()

        logger.debug("swarp stdout:\n"+swarp_stdout)
        if (len(swarp_stderr) > 0 and ret.returncode != 0):
            logger.warning("swarp stderr:\n"+swarp_stderr)
        else:
            logger.debug("swarp stderr:\n"+swarp_stderr)
        #print "\n".join(swarp_stderr)
        logger.info("done, swarp returned (ret-code: %d)!" % ret.returncode)
    except OSError as e:
        podi_logging.log_exception()
        print >>sys.stderr, "Execution failed:", e

    return

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        import podi_swarpstack as me
        print me.__doc__
        sys.exit(0)

    # Setup everything we need for logging
    options = set_default_options()
    podi_logging.setup_logging(options)

    try:
        swarpstack()
    except KeyboardInterrupt, SystemExit:
        pass
    except:
        podi_logging.log_exception()
        pass
    finally:
        podi_logging.shutdown_logging(options)
