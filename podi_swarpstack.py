#!/usr/bin/env python


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
            log_exception()
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
                 -c $(swarp_default)s \
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
                log_exception()
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
                 %(prepared_files)s \
                 """ % dic

    logger = logging.getLogger("SwarpStack - FinalStack")
    logger.info("Starting final stacking...")
    # print swarp_opts

    swarp_cmd = "%s %s" % (swarp_exec, swarp_opts)
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
        log_exception()
        print >>sys.stderr, "Execution failed:", e

    return

if __name__ == "__main__":

    # Setup everything we need for logging
    options = set_default_options()
    setup_logging(options)

    try:
        swarpstack()
    except KeyboardInterrupt, SystemExit:
        pass
    except:
        log_exception()
        pass
    finally:
        shutdown_logging(options)
