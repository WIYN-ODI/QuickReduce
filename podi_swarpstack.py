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


import os
import sys
import pyfits
import subprocess
import math

from podi_commandline import *
from podi_definitions import *
#from podi_collectcells import *
import podi_sitesetup as sitesetup
import multiprocessing

import podi_logging
import logging
import socket
import tempfile
import shutil
import warnings

try:
    dir, _ = os.path.split(os.path.abspath(sys.argv[0]))
    sys.path.append(dir+"/test")
    import ephemerides
    import podi_ephemerides
except:
    pass

def mp_prepareinput(input_queue, output_queue, swarp_params, options):

    while (True):

        cmd = input_queue.get()
        if (cmd == None):
            input_queue.task_done()
            break

        input_file = cmd
        logger = logging.getLogger("MP-Prep( %s )" % (os.path.basename(input_file)))


        ret = {
            "master_reduction_files": {},
            "corrected_file": None,
            'exptime': 0,
            'mjd_obs_start': 0,
            'mjd_obs_end': 0,
        }
            

        try:
            hdulist = pyfits.open(input_file)
        except IOError:
            logger.error("Can't open file %s" % (inputlist[i]))
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
        master_reduction_files_used = collect_reduction_files_used(
            {}, {"calibrated": input_file} ) #ret['corrected_file']})


        # Assemble the temporary filename for the corrected frame
        suffix = None
        # Now construct the output filename
        if (corrected_filename and swarp_params['use_nonsidereal']):
            suffix = "nonsidereal"
        if (corrected_filename == None and swarp_params['use_ephemerides']):
            suffix = "ephemerides"
        if (corrected_filename == None and options['skip_otas'] != []):
            suffix = "otaselect"
        if (corrected_filename == None and not options['bpm_dir'] == None):
            suffix = "bpmfixed"

        if (not suffix == None):
            corrected_filename = "%(single_dir)s/%(obsid)s.%(suffix)s.fits" % {
                "single_dir": swarp_params['unique_singledir'],
                "obsid": hdulist[0].header['OBSID'],
                "suffix": suffix,
            }

        #
        # Check if we need to apply any corrections
        #
        if (corrected_filename == None or 
            (os.path.isfile(corrected_filename) and swarp_params['reuse_singles'])
        ):
            # Either we don't need to apply any corrections or we can re-use an 
            # older file with these corrections already applied

            ret['corrected_file'] = input_file
        else:

            if (swarp_params['use_nonsidereal']):
                logger.info("Applying the non-sidereal motion correction")

                # First apply the non-sidereal correction to all input frames
                # Keep frames that are already modified
                from podi_collectcells import apply_nonsidereal_correction
                # print options['nonsidereal']
                apply_nonsidereal_correction(hdulist, options, logger)

                try:
                    if (os.path.isfile(options['nonsidereal']['ref'])):
                        master_reduction_files_used = collect_reduction_files_used(
                            master_reduction_files_used, 
                            {"nonsidereal-reference": options['nonsidereal']['ref']})
                except:
                    pass

            if (swarp_params['use_ephemerides']):
                # get MJD of current frame
                mjd_thisframe = hdulist[0].header['MJD-OBS']
                mjd_ref = swarp_params['ephemerides']['ref-mjd']

                # print "\n"*10
                # print "ref-mjd",mjd_ref
                # print "this mjd",mjd_thisframe
                # print "min mjd", numpy.min(swarp_params['ephemerides']['data'][:,0])
                # print "max mjd", numpy.max(swarp_params['ephemerides']['data'][:,0])
                # print "\n"*10

                # now compute the Ra/Dec of the target in both the reference 
                # frame and in this frame
                ra_ref = swarp_params['ephemerides']['ra'](mjd_ref)
                ra_this = swarp_params['ephemerides']['ra'](mjd_thisframe)

                dec_ref = swarp_params['ephemerides']['dec'](mjd_ref)
                dec_this = swarp_params['ephemerides']['dec'](mjd_thisframe)

                # print "\n"*5, "ra//dec = ", ra_ref, dec_ref, "\n"*5
                # The Ra/Dec correction is how much the object has moved 
                # (as derived from the ephemerides) between the reference MJD 
                # and the timestamp of this frame
                d_ra = ra_ref - ra_this
                d_dec = dec_ref - dec_this
                # print d_ra, d_dec
                d_days = mjd_thisframe - mjd_ref
                logger.info("Applying ephemerid correction (%+.6f deg, %+.6f deg, dT=%.3f days)" % (
                    d_ra, d_dec, d_days))

                # Now apply these corrections to all extensions with an
                # apparently valid WCS system
                orig_ra, orig_dec = None, None
                for ext in hdulist:
                    if ('CRVAL1' in ext.header and
                        'CRVAL2' in ext.header):
                        # print ext.header['CRVAL1'], ext.header['CRVAL2'], (ext.header['EXTNAME'] if 'EXTNAME' in ext.header else "??")
                        orig_ra = ext.header['CRVAL1'] if orig_ra == None else orig_ra
                        orig_dec = ext.header['CRVAL2'] if orig_dec == None else orig_dec
                        ext.header['CRVAL1'] += d_ra
                        ext.header['CRVAL2'] += d_dec
                        # print ext.header['CRVAL1'], ext.header['CRVAL2'], (ext.header['EXTNAME'] if 'EXTNAME' in ext.header else "??")

                logger.debug("Pre-correction Ra/Dec was: %12.7f  %+12.7f" % (orig_ra, orig_dec))
                logger.debug("Post-corrected Ra/Dec is: %12.7f %+12.7f" % (orig_ra + d_ra, orig_dec + d_dec))



            if (options['skip_otas'] != []):
                logger.info("Skipping some OTAs")
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

            if (not options['bpm_dir'] == None):
                logger.info("Applying bad-pixel masks")
                for ext in range(len(hdulist)):
                    if (not is_image_extension(hdulist[ext])):
                        continue

                    fppos = None
                    if ('FPPOS' in hdulist[ext].header):
                        fppos = hdulist[ext].header['FPPOS']
                    if (not fppos == None):
                        region_file = "%s/bpm_%s.reg" % (options['bpm_dir'], fppos)
                        if (os.path.isfile(region_file)):
                            mask_broken_regions(hdulist[ext].data, region_file)
                            master_reduction_files_used = collect_reduction_files_used(
                                master_reduction_files_used, {"bpm": region_file})


            # Loop over all extensions and only select those that are not marked as guide chips
            if (True): #options['skip_otas'] != []):
                logger.info("Sorting out guide-OTAs")
                ota_list = []
                for ext in hdulist:
                    if ('CELLMODE' in ext.header and
                        ext.header['CELLMODE'].find("V") >= 0):
                        logger.debug("skipping ota %s as requested" % (ext.header['EXTNAME']))
                        continue
                    ota_list.append(ext)

                # Save the modified OTA list for later
                hdulist = pyfits.HDUList(ota_list)

            # Check if the corrected file already exists - if not create it
            #if (not os.path.isfile(corrected_filename)):
            logger.info("Writing correctly prepared file--> %s" % (corrected_filename))

            clobberfile(corrected_filename)
            hdulist.writeto(corrected_filename, clobber=True)

            # Now change the filename of the input list to reflect 
            # the corrected file
            ret['corrected_file'] = corrected_filename
    
        #
        # Now we have the filename of the file to be used for the swarp-input
        #
        ret["master_reduction_files"] = master_reduction_files_used
        logger.debug("Sending return value to master process")
        output_queue.put(ret)
        input_queue.task_done()

        # print input_file, "\n", ret["master_reduction_files"]

    # end of routine



def prepare_input(inputlist, swarp_params, options):

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
    for i in range(len(inputlist)):
        if (not os.path.isfile(inputlist[i])):
            continue
        try:
            hdulist = pyfits.open(inputlist[i])
        except IOError:
            logger.error("Can't open file %s" % (inputlist[i]))
            inputlist[i] = None
            continue

        in_queue.put(inputlist[i])
        n_jobs += 1

    logger.info("Queued %d jobs ..." % (n_jobs))

    #
    # Start worker processes
    #
    worker_args = (in_queue, out_queue, swarp_params, options)
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

    stack_start_time = 1e9
    stack_end_time = -1e9
    stack_total_exptime = 0
    stack_framecount = 0


    for i in range(n_jobs):
        ret = out_queue.get()
        logger.debug("Received results from job %d" % (i+1))


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

        master_reduction_files_used = collect_reduction_files_used(master_reduction_files_used, 
                                                                   ret['master_reduction_files'])

        corrected_file_list.append(ret['corrected_file'])

    #
    # By now all frames have all corrections applied,
    # so we can go ahead and stack them as usual
    #
    
    # Make sure to join/terminate all processes
    for p in processes:
        p.join()
        
    logger.info("All files prepared!")

    return corrected_file_list, stack_total_exptime, stack_framecount, stack_start_time, stack_end_time, master_reduction_files_used 


def mp_swarp_single(sgl_queue, dum):

    while(True):
        cmd = sgl_queue.get()
        if (cmd == None):
            sgl_queue.task_done()
            break

        swarp_cmd, prepared_file, single_file = cmd
        logger = logging.getLogger("MPSwarpSgl(%s)" % (os.path.basename(single_file)))

        hdulist = pyfits.open(prepared_file)
        obsid = hdulist[0].header['OBSID']

        try:
            logger.debug(" ".join(swarp_cmd.split()))
            ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
            (swarp_stdout, swarp_stderr) = ret.communicate()
            logger.debug("swarp stdout:\n"+swarp_stdout)
            if (len(swarp_stderr) > 0 and ret.returncode != 0):
                logger.warning("swarp stderr:\n"+swarp_stderr)
            else:
                logger.debug("swarp stderr:\n"+swarp_stderr)

            # Add some basic headers from input file to the single file
            # this is important for the differencing etc.
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
            hdu_single.flush()

            #print "\n".join(swarp_stderr)
            # single_prepared_files.append(single_file)
        except OSError as e:
            podi_logging.log_exception()
            print >>sys.stderr, "Execution failed:", e

        sgl_queue.task_done()



def swarpstack(outputfile, inputlist, swarp_params, options, keep_intermediates=False, unique_dir=None):

    logger = logging.getLogger("SwarpStack - Prepare")

    # Figure out the config path
    swarp_default = "%s/.config/swarp.default" % (sitesetup.exec_dir)
    logger.debug("Using swarp-default in %s" % (swarp_default))

    if (len(inputlist) <= 0):
        logger.error("No input files specified!")
        return

    master_reduction_files_used = {}
    
    logger.info("Removing some OTAs from input: %s" % (str(options['skip_otas'])))

    #
    # Construct a unique name to hold intermediate files
    #
    if (not unique_dir == None and os.path.isdir(unique_dir)):
        unique_singledir = unique_dir
        keep_intermediates = True
    else:
        process_id = os.getpid()
        hostname = socket.gethostname()
        # Create a temporary directory
        unique_singledir = tempfile.mkdtemp(dir=sitesetup.swarp_singledir,
                                            prefix="%s-%05d----" % (hostname, process_id))


    ############################################################################
    #
    # Prepare all QR'ed input files, applying additional corrections where needed
    #
    ############################################################################

    logger.info("Storing intermediate files in %s ..." % (unique_singledir))
    swarp_params['unique_singledir'] = unique_singledir

    stack_start_time = 1e9
    stack_end_time = -1e9
    stack_total_exptime = 0
    stack_framecount = 0

    # print "input=",inputlist
    # print "output=",outputfile

    modified_files, stack_total_exptime, stack_framecount, \
        stack_start_time, stack_end_time, master_reduction_files_used = \
                                prepare_input(inputlist, swarp_params, options)
    # print modified_files
    inputlist = modified_files


    add_only = swarp_params['add'] and os.path.isfile(outputfile)
    if (add_only):
        logger.info("Activating ADD mode")

    if (outputfile.endswith(".fits")):
        outputfile = outputfile[:-5]

    header_only_file = "%s/preswarp.fits" % (unique_singledir)
    logger.debug("Using header-only-file: %s" % (header_only_file))

    # Make sure the reference file is a valid file
    if (not swarp_params['reference_file'] == None and os.path.isfile(swarp_params['reference_file'])):
        logger.info("Using %s as reference file" % (swarp_params['reference_file']))
    else:
        swarp_params['reference_file'] = None

    logging.debug("Using modified input list: %s" % (str(inputlist)))

    ############################################################################
    #
    # Figure out the pixel-grid and sky-coverage of the final stack 
    #
    ############################################################################
    if (add_only or not swarp_params['reference_file'] == None):
        #
        # This is the simpler add-only mode
        #

        if (not swarp_params['reference_file'] == None):
            output_info = pyfits.open(swarp_params['reference_file'])
        else:
            # Open the existing output header and get data from there
            output_info = pyfits.open(outputfile+".fits")

        logger.info("Stack information...")
        logger.info("   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header))

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']

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
        if (swarp_params['no-fluxscale']):
            swarp_opts += " -FSCALE_KEYWORD none "

        swarp_opts += " -SUBTRACT_BACK %s " % ("Y" if swarp_params['subtract_back'] else "N")

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
        
        if (swarp_params['pixelscale'] <= 0):
            swarp_params['pixelscale'] = math.fabs(output_info[0].header['CD1_1']) * 3600.
            #pixelscale = (output_info[0].header['CD1_1'] * output_info[0].header['CD2_2'] \
            #             - output_info[0].header['CD1_2'] * output_info[0].header['CD2_1']) * 3600.
            logger.info("Computing pixelscale from data: %.4f arcsec/pixel" % (swarp_params['pixelscale']))
    
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

    for prepared_file in inputlist:
        hdulist = pyfits.open(prepared_file)
        obsid = hdulist[0].header['OBSID']

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
               'fluxscale': 'none' if swarp_params['no-fluxscale'] else 'FLXSCALE'
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
                 -FSCALE_KEYWORD %(fluxscale)s \
                 %(inputfile)s \
                 """ % dic

        single_file = "%(singledir)s/%(obsid)s.fits" % dic

        # print swarp_opts
        swarp_cmd = "%s %s" % (sitesetup.swarp_exec, swarp_opts)

        if (add_only and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, skipping it" % (single_file))
        elif (swarp_params['reuse_singles'] and os.path.isfile(single_file)):
            logger.info("This single-swarped file (%s) exist, re-using it" % (single_file))
            single_prepared_files.append(single_file)
        else:
            logger.info("Preparing file %s, please wait ..." % (prepared_file))
            logger.debug(" ".join(swarp_cmd.split()))

            sgl_queue.put( (swarp_cmd, prepared_file, single_file) )
            single_prepared_files.append(single_file)

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
    if (swarp_params['target_dimension'] > 0 and swarp_params['subtract_back']):
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


    if (swarp_params['subtract_back']):
        logger = logging.getLogger("SwarpStack - SkySub")
        logger.info("Performing sky-subtraction on all frames")

        final_prepared_files = []

        # Prepare the worker queue
        sgl_queue = multiprocessing.JoinableQueue()

        for prepared_file in single_prepared_files:
            hdulist = pyfits.open(prepared_file)
            obsid = hdulist[0].header['OBSID']

            
            # assemble all swarp options for that run
            bgsub_file = "%(singledir)s/%(obsid)s.bgsub.fits" % {
                'singledir': unique_singledir,
                'obsid': obsid,}
            bgsub_weight_file = "%(singledir)s/%(obsid)s.bgsub.weight.fits" % {
                'singledir': unique_singledir,
                'obsid': obsid,}

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
                   'bgsub': "Y" if swarp_params['subtract_back'] else "N",
                   'bgsub_file': bgsub_file,
                   'bgsub_weight_file': bgsub_weight_file,
                   'bgopts': bg_opts,
                   'inputfile': prepared_file,
               }
            
            swarp_opts = """\
                     -c %(swarp_default)s \
                     -IMAGEOUT_NAME %(singledir)s/%(obsid)s.bgsub.fits \
                     -WEIGHTOUT_NAME %(singledir)s/%(obsid)s.bgsub.weight.fits \
                     -PIXEL_SCALE %(pixelscale)f \
                     -PIXELSCALE_TYPE %(pixelscale_type)s \
                     -COMBINE Y \
                     -COMBINE_TYPE AVERAGE \
                     -CENTER_TYPE MANUAL \
                     -CENTER %(center_ra)f,%(center_dec)f \
                     -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
                     -RESAMPLE_DIR %(resample_dir)s \
                     -RESAMPLE Y \
                     -SUBTRACT_BACK %(bgsub)s \
                     -WEIGHT_TYPE MAP_WEIGHT \
                     -DELETE_TMPFILES Y \
                     -WRITE_FILEINFO Y
                     %(bgopts)s \
                     %(inputfile)s \
                     """ % dic

            
            # print swarp_opts
            swarp_cmd = "%s %s" % (sitesetup.swarp_exec, swarp_opts)
            # print swarp_cmd

            if (add_only and os.path.isfile(bgsub_file)):
                logger.info("This single-swarped file (%s) exist, skipping it" % (bgsub_file))
            elif (swarp_params['reuse_singles'] and os.path.isfile(bgsub_file)):
                logger.info("This single-swarped file (%s) exist, re-using it" % (bgsub_file))
                final_prepared_files.append(bgsub_file)
            else:
                logger.info("Preparing file %s, please wait ..." % (prepared_file))
                logger.debug(" ".join(swarp_cmd.split()))

                sgl_queue.put( (swarp_cmd, prepared_file, bgsub_file) )
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

    #############################################################################
    #
    # Now all single files are prepared, go ahead and produce the actual stack
    # Use the background-subtracted or single files from above as direct input, 
    # i.e. do not resample these files as they are already on the right pixel grid.
    #
    #############################################################################
    dic['combine_type'] = swarp_params['combine-type'] #"AVERAGE"
    dic['imageout'] = outputfile+".fits"
    dic['weightout'] = outputfile+".weight.fits"
    dic['prepared_files'] = " ".join(single_prepared_files)
    dic['bgsub'] = "N" # as this was done before if swarp_params['subtract_back'] else "N"
    dic['clip-sigma'] = swarp_params['clip-sigma']
    dic['clip-ampfrac'] = swarp_params['clip-ampfrac']

    swarp_opts = """\
                 -c %(swarp_default)s \
                 -IMAGEOUT_NAME %(imageout)s \
                 -WEIGHTOUT_NAME %(weightout)s \
                 -COMBINE_TYPE %(combine_type)s \
                 -PIXEL_SCALE %(pixelscale)f \
                 -PIXELSCALE_TYPE %(pixelscale_type)s \
                 -COMBINE Y \
                 -COMBINE_TYPE %(combine_type)s \
                 -CLIP_AMPFRAC %(clip-ampfrac)f \
                 -CLIP_SIGMA %(clip-sigma)f \
                 -CENTER_TYPE MANUAL \
                 -CENTER %(center_ra)f,%(center_dec)f \
                 -IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
                 -RESAMPLE N \
                 -RESAMPLE_DIR %(singledir)s \
                 -SUBTRACT_BACK %(bgsub)s \
                 -WEIGHT_TYPE MAP_WEIGHT \
                 -DELETE_TMPFILES N \
                 -WRITE_FILEINFO Y
                 """ % dic


    logger = logging.getLogger("SwarpStack - FinalStack")
    logger.info("Starting final stacking...")
    # print swarp_opts

    swarp_cmd = "%s %s %s" % (sitesetup.swarp_exec, swarp_opts, " ".join(final_prepared_files))
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

    logger.info("Stack (%s) complete, adding headers" % (dic['imageout']))

    # Finally, open the output file and copy a bunch of headers into it
    hdustack = pyfits.open(dic['imageout'], mode='update')
    # Also open the first frame in the stack to be used as data source
    firsthdu = pyfits.open(inputlist[0])

    for hdrkey in [
            'TARGRA', 'TARGDEC',
            'FILTER', 'FILTERID', 'FILTDSCR', 
            'OBSID', 'OBJECT', 
            'DATE-OBS', 'TIME-OBS', 'MJD-OBS']:
        if (hdrkey in firsthdu[0].header):
            key, val, com = firsthdu[0].header.cards[hdrkey]
            hdustack[0].header[key] = (val, com)

    hdustack[0].header['EXPTIME'] = (stack_total_exptime, "total exposure time in stack")
    hdustack[0].header['MJD-STRT'] = (stack_start_time, "MJD at start of earliest exposure")
    hdustack[0].header['MJD-END'] = (stack_end_time, "MJD at end of last exposure")
    hdustack[0].header['NCOMBINE'] = (stack_framecount, "number of exposures in stack")

    # Add some additional headers
    hdustack[0].header['MAGZERO'] = 25.
    hdustack[0].header['_BGSUB'] = "yes" if swarp_params['subtract_back'] else "no"
    hdustack[0].header['_PXLSCLE'] = swarp_params['pixelscale']
    hdustack[0].header['_RESGL'] = ("yes" if swarp_params['reuse_singles'] else "no", "reuse singles?")
    hdustack[0].header['_NONSDRL'] = ("yes" if swarp_params['use_nonsidereal'] else "no", "using non-sidereal correction?")
    valid_nonsidereal_reference = False
    try:
        hdustack[0].header['_NSID_RA'] = options['nonsidereal']['dra']
        hdustack[0].header['_NSID_DE'] = options['nonsidereal']['ddec']
        hdustack[0].header['_NSID_RT'] = options['nonsidereal']['ref_mjd']
        hdustack[0].header['_NSID_RF'] = options['nonsidereal']['ref']
        if (os.path.isfile(options['nonsidereal']['ref'])):
            valid_nonsidereal_reference = True
    except:
        pass
    firsthdu.close()

    if (valid_nonsidereal_reference):
        firsthdu = pyfits.open(options['nonsidereal']['ref'])
        try:
            hdustack[0].header['NSR-OBST'] = firsthdu[0].header['TIME-OBS']
            hdustack[0].header['NSR-OBSD'] = firsthdu[0].header['DATE-OBS']
            hdustack[0].header['NSR-OBSM'] = firsthdu[0].header['MJD-OBS']

            hdustack[0].header['NSR-MIDT'] = firsthdu[0].header['TIME-MID']
            hdustack[0].header['NSR-MIDD'] = firsthdu[0].header['DATE-MID']
            hdustack[0].header['NSR-MIDM'] = firsthdu[0].header['MJD-MID']

            hdustack[0].header['NSR-ENDT'] = firsthdu[0].header['TIME-END']
            hdustack[0].header['NSR-ENDD'] = firsthdu[0].header['DATE-END']
            hdustack[0].header['NSR-ENDM'] = firsthdu[0].header['MJD-END']
        except:
            pass

    # Add the user-defined keywords to the stacked file. This is required for
    # proper integration with the PPA framework.
    # print options['additional_fits_headers']
    for key, value in options['additional_fits_headers'].iteritems():
        hdustack[0].header[key] = (value, "user-added keyword")

    #
    # Create an association table from the master reduction files used.
    # 
    # print master_reduction_files_used
    assoc_table = create_association_table(master_reduction_files_used)
    hdustack.append(assoc_table)

    hdustack.flush()
    hdustack.close()
    
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




def read_swarp_params(filelist):

    params = {}

    logger = logging.getLogger("SwarpStack - Config")

    params['pixelscale'] = float(cmdline_arg_set_or_default("-pixelscale", 0))
    params['subtract_back'] = cmdline_arg_isset("-bgsub")
    params['reuse_singles'] = cmdline_arg_isset("-reusesingles")
    params['use_nonsidereal'] = cmdline_arg_isset("-nonsidereal")
    params['target_dimension'] = float(cmdline_arg_set_or_default('-dimension', -1))
    params['add'] = cmdline_arg_isset("-add")
    params['reference_file'] = cmdline_arg_set_or_default("-reference", None)
    params['no-fluxscale'] = cmdline_arg_isset('-nofluxscale')

    combine_method = cmdline_arg_set_or_default('-combine', 'average')
    if (not combine_method.lower() in ['average', 'median', 'sum', 'min', 'max', 'weighted', 'chi2',
                                       'chi-old', 'chi-mode', 'chi-mean', 'clipped',
                                       'weighted_weight', 'median_weight', 'and', 'nand', 'or', 'nor']):
        logger = logging.getLogger("Setup")
        logger.error("The specified combine method (%s) is not supported, using average instead" % (combine_method))
        combine_method = 'average'
    params['combine-type'] = combine_method.upper()

    params['use_ephemerides'] = cmdline_arg_isset("-ephemerides")
    if (params['use_ephemerides']):
        opts = cmdline_arg_set_or_default('-ephemerides', None)
        # print opts
        if (opts == None):
            params['use_ephemerides'] = False
        else:
            items = opts.split(',')
            # print items

            # See if the first parameter is a number, if so it's the reference MJD,
            # if not, assume it's a FITS file and we are to read MJD-OBS from the header
            ref_mjd = None
            try:
                ref_mjd = float(items[0])
            except:
                hdulist = pyfits.open(items[0])
                for ext in hdulist:
                    if ('MJD-OBS' in ext.header):
                        ref_mjd = ext.header['MJD-OBS']
                        break
            if (ref_mjd == None):
                params['use_ephemerides'] = False
            else:

                if (items[1] == "NASA"):
                    object_name = items[2]
                    results = podi_ephemerides.get_ephemerides_for_object_from_filelist(
                        object_name=object_name, 
                        filelist=filelist,
                        session_log_file="swarpstack_horizon.session",
                        verbose=False
                    )

                    params['ephemerides'] = {
                        'ref': items[0],
                        'ref-mjd': ref_mjd,
                        'datafile': "NASA-TELNET",
                        'data': results['data'],
                        'ra': results['ra'],
                        'dec': results['dec'],
                    }
                else:
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
                    }

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
            
    return params

if __name__ == "__main__":
    if (len(sys.argv) <= 1):
        import podi_swarpstack as me
        print me.__doc__
        sys.exit(0)

    # Setup everything we need for logging
    options = set_default_options()
    podi_logging.setup_logging(options)

    logger = logging.getLogger("SwarpStack-Startup")
    if (cmdline_arg_isset("-fromfile")):
        configfile = get_cmdline_arg("-fromfile")
        if (os.path.isfile(configfile)):
            logger.info("Reading additional command line parameters from file (%s)" % (configfile))
            conf = open(configfile, "r")
            lines = conf.readlines()
            for line in lines:
                line = line.strip()
                if (len(line) <= 0):
                    continue
                # print "Adding ",line.strip(),"to command line"
                if (not line.startswith("#")):
                    sys.argv.append(line)
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
        inputlist = get_clean_cmdline()[2:]
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

    except KeyboardInterrupt, SystemExit:
        pass
    except:
        podi_logging.log_exception()
        pass
    finally:
        podi_logging.shutdown_logging(options)
