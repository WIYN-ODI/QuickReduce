#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
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


import Queue
import threading
import multiprocessing
import ctypes
import matplotlib.pyplot
import time

gain_correct_frames = False
from podi_definitions import *
import podi_findstars
import podi_search_ipprefcat
import podi_fixwcs
import podi_sitesetup as sitesetup
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

def read_reduction_directories(start=1, warn=True, verbose=True):
    #
    # Read other parameters, specifying the directories for the 
    # flatfields, darks and biases
    #
    # Set all reduction folder to None to mask them as not set
    flatfield_dir = None
    bias_dir = None
    dark_dir = None
    bpm_dir = None

    if (cmdline_arg_isset("-cals")):
        bias_dir = get_cmdline_arg("-cals")
        dark_dir = get_cmdline_arg("-cals")
        flatfield_dir = get_cmdline_arg("-cals")

    bias_dir = cmdline_arg_set_or_default("-bias", bias_dir)
    dark_dir = cmdline_arg_set_or_default("-dark", dark_dir)
    flatfield_dir = cmdline_arg_set_or_default("-flat", flatfield_dir)

    bpm_dir = cmdline_arg_set_or_default("-bpm", bpm_dir)
    if (bpm_dir == "auto"):
        full_path = os.path.abspath(sys.argv[0])
        bpm_dir, dummy = os.path.split()

    # Output some summary on the reduction
    if (verbose):
        print """
Calibration data:
            Bias: %s
            Dark: %s
      Flatfields: %s
  Bad pixel mask: %s
""" % (bias_dir, dark_dir, flatfield_dir, bpm_dir)

    i = 0
    return bias_dir, dark_dir, flatfield_dir, bpm_dir, i

def collect_reduce_ota(filename,
                       bias_dir, dark_dir, flatfield_dir, bpm_dir,
                       offset_pointing=[0,0], offset_dither=[0,0], target_coords=None,
                       pixelvalue_indef=numpy.NaN,
                       wcs_solution=None,
                       prepare_fixwcs=False,
                       hardcoded_detsec=False,
                       verbose=False):

    if (not os.path.isfile(filename)):
        stdout_write("Couldn't find file %s ..." % (filename))
    else:
        # Create an fits extension to hold the output
        hdu = pyfits.ImageHDU()

        try:
            hdulist = pyfits.open(filename, memmap=False)
        except:
            # This happed with corrupt files
            return None, None

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

        #
        # Allocate memory for the merged frame, and set all pixels by default to NaN.
        # Valid pixels will subsequently be overwritten with real numbers
        #
        merged = numpy.ones(shape=(size_x, size_y), dtype=numpy.float32)
        merged[:,:] = pixelvalue_indef
        
        for cell in range(1,65):
            stdout_write("\r%s:   OTA %02d, cell %s ..." % (obsid, ota, hdulist[cell].header['EXTNAME']))

            # Check if this is one of the broken cells
            wm_cellx, wm_celly = hdulist[cell].header['WN_CELLX'], hdulist[cell].header['WN_CELLY']
            broken = False
            list_of_broken_cells = broken_ota_cells[ota_name]
            for broken_cell in list_of_broken_cells:
                x,y = broken_cell
                #print x,y
                if (wm_cellx == x and wm_celly == y):
                    broken = True
                    #print "found broken cell", hdulist[cell].header['EXTNAME'],broken_cell
                    break

            # If not, overscan subtract and insert into large frame
            if (not broken):
                overscan_region = extract_region(hdulist[cell].data, '[500:530,1:494]')
                overscan_level = numpy.median(overscan_region)

                hdulist[cell].data -= overscan_level

                if (gain_correct_frames):
                # Correct for the gain variations in each cell
                    try:
                        gain = float(hdulist[cell].header['GAIN'])
                        hdulist[cell].data *= gain
                    except:
                        print "Couldn't find the GAIN header!"
                        pass

                if (True): #hardcoded_detsec):
                    cell_xy = hdulist[cell].header['EXTNAME'][2:4]
                    x, y = int(cell_xy[0]), int(cell_xy[1])

                    datasec = '[1:480,1:494]'
                    _y = 7 - y
                    y1 = 1 + (505*_y)  #was 503 
                    #taken from ODI OTA Technical Datasheet (det area 480x494, streets 11/28 px)
                    y2 = y1 + 493
                    x1 = 1 + 508 * x
                    x2 = x1 + 479
                    detsec = '[%d:%d,%d:%d]' % (x1, x2, y1, y2)

                    insert_into_array(hdulist[cell].data, 
                                      datasec, merged, detsec)
                else:
                    insert_into_array(hdulist[cell].data, 
                                      hdulist[cell].header['DATASEC'],
                                      merged,
                                      hdulist[cell].header['DETSEC'])

            #
            # Special case for cell 0,7 (the one in the bottom left corner):
            # Copy the CRPIX values into the merged image header 
            #
            if (hdulist[cell].header['EXTNAME'] == "XY07"):
                # print "Setting CRPIXs", hdulist[cell].header['CRPIX1'], hdulist[cell].header['CRPIX2']
                hdu.header.update("CRPIX1", hdulist[cell].header['CRPIX1'], "Ref. pixel RA")
                hdu.header.update("CRPIX2", hdulist[cell].header['CRPIX2'], "Ref. pixel DEC")
                
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
        if (not bias_dir == None):
            bias_filename = "%s/bias.fits" % (bias_dir)
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
        if (not dark_dir == None):

            # For now assume all detectors are switched on
            detectorglow = "yes"

            dark_filename = "%s/dark_%s.fits" % (dark_dir, detectorglow)
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
        if (not flatfield_dir == None):
            flatfield_filename = "%s/flat_%s.fits" % (flatfield_dir, filter_name)
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
        # Determine which region file we need
        if (not bpm_dir == None):
            region_file = "%s/bpm_%s.reg" % (bpm_dir, fppos)
            if (os.path.isfile(region_file)):
                # Apply the bad pixel regions to file, marking
                # all bad pixels as NaNs
                mask_broken_regions(merged, region_file)
                hdu.header.add_history("CC-BPM: %s" % (os.path.abspath(region_file)))

        # Insert the DETSEC header so IRAF understands where to put the extensions
        start_x = ota_c_x * 4096
        start_y = ota_c_y * 4096        
        end_x = start_x + det_x2 - det_x1
        end_y = start_y + det_y2 - det_y1
        detsec_str = "[%d:%d,%d:%d]" % (start_x, end_x, start_y, end_y)
        hdu.header.update("DETSEC", detsec_str, "position of OTA in focal plane")
                
        if (cmdline_arg_isset("-simplewcs") or wcs_solution != None):
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
        if (not target_coords == None):
            ra, dec = target_coords
            coord_j2000 = ephem.Equatorial(ra, dec, epoch=ephem.J2000)

        # Write the CRVALs with the pointing information
        #print numpy.degrees(coord_j2000.ra), numpy.degrees(coord_j2000.dec)  
        hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  
        hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) 

        # Compute total offsets as the sum from pointing and dither offset
        offset_total = numpy.array(offset_pointing) + numpy.array(offset_dither)

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

        # Now add the canned WCS solution
        if (wcs_solution != None):
            #print "Adding header from WCS minifits (%s)" % (extname)
            wcs = pyfits.open(wcs_solution)
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

        if (prepare_fixwcs):
            # Create source catalog
            
            if (False):
                source_cat = podi_findstars.find_stars(hdu, binning=4, boxsize=24, dumpfile=None, verbose=False,
                                                   detect_threshold=1., detect_minarea=4, roundness_limit=[-0.2,+0.2],
                                                   max_starcount=150,
                                                   extension_id=ota)
            else:
                hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdu.header, data=hdu.data)])
                fitsfile = "/scratch/tmp_OTA%02d.fits" % (ota)
                catfile = "/scratch/tmp_OTA%02d.cat" % (ota)
                hdulist.writeto(fitsfile, clobber=True)
                sexcmd = "sex -c /work/podi_devel/.config/wcsfix.sex -CATALOG_NAME %s %s" % (catfile, fitsfile)
                os.system(sexcmd)
                try:
                    source_cat = numpy.loadtxt(catfile)
                    if (source_cat.shape[0] == 0):
                        source_cat == None
                    else:
                        print "source-cat:", source_cat.shape
                        source_cat[:,12] = ota
                        #valid = source_cat[:,11] == 0
                        #source_cat = source_cat[valid]
                except:
                    source_cat == None
                clobberfile(fitsfile)

            fixwcs_data = None
            if (source_cat != None):
                if (source_cat.shape[0] > 0):
                    odi_ra = source_cat[:,0]
                    odi_dec = source_cat[:,1]
                    odi_mag = -2.5 * numpy.log10(source_cat[:,6]) + 30

                    # Read the reference catalog
                    center_ra, center_dec = center_coords(hdu.header)
                    search_size = (8+3) * (1./60.)
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
                    dx, dy, n = podi_fixwcs.shift_align_wcs(ota_odi_ra, ota_odi_dec, ota_ref_ra, ota_ref_dec)
                    print "WCSFIX dx/dy =", dx, dy
                    fixwcs_data = (odi_ra, odi_dec, ref_ra, ref_dec, dx, dy, source_cat)
        else:
            fixwcs_data = None
            
    return hdu, fixwcs_data
    


#########
#
# This routine is a wrapper around the actual collect_reduce_ota routine,
# mainly dealing with the parallelization and inter-process communication
#
#########
def parallel_collect_reduce_ota(queue, return_queue,
                                bias_dir, dark_dir, flatfield_dir, bpm_dir,
                                offset_pointing=[0,0], offset_dither=[0,0], target_coords=None,
                                pixelvalue_indef=numpy.NaN,
                                wcs_solution=None, prepare_fixwcs=False,
                                hardcoded_detsec=False):

    while (True):
        cmd_quit, filename, ota_id = queue.get()
        if (cmd_quit):
            queue.task_done()
            return

        # Do the work
        hdu, wcsfix_data = collect_reduce_ota(filename, 
                           bias_dir, dark_dir, flatfield_dir, bpm_dir,
                           offset_pointing=offset_pointing,
                           offset_dither=offset_dither,
                           target_coords=target_coords,
                           pixelvalue_indef=pixelvalue_indef,
                           wcs_solution=wcs_solution,
                           prepare_fixwcs=prepare_fixwcs,
                           hardcoded_detsec=hardcoded_detsec,
            )

        # Add the results to the return_queue so the master process can assemble the result file
        # print "Adding results for OTA",ota_id,"to return queue"
        return_queue.put( (hdu, ota_id, wcsfix_data) )
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
                 bias_dir=None, dark_dir=None, flatfield_dir=None, bpm_dir=None,
                 wcs_solution=None,
                 batchmode=False,
                 fixwcs=False,
                 hardcoded_detsec=False,
                 clobber_mode=True):

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

    if (outputfile.find("%") >= 0):
        # The output filename contains special tags that should 
        # be replaced by values from the file header

        filename = "%s/%s.%02d.fits" % (directory, filebase, 33)
        if (not os.path.isfile(filename)):
            filename = "%s/%s.%02d.fits.fz" % (directory, filebase, 33)
            if (not os.path.isfile(filename)):
                print "Error opening OTA 3,3 (%s)" % (filename)
                return -1
        try:
            hdulist = pyfits.open(filename)
        except:
            stdout_write("Problem opening a fits-file (%s), aborting!\n" % filename)
            return None

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
            elif (outputfile[start:start+8] == "%EXPTIME"):
                outputfile = "%s%.1f%s" % (outputfile[:start], header['EXPTIME'], outputfile[start+8:])
            else:
                stdout_write("found unknown tag in %s\n" % outputfile)
                break

        hdulist.close()
        del hdulist
        del header

        stdout_write("Replaced some keywords, new output filename: ---> %s\n" % (outputfile))

    if (os.path.isfile(outputfile) and not clobber_mode):
        print "#####################################################"
        print "#"
        print "# File %s already exists, skipping!" % (outputfile)
        print "#"
        print "#####################################################"
        print "\n\n\n\n"

    #
    # Read all offsets from command line
    # For convenience, there are two sets of offset parameters, that internally simply 
    # get added up. The reason for this is to make specifying them on the command line 
    # easier, since the pointing offsets stay constant across a dither pattern, while 
    # the dither offsets change.
    #
    _offset_pointing = cmdline_arg_set_or_default("-pointing", "0,0")
    dx,dummy,dy = _offset_pointing.partition(",")
    offset_pointing = [float(dx), float(dy)]

    _offset_dither = cmdline_arg_set_or_default("-dither", "0,0")
    dx,dummy,dy = _offset_dither.partition(",")
    offset_dither = [float(dx), float(dy)]

    target_coords = None
    if (cmdline_arg_isset("-target")):
        _target_coords = cmdline_arg_set_or_default("-target", "0,0")
        ra,dummy,dec = _target_coords.partition(",")
        target_coords = (ra, dec)


    # Start new list of HDUs
    ota_list = [None] * (len(available_ota_coords)+1)

    # And add the primary HDU to make the fits file a valid one
    primhdu = pyfits.PrimaryHDU()
    ota_list[0] = primhdu
    
    # Set the fallback value for undefined pixels (mostly the gaps between the OTA cells)
    # This works perfectly fine in ds9, but not SExtractor
    pixelvalue_indef = numpy.NaN
    if (cmdline_arg_isset("-prep4sex")):
        # Check if the user requested us to prepare the frame for SExtractor
        # SExtractor doesn't like NaNs, so replace all of them with something
        # more negative than -1e30 (that's -1 times SEx's BIG variable)
        pixelvalue_indef = -1e31
    

    #
    # Set up the parallel processing environment
    #
    #result_buffer = numpy.zeros(shape=(buffer.shape[0], buffer.shape[1]), dtype=numpy.float32)
    queue = multiprocessing.JoinableQueue()
    return_queue = multiprocessing.Queue()
    
    processes = []
    #for i in range(number_cpus):

    worker_args = (queue, return_queue,
                   bias_dir, dark_dir, flatfield_dir, bpm_dir,
                   offset_pointing,
                   offset_dither,
                   target_coords,
                   pixelvalue_indef,
                   wcs_solution,
                   fixwcs,
                   hardcoded_detsec,
                   )

    number_extensions = 0
    for ota_id in range(len(available_ota_coords)):
        ota_c_x, ota_c_y = available_ota_coords[ota_id]        
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
        queue.put( (False, filename, ota_id+1) )

    # Tell all workers to shut down when no more data is left to work on
    for i in range(len(processes)):
        stdout_write("Sending quit command!\n")
        queue.put((True,None,None))

    # Create all processes to handle the actual reduction and combination
    #print "Creating",number_cpus,"worker processes"
    for i in range(number_cpus):
        p = multiprocessing.Process(target=parallel_collect_reduce_ota, args=worker_args)
        p.start()
        processes.append(p)
        time.sleep(0.1)

    #
    # By now all workers have computed their HDUs or are busy doing so,
    # let's extract their results from the return queue and assemble the ota list
    #
    fixwcs_odi_ra, fixwcs_odi_dec = numpy.array([]), numpy.array([])
    fixwcs_ref_ra, fixwcs_ref_dec = numpy.array([]), numpy.array([])
    fixwcs_odi_sourcecat = None
    fixwcs_bestguess = numpy.ones(shape=(len(available_ota_coords),2)) * -1000

    for i in range(len(available_ota_coords)):
        hdu, ota_id, wcsfix_data = return_queue.get()
        if (hdu == None):
            continue

        ota_list[ota_id] = hdu

        if (fixwcs):
            
            if (wcsfix_data != None):
                odi_ra, odi_dec, ref_ra, ref_dec, dx, dy, source_cat = wcsfix_data

                fixwcs_odi_ra  = numpy.append(fixwcs_odi_ra,  odi_ra,  axis=0)
                fixwcs_odi_dec = numpy.append(fixwcs_odi_dec, odi_dec, axis=0)
                fixwcs_ref_ra  = numpy.append(fixwcs_ref_ra,  ref_ra,  axis=0)
                fixwcs_ref_dec = numpy.append(fixwcs_ref_dec, ref_dec, axis=0)

                if (fixwcs_odi_sourcecat == None):
                    fixwcs_odi_sourcecat = source_cat
                else:
                    #print "Adding some entries to source catalog"
                    fixwcs_odi_sourcecat = numpy.append(fixwcs_odi_sourcecat, source_cat, axis=0)

                fixwcs_bestguess[i,:] = [dx, dy]
            #add_to_bestguess = numpy.array([dx, dy]).reshape((1,2))
            #print fixwcs_bestguess.shape, add_to_bestguess.shape
            #continue
            #fixwcs_bestguess = numpy.append(fixwcs_bestguess, add_to_bestguess, axis=0)

    # Now all processes have returned their results, terminate them 
    # and delete all instances to free up memory
    for cur_process in processes:
        cur_process.terminate()
        del cur_process


    if (fixwcs):
        print fixwcs_bestguess.shape
        print fixwcs_bestguess
    
    #
    # Now do some post-processing:
    # 1) Add or overwrite some headers with values from an external wcs minifits file
    #    to improve the wcs accuracy.
    # 2) Move a couple of headers out of each individual extension and put it in the 
    #    primary extension instead (defined in headers_to_inherit, see podi_definitions)
    # 3) Delete a bunch of headers that are no longer necessary (defined in 
    #    headers_to_delete_from_otas, see podi_definitions)
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
    # Fix the WCS if requested
    #
    if (fixwcs):

        #
        # Combine the most-matches shifts from all available OTAs, 
        # reject outliers to get a more robust result
        #
        wcs_shift_guess = podi_fixwcs.get_overall_best_guess(fixwcs_bestguess)

        #
        # Now go back, match the ODI and reference catalogs and refine the solution
        #
        wcs_shift_refinement, matched_cat, matches = \
            podi_fixwcs.refine_wcs_shift(ref_x=fixwcs_ref_ra, ref_y=fixwcs_ref_dec, 
                                         ota_x=fixwcs_odi_ra, ota_y=fixwcs_odi_dec, 
                                         best_guess=wcs_shift_guess, 
                                         matching_radius=3, #arcsec
                                         alignment_checkfile=None,
                                         verbose=False)

        # Add the previous (best-guess) shift and the new refinement
        wcs_shift = wcs_shift_guess + wcs_shift_refinement

        # Create some plots for WCS diagnosis
        fig = matplotlib.pyplot.figure()
        matches_zeroed = matches - wcs_shift_refinement

        count, xedges, yedges = numpy.histogram2d(matches_zeroed[:,0]*3600., matches_zeroed[:,1]*3600.,
                                                  bins=[60,60], range=[[-3,3], [-3,3]])
        img = matplotlib.pyplot.imshow(count, interpolation='nearest', extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),origin='lower')
        fig.colorbar(img)

        matplotlib.pyplot.plot(matches_zeroed[:,0]*3600., matches_zeroed[:,1]*3600., "b,", linewidth=0)
        matplotlib.pyplot.title("WCS Scatter\n%s" % outputfile)
        matplotlib.pyplot.xlabel("error RA ['']")
        matplotlib.pyplot.ylabel("error DEC ['']")
        matplotlib.pyplot.xlim((-3,3))
        matplotlib.pyplot.ylim((-3,3))
        matplotlib.pyplot.grid(True)
        matplotlib.pyplot.axes().set_aspect('equal')
        png_wcsscatter = outputfile[:-5]+".wcs1.png"
        fig.savefig(png_wcsscatter)
        fig.savefig(png_wcsscatter[:-4]+".eps")
        matplotlib.pyplot.close()

        fig = matplotlib.pyplot.figure()
        matched_zeroed = matched_cat
        matched_zeroed[:,0:2] -= wcs_shift_refinement
        #matplotlib.pyplot.plot(matched_zeroed[:,0], matched_zeroed[:,1], ",", color=(1,1,1))
        matching_radius_arcsec = 3. / 3600.
        valid = (numpy.fabs(matched_zeroed[:,0]-matched_zeroed[:,2]) < matching_radius_arcsec) & \
            (numpy.fabs(matched_zeroed[:,1]-matched_zeroed[:,3]) < matching_radius_arcsec)
        matched_zeroed = matched_zeroed[valid]
        matplotlib.pyplot.quiver(matched_zeroed[:,0], matched_zeroed[:,1],
                                 (matched_zeroed[:,0]-matched_zeroed[:,2]), 
                                 (matched_zeroed[:,1]-matched_zeroed[:,3]),
                                 linewidth=0
                                 )
        # Determine min and max values
        ramin, ramax = numpy.min(matched_zeroed[:,0]), numpy.max(matched_zeroed[:,0])
        decmin, decmax = numpy.min(matched_zeroed[:,1]), numpy.max(matched_zeroed[:,1])
        matplotlib.pyplot.title("WCS misalignment\n%s" % outputfile)
        matplotlib.pyplot.xlim((ramin-0.02, ramax+0.02))
        matplotlib.pyplot.ylim((decmin-0.02, decmax+0.02))
        matplotlib.pyplot.xlabel("RA [degrees]")
        matplotlib.pyplot.xlabel("DEC [degrees]")
        png_wcsdirection = outputfile[:-5]+".wcs2.png"
        fig.savefig(png_wcsdirection)
        fig.savefig(png_wcsdirection[:-4]+".eps")
        #fig.show(block=True)
        matplotlib.pyplot.close()
        


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
            podi_fixwcs.apply_wcs_shift(wcs_shift, ota_list[extension].header)

        # Save the two catalogs to the output file
        ota_list.append( sexcat_to_tablehdu(fixwcs_ref_ra, fixwcs_ref_dec) )

        ota_list.append( odi_sources_to_tablehdu(ota_list, fixwcs_odi_sourcecat) )

        print >>x, "\n\n\n\n\n"
        dummy = numpy.ones(shape=(fixwcs_odi_ra.shape[0],2))
        dummy[:,0] = fixwcs_odi_ra[:] + wcs_shift[0]
        dummy[:,1] = fixwcs_odi_dec[:] + wcs_shift[1]
        numpy.savetxt(x, dummy)

        print >>x, "\n\n\n\n\n"
        numpy.savetxt(x, matched_cat)

        print >>x, "\n\n\n\n\n"
        numpy.savetxt(x, fixwcs_odi_sourcecat)

        x.close()

    hdulist = pyfits.HDUList(ota_list)
    if (not batchmode):
        stdout_write(" writing ...")
        clobberfile(outputfile)
        hdulist.writeto(outputfile, clobber=True)
    else:
        stdout_write(" continuing ...")
        return hdulist

    stdout_write(" done!\n")
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

    # First of all update all Ra/Dec values
    final_cat = None
    for ext in range(1, len(ota_list)):
        extname = ota_list[ext].header['EXTNAME']
        if (extname[0:3] == "OTA" and extname[5:] == ".SCI"):
            ota = int(ota_list[ext].header['EXTNAME'][3:5])
            wcs = astWCS.WCS(ota_list[ext].header, mode="pyfits")

            in_this_ota = fixwcs_odi_sourcecat[:,12] == ota
            sources_ota = fixwcs_odi_sourcecat[in_this_ota]

            xy = sources_ota[:,2:4]
            
            sources_ota[:,0:2] = odi2wcs(xy, wcs)
            
            if (final_cat == None):
                final_cat = sources_ota
            else:
                final_cat = numpy.append(final_cat, sources_ota, axis=0)

    # source_info = [ ra, dec, center_x, center_y, fwhm_x, fwhm_y, amplitude, peak, 
    #            bg_level, bg_variance, s_n, area, extension_id]

    columns = [\
        pyfits.Column(name='Ra',                 format='E', array=final_cat[:, 0]),
        pyfits.Column(name='Dec',                format='E', array=final_cat[:, 1]),
        pyfits.Column(name='X',                  format='E', array=final_cat[:, 2]),
        pyfits.Column(name='Y',                  format='E', array=final_cat[:, 3]),
        pyfits.Column(name='FWHM_X',             format='E', array=final_cat[:, 4]),
        pyfits.Column(name='FWHM_Y',             format='E', array=final_cat[:, 5]),
        pyfits.Column(name='Amplitude',          format='E', array=final_cat[:, 6]),
        pyfits.Column(name='Peak',               format='E', array=final_cat[:, 7]),
        pyfits.Column(name='Background',         format='E', array=final_cat[:, 8]),
        pyfits.Column(name='BackgroundNoise',    format='E', array=final_cat[:, 9]),
        pyfits.Column(name='SignalToNoise',      format='E', array=final_cat[:,10]),
        pyfits.Column(name='NSignificantPixels', format='E', array=final_cat[:,11]),
        pyfits.Column(name='OTA',                format='E', array=final_cat[:,12]),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.update_ext_name("CAT.ODI", comment=None)
    return tbhdu



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

    # For now assume that the WCS template file is located in the same directory as the executable
    root_dir, py_exe = os.path.split(os.path.abspath(sys.argv[0]))
    wcs_solution = root_dir + "/wcs_distort2.fits"
    wcs_solution = cmdline_arg_set_or_default("-wcs", wcs_solution)
    if (not os.path.isfile(wcs_solution)):
        wcs_solution = None

    fixwcs = cmdline_arg_isset("-fixwcs")
    
    hardcoded_detsec = cmdline_arg_isset("-hard_detsec")
    clobber_mode = not cmdline_arg_isset("-noclobber")

    # Handle all reduction flags from command line
    bias_dir, dark_dir, flatfield_dir, bpm_dir, start = read_reduction_directories()
    
    # Collect all cells, perform reduction and write result file
    try:
        collectcells(input, outputfile,
                     bias_dir, dark_dir, flatfield_dir, bpm_dir,
                     wcs_solution=wcs_solution,
                     fixwcs=fixwcs,
                     hardcoded_detsec=hardcoded_detsec,
                     clobber_mode=clobber_mode)
    except:
        stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
        etype, error, stackpos = sys.exc_info()
        stdout_write("# Exception report:")
        stdout_write("#  ==> %s\n" % (error))
        stdout_write("#\n##############################\n")

