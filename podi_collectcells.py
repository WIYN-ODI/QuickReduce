#!/usr/bin/env python

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
import ephem

gain_correct_frames = False
from podi_definitions import *

def break_region_string(str_region):
    reg = str_region[1:-1]
    x,dummy,y = reg.partition(",")
    x1,dummy,x2 = x.partition(":")
    y1,dummy,y2 = y.partition(":")
    return int(x1)-1, int(x2)-1, int(y1)-1, int(y2)-1

def extract_region(data, str_region):
    x1,x2,y1,y2 = break_region_string(str_region)
    return data[y1:y2+1, x1:x2+1]


def insert_into_array(data, from_region, target, target_region):

    fx1, fx2, fy1, fy2 = break_region_string(from_region)
    tx1, tx2, ty1, ty2 = break_region_string(target_region)

    if (fx2-fx1 != tx2-tx1 or fy2-fy1 != ty2-ty1):
        print "Dimensions do not match, doing nothing"
    else:
        target[ty1:ty2+1, tx1:tx2+1] = data[fy1:fy2+1, fx1:fx2+1]

    return 0

def mask_broken_regions(datablock, regionfile, verbose=False):

    counter = 0
    file = open(regionfile)
    for line in file:
        if (line[0:3] == "box"):
            coords = line[4:-2]
            coord_list = coords.split(",")
                        
            if (not datablock == None):
                x, y = int(float(coord_list[0])), int(float(coord_list[1]))
                dx, dy = int(0.5*float(coord_list[2])), int(0.5*float(coord_list[3]))
                #mask[y-dy:y+dy,x-dx:x+dx] = 1

                x1 = numpy.max([0, x-dx])
                x2 = numpy.min([datablock.shape[1], x+dx])
                y1 = numpy.max([0, y-dy])
                y2 = numpy.min([datablock.shape[0], y+dy])
                datablock[y1:y2, x1:x2] = numpy.NaN

                # print x,x+dx,y,y+dy
            counter += 1

    file.close()
    if (verbose):
        print "Marked",counter,"bad pixel regions"
    return datablock

def read_reduction_directories(start=1, warn=True):
    #
    # Read other parameters, specifying the directories for the 
    # flatfields, darks and biases
    #
    # Set all reduction folder to None to mask them as not set
    flatfield_dir = None
    bias_dir = None
    dark_dir = None
    bpm_dir = None
    #
    # Then deal with the user-specified values from the command line
    #
    i=start
    while i<len(sys.argv):
        cmd = sys.argv[i]
        if (cmd == "-flat" and i+1<len(sys.argv)):
            flatfield_dir = sys.argv[i+1]
            i += 1
        elif (cmd == "-bias" and i+1<len(sys.argv)):
            bias_dir = sys.argv[i+1]
            i += 1
        elif (cmd == "-dark" and i+1<len(sys.argv)):
            dark_dir = sys.argv[i+1]
            i += 1
        elif (cmd == "-bpm" and i+1<len(sys.argv)):
            bpm_dir = sys.argv[i+1]
            i += 1
            if (bpm_dir == "auto"):
                bpm_dir, py = os.path.split(sys.argv[0])
                if (bpm_dir == ""):
                    bpm_dir = "."
        else:
            if (warn):
                print "Don't understand parameter %d: %s" % (i, sys.argv[i])
            break
        i += 1


    # Output some summary on the reduction
    print """
Calibration data:
            Bias: %s
            Dark: %s
      Flatfields: %s
  Bad pixel mask: %s
""" % (bias_dir, dark_dir, flatfield_dir, bpm_dir)

    return bias_dir, dark_dir, flatfield_dir, bpm_dir, i

def collect_reduce_ota(filename,
                       bias_dir, dark_dir, flatfield_dir, bpm_dir):

    # Create an fits extension to hold the output
    hdu = pyfits.ImageHDU()

    if (not os.path.isfile(filename)):
        stdout_write("Couldn't find file %s ..." % (filename))
    else:
        hdulist = pyfits.open(filename)

        detsize = break_region_string(hdulist[0].header['DETSIZE'])
        det_x1, det_x2, det_y1, det_y2 = detsize
        #print det_x1, det_x2, det_y1, det_y2

        size_x, size_y = det_x2 - det_x1 + 1, det_y2 - det_y1 + 1
        #print size_x, size_y
        size_x, size_y = 4096, 4096
        #print size_x, size_y

        obsid = hdulist[0].header["OBSID"]
        ota = int(hdulist[0].header['FPPOS'][2:])
        try:
            ota_id = all_otas.index(ota)
        except:
            stdout_write("Something is wrong with this OTA, it's not among those listed as available")
            sys.exit(-1)
        ota_c_x, ota_c_y = available_ota_coords[ota_id]

	merged = numpy.ones(shape=(size_x, size_y), dtype=numpy.float32) * -1e8
        
        for cell in range(1,65):
            stdout_write("\r%s:   OTA %02d, cell %s ..." % (obsid, ota, hdulist[cell].header['EXTNAME']))

            # Check if this is one of the broken cells
            wm_cellx, wm_celly = hdulist[cell].header['WN_CELLX'], hdulist[cell].header['WN_CELLY']
            broken = False
            list_of_broken_cells = broken_ota_cells[ota_id]
            for broken_cell in list_of_broken_cells:
                x,y = broken_cell
                #print x,y
                if (wm_cellx == x and wm_celly == y):
                    broken = True
                    #print "found broken cell", hdulist[cell].header['EXTNAME'],broken_cell
                    break

            # If not, overscan subtract and insert into large frame
            if (not broken):
                overscan_region = extract_region(hdulist[cell].data, hdulist[cell].header['BIASSEC'])
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
                
        # Set all unset pixels to NaN 
        merged[merged < -1e5] = numpy.NaN

        #
        # Get some information for the OTA
        #
        fppos = hdulist[0].header['FPPOS']
        filter_name = hdulist[0].header['FILTER']
        exposure_time = hdulist[0].header['EXPTIME']

	# If we are to do some bias subtraction:
	if (not bias_dir == None):
       	    bias_filename = "%s/bias.fits" % (bias_dir)
	    if (os.path.isfile(bias_filename)):
                bias = pyfits.open(bias_filename)
	   	# Search for the flatfield data for the current OTA
            	for bias_ext in bias[1:]:
                    fppos_bias = bias_ext.header['FPPOS']

            	    if (fppos_bias == fppos):
                    	# This is the one
                    	merged -= bias_ext.data
                	break

                bias.close()
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
	   	# Search for the flatfield data for the current OTA
            	for dark_ext in dark[1:]:
                    fppos_dark = dark_ext.header['FPPOS']

            	    if (fppos_dark == fppos):
                        darktime = dark_ext.header['EXPTIME']
                    	# This is the one
                    	merged -= (dark_ext.data * exposure_time / darktime)
                	break

                dark.close()
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
                del flatfield

        # Finally, apply bad pixel masks 
        # Determine which region file we need
        if (not bpm_dir == None):
            region_file = "%s/bpm_%s.reg" % (bpm_dir, fppos)
            if (os.path.isfile(region_file)):
                # Apply the bad pixel regions to file, marking
                # all bad pixels as NaNs
                mask_broken_regions(merged, region_file)

        # Now copy the headers from the original file into the new one
        cards = hdulist[0].header.ascardlist()
        for c in cards:
            hdu.header.update(c.key, c.value, c.comment)

        # Insert the DETSEC header so IRAF understands where to put the extensions
	start_x = ota_c_x * 4100
	start_y = ota_c_y * 4100        
	end_x = start_x + det_x2 - det_x1
	end_y = start_y + det_y2 - det_y1
	detsec_str = "[%d:%d,%d:%d]" % (start_x, end_x, start_y, end_y)
	hdu.header.update("DETSEC", detsec_str, "position of OTA in focal plane")
                
        #
        # Fudge with the WCS headers, largely undoing what's in the fits file right now,
        # and replacing it with a simpler version that hopefully works better
        #
        hdu.header['CTYPE1'] = "RA---TAN"
        hdu.header['CTYPE2'] = "DEC--TAN"
        hdu.header.update("CUNIT1", "deg", "")
        hdu.header.update("CUNIT2", "deg", "")
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
        
        #print hdu.header['RA'], hdu.header['DEC'], hdu.header['EQUINOX']
        #coord = ephem.Equatorial(hdu.header['RA'], hdu.header['DEC'], epoch=str(hdu.header['EQUINOX']))
        #coord_j2000 = ephem.Equatorial(coord, epoch=ephem.J2000)
        coord_j2000 = ephem.Equatorial(hdu.header['RA'], hdu.header['DEC'], epoch=ephem.J2000)
        
        # Add some offsets, again from the PPA pipeline, to correct the pODI pointing
        #hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  #+ 0.2630/numpy.cos(coord_j2000.dec)
        #hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) #- 0.0435

        # This works for NGC2146
        #hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  + 0.1592/numpy.cos(coord_j2000.dec)
        #hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) - 0.0386

        hdu.header.update('XXVAL1', numpy.degrees(coord_j2000.ra),  "CRPIX1 w/o correction")
        hdu.header.update('XXVAL2', numpy.degrees(coord_j2000.dec), "CRPIX2 w/o correction")
        
        del hdu.header['EQUINOX']

        # THis works for NGC2146
        #hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  
        #hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) 
        #print hdu.header['CRPIX1'], hdu.header['CRPIX2']
        #hdu.header['CRPIX1'] -= 5123
        #hdu.header['CRPIX2'] += 1254

        # This works for the Crab Nebula
        hdu.header['CRVAL1'] = numpy.degrees(coord_j2000.ra)  
        hdu.header['CRVAL2'] = numpy.degrees(coord_j2000.dec) 

        # Also update the CD matric with data from the official pipeline to
        # account for small misalignments between the different OTAs
        hdu.header['CD1_1'] = wcs_cd[ota_id][0]
        hdu.header['CD2_2'] = wcs_cd[ota_id][1]
        hdu.header['CD1_2'] = wcs_cd[ota_id][2]
        hdu.header['CD2_1'] = wcs_cd[ota_id][3]

        # Insert the new image data. This also makes sure that the headers
        # NAXIS, NAXIS1, NAXIS2 are set correctly
        hdu.data = merged
    return hdu
    

def collectcells(input, outputfile,
                 bias_dir, dark_dir, flatfield_dir, bpm_dir):

    # As a safety precaution, if the first parameter is the directory containing 
    # the files, extract just the ID string to be used for this script
    if (input[-1] == "/"):
	input = input[:-1]

    directory,basename = os.path.split(input)
    if (directory == ""):
        directory = "."
    print "Merging cells for frame %s" % (basename)

    if (outputfile == None):
        outputfile = "%s/%s.fits" % (directory, basename)

    # Check if the user requested us to prepare the frame for SExtractor
    prep_for_sextractor = cmdline_arg_isset("-prep4sex")

    # Start new list of HDUs
    ota_list = []

    # And add the primary HDU to make the fits file a valid one
    primhdu = pyfits.PrimaryHDU()
    ota_list.append(primhdu)
    
    for ota_id in range(len(available_ota_coords)):
        ota_c_x, ota_c_y = available_ota_coords[ota_id]        
        ota = ota_c_x * 10 + ota_c_y

        filename = "%s/%s/%s.%02d.fits" % (directory, basename, basename, ota)
                
        hdu = collect_reduce_ota(filename,
                                 bias_dir, dark_dir, flatfield_dir, bpm_dir)

        # SExtractor doesn't like NaNs, so replace all of them with something
        # more negative than -1e30 (that's -1 times SEx's BIG variable)
        if (prep_for_sextractor):
            hdu.data[numpy.isnan(hdu.data)] = -1e31 

        # Insert into the list to be used later
        ota_list.append(hdu)

    stdout_write(" writing ...")
    
    hdulist = pyfits.HDUList(ota_list)
    clobberfile(outputfile)
    hdulist.writeto(outputfile, clobber=True)

    stdout_write(" done!\n")
    return 0

if __name__ == "__main__":

    # Read the input directory that contains the individual OTA files
    input = sys.argv[1]

    # Assign a fallback output filename if none is given 
    if (len(sys.argv)>2):
        outputfile = sys.argv[2]
    else:
        print "No output filename has been given, setting to default mergedcells.fits"
        outputfile = "mergedcells.fits"
    print "Writing results into",outputfile

    # Handle all reduction flags from command line
    bias_dir, dark_dir, flatfield_dir, bpm_dir, start = read_reduction_directories(start=3)

    # Collect all cells, perform reduction and write result file
    collectcells(input, outputfile,
                 bias_dir, dark_dir, flatfield_dir, bpm_dir)
    
