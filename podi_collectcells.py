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
                       pixelvalue_indef=numpy.NaN):

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

        # Now copy the headers from the original file into the new one
        cards = hdulist[0].header.ascardlist()
        for c in cards:
            hdu.header.update(c.key, c.value, c.comment)

        ota_c_x, ota_c_y = available_ota_coords[ota_id]

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
                
        if (cmdline_arg_isset("-simplewcs") or cmdline_arg_isset("-scamp")):
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
        hdu.header['CRVAL1'] += offset_total[0] / 3600. / numpy.cos(numpy.radians(hdu.header['CRVAL2']))
        hdu.header['CRVAL2'] += offset_total[1] / 3600.
        #
        # To do:
        # =========================================================
        # Check if the above still makes sense !!!!
        # In particular the addition of the telescope offsets 
        # should be included in RA/DEC already !!!
        # =========================================================
        #
        
        # Insert the new image data. This also makes sure that the headers
        # NAXIS, NAXIS1, NAXIS2 are set correctly
        hdu.data = merged
    return hdu
    

def read_scamp_header(filename, dump_header=False):

    if (not os.path.isfile(filename)):
        return None

    headfile = open(filename, "r")
    lines = headfile.readlines()

    values = {}
    values_list = []
    comments = {}
    comments_list = []

    for line in lines:
        #print line

        key, value, comment = line[0:8].strip(), line[9:30].strip(), line[32:].strip()
        if (key in ("HISTORY", "COMMENT",
                    ) ):
            # Don't know what to do with those, so skip'em
            continue

        elif (key in ("FGROUPNO", "FLXSCALE", "MAGZEROP", 
                      "ASTINST",
                      "PHOTIRMS", "PHOTINST", "PHOTLINK",
                      ) ):
            # These are some scamp-specific headers, let's not copy them
            continue

        elif (key in ("CRVAL1", "CRVAL2",
                      "CRPIX1", "CRPIX2", 
                      "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                      "PV1_0", "PV1_1", "PV1_2", "PV1_4", "PV1_5", "PV1_6",
                      "PV2_0", "PV2_1", "PV2_2", "PV2_4", "PV2_5", "PV2_6",
                      "EQUINOX",
                      ) ):
            value = float(value)

        elif (key in ("RADECSYS", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2") ):
            # Strip the unnecessary quotes and spaces from the end
            value = value[1:-1].strip()

        elif (key == "END"):
            # This concludes one extension, add it to list and start new 
            # list for the next OTA
            values_list.append(values)
            comments_list.append(comments)
            values = {}
            comments = {}
            continue

        values[key] = value
        comments[key] = comment

    #print head_list
    #if (dump_header):
    #    for ota in head_list:
    #        for key,value,comment in ota:
    #            print key,"=",value,"(",comment,")"
    #        print "\n\n\n"

    return values_list, comments_list

def apply_scamp_solution(scampfile, ota_list):

    if (not os.path.isfile(scampfile)):
        stdout_write("Can not find scamp header file %s,skipping ...\n" % scampfile)
        return

    scamp_value, scamp_comment = read_scamp_header(scampfile)

    # First of all, figure out the best offset we need to apply 
    # to the CRVAL headers from the scamp solution

    # Important: The solution from scamp has to have 13 blocks for 13 OTAs.
    if (len(scamp_value) != len(available_otas)):
        stdout_write("Illegal scamp solution!\n")
        return -1

    # Ok, now we know we have a valid solution
    ref_crval1, ref_crval2 = scamp_value[6]["CRVAL1"], scamp_value[6]["CRVAL2"]
    #print "\nREF=",ref_crval1,ref_crval2

    # Also compute where the telescope if pointing at, i.e. what CRVAL should be
    target_crval1, target_crval2 = ota_list[6].header['CRVAL1'], ota_list[6].header['CRVAL2']

    # Adjust the solution for the known shift in CRPIX position
    # In pODI: True pointing of telescope is ~ at pixel 4200,4200 of OTA 3,3
    # From scamp: Typically somewhere around 2000,2000 in OTA 3,3
    # NB: 0.11 / 3600. is the pixel scale in degrees/pixel, ignoring rotation and distortion
    dx_crpix = 4200 - scamp_value[6]["CRPIX1"]
    dy_crpix = 4200 - scamp_value[6]["CRPIX2"]
    declination = target_crval2
    d_ra  = dx_crpix * 0.11 / 3600. / math.cos(math.radians(declination))
    d_dec = dy_crpix * 0.11 / 3600.
    target_crval1 -= d_ra
    target_crval2 -= d_dec

    # With this data at hand, work out the shift we need to apply to the scamp solution
    for ext in range(1,len(ota_list)):
        ota = ext - 1

        scamp_value[ota]['CRVAL1'] += target_crval1 - ref_crval1
        scamp_value[ota]['CRVAL2'] += target_crval2 - ref_crval2

        # As a last step, simple copy all headers to the ota_list
        for key in scamp_value[ota]:
            
            # If the header already exists, simple change the value
            if (key in ota_list[ext].header):
                ota_list[ext].header[key] = scamp_value[ota][key]
                
            # If the header doesn't exist yet, create it and use the comment from SCAMP
            else:
                ota_list[ext].header.update(key, scamp_value[ota][key], scamp_comment[ota][key])

    # That's it, all done!
    return 0


def collectcells(input, outputfile,
                 bias_dir=None, dark_dir=None, flatfield_dir=None, bpm_dir=None,
                 batchmode=False):

    if (os.path.isfile(input)):
        # Assume this is one of the fits files in the right directory
        # In that case, extract the FILENAME header and convert it into 
        # the filebase we need to construct the filenames of all OTA fits files.

        hdulist = pyfits.open(input)
        filebase = hdulist[0].header['FILENAME'][:18]
        hdulist.close()
        del hdulist

        # Split the input filename to extract the directory part
        directory,dummy = os.path.split(input)

    elif (os.path.isdir(input)):
        # As a safety precaution, if the first parameter is the directory containing 
        # the files, extract just the ID string to be used for this script
        if (input[-1] == "/"):
            input = input[:-1]

        basedir, filebase = os.path.split(input)
        directory = input
        

    #print "Merging cells for frame %s" % (basename)

    if (outputfile == None):
        outputfile = "%s/%s.fits" % (directory, filebase)

    if (outputfile.find("%") >= 0):
        # The output filename contains special tags that should 
        # be replaced by values from the file header

        filename = "%s/%s.%02d.fits" % (directory, filebase, 33)
        hdulist = pyfits.open(filename)
        header = hdulist[0].header
        
        while (outputfile.find("%") >= 0):
            start = outputfile.find("%") 
            if (outputfile[start:start+7] == "%FILTER"):
                outputfile = outputfile[:start] + header['FILTER'] + outputfile[start+7:]
            elif (outputfile[start:start+7] == "%OBJECT"):
                # The object name might contain spaces, replace them with underscores
                objectname = header['OBJECT'].replace(' ', '_')
                outputfile = outputfile[:start] + objectname  + outputfile[start+7:]
            elif (outputfile[start:start+6] == "%OBSID"):
                outputfile = outputfile[:start] + header['OBSID'] + outputfile[start+6:]
            else:
                stdout_write("found unknown tag in %s\n" % outputfile)
                break

        hdulist.close()
        del hdulist
        del header

        stdout_write("Replaced some keywords, new output filename:\n ---> %s\n" % (outputfile))


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
    ota_list = []

    # And add the primary HDU to make the fits file a valid one
    primhdu = pyfits.PrimaryHDU()
    ota_list.append(primhdu)
    
    # Set the fallback value for undefined pixels (mostly the gaps between the OTA cells)
    # This works perfectly fine in ds9, but not SExtractor
    pixelvalue_indef = numpy.NaN
    if (cmdline_arg_isset("-prep4sex")):
        # Check if the user requested us to prepare the frame for SExtractor
        # SExtractor doesn't like NaNs, so replace all of them with something
        # more negative than -1e30 (that's -1 times SEx's BIG variable)
        pixelvalue_indef = -1e31
    
    for ota_id in range(len(available_ota_coords)):
        ota_c_x, ota_c_y = available_ota_coords[ota_id]        
        ota = ota_c_x * 10 + ota_c_y

        if (cmdline_arg_isset('-singleota')):
            single_ota = int(get_cmdline_arg("-singleota"))
            if (ota != single_ota):
                continue

        filename = "%s/%s.%02d.fits" % (directory, filebase, ota)
                
        hdu = collect_reduce_ota(filename,
                                 bias_dir, dark_dir, flatfield_dir, bpm_dir,
                                 offset_pointing=offset_pointing,
                                 offset_dither=offset_dither,
                                 target_coords=target_coords,
                                 pixelvalue_indef=pixelvalue_indef)

        # Insert into the list to be used later
        ota_list.append(hdu)

    #
    # Now do some post-processing:
    # 1) Add or overwrite some headers with values from an external scamp header file
    #    to improve the wcs accuracy.
    # 2) Move a couple of headers out of each individual extension and put it in the 
    #    primary extension instead (defined in headers_to_inherit, see podi_definitions)
    # 3) Delete a bunch of headers that are no longer necessary (defined in 
    #    headers_to_delete_from_otas, see podi_definitions)
    #

    # If the user specified to overwrite the WCS with a SCAMP solution,
    # Read the solution and store it for later use
    scamp_solution = cmdline_arg_set_or_default("-scamp", None)
    scamp_header = None
    if (not scamp_solution == None and not cmdline_arg_isset("-singleota")):
        apply_scamp_solution(scamp_solution, ota_list)

    # Now update the headers in all OTA extensions.
    for extension in range(1, len(ota_list)):
        ota = ota_list[extension]

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
    
