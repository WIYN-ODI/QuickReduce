#!/usr/bin/env python


import sys
import numpy
import os
from query_usno import query_usno
from podi_definitions import *
import pyfits
#import date
import datetime
#import pywcs  
#from astLib import astWCS
import pdb
import scipy
import scipy.stats


arcsec = 1./3600.
number_bright_stars = 100
max_offset = 0.1


N = 200
N_brightest_ota = 70
N_brightest_ref = 150

# Compute matches in smaller blocks to not blow up memory
# Improve: Change execution to parallel !!!
blocksize = 100



def photcalib(fitsfile, output_filename, calib_directory):

    tmp, dummy = os.path.split(sys.argv[0])
    dot_config_dir = tmp + "/.config/"
    print dot_config_dir
    
    hdulist = pyfits.open(fitsfile)
    ra, dec = hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

    #
    # Run Sextractor on the input frame
    #
    sex_catalogfile = fitsfile[:-5]+".photcalib.cat"
    sex_config_file = dot_config_dir+"fixwcs.conf"
    sex_param_file = dot_config_dir+"fixwcs.param"
    sex_logfile = fitsfile[:-5]+".sextractor.log"
    sex_cmd = "sex -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s >& %s" % (sex_config_file, sex_param_file, sex_catalogfile, fitsfile, sex_logfile)
    print sex_cmd
    stdout_write("Running SExtractor to search for stars, be patient (logfile: %s) ..." % (sex_logfile))
    if ((not cmdline_arg_isset("-skip_sex")) or (not os.path.isfile(sex_catalogfile))):
        os.system(sex_cmd)
    else:
        stdout_write("Re-using previous source catalog\n")
    stdout_write(" done!\n")

    stdout_write("\nPreparing work ...\n\n")
    
    # Read the Sextractor output catalog
    sex_cat = numpy.loadtxt(sex_catalogfile)
    stdout_write("Reading %d stars from Sextractor catalog\n" % sex_cat.shape[0])

    # Eliminate all stars with flags
    sex_cat = sex_cat[sex_cat[:,10] == 0]
    stdout_write("%d stars left after eliminating flags\n" % sex_cat.shape[0])


    # Figure out which SDSS to use for calibration
    filter = hdulist[0].header['FILTER']
    sdss_filter = sdss_equivalents[filter]

    # Load the standard star catalog
    sdss_cat_filename = "%s/stripe82__%s.cat.npy" % (calib_directory, sdss_filter)
    stdout_write("Loading standard stars (%s)..." % (sdss_cat_filename))
    sdss_cat = numpy.load(sdss_cat_filename)
    stdout_write(" %d stars loaded!\n" % (sdss_cat.shape[0]))

    # Now take take of pointings close to RA of 0 hours
    if (ra < 5):
        stdout_write("Close to RA=0h, wrapping coordinates ...\n")
        sdss_cat[:,0][sdss_cat[:,0] > 180] -= 360.
        sex_cat[:,6][sex_cat[:,6] > 180] -= 360.

    # Break down the full catalog to stars nearby
    nearby = (sdss_cat[:,0] > ra-1) & (sdss_cat[:,0] < ra+1) & (sdss_cat[:,1] > dec-1) & (sdss_cat[:,1] < dec+1) 
    std_stars = sdss_cat[nearby]
    stdout_write("Found %d nearby (RA=%.1f, DEC=%.1f, +/- 1deg) standard stars!\n" % (std_stars.shape[0], ra, dec))
    dump = open("dump.cat", "w")
    for i in range(std_stars.shape[0]):
        print >>dump, std_stars[i,0], std_stars[i,1]
    print >>dump, "\n\n\n\n\n\n\n\n"
    for i in range(sex_cat.shape[0]):
        print >>dump, sex_cat[i,6], sex_cat[i,7]
    
    dump.close()
    
    #
    # Now go through each of the extension
    # Improve: Change execution to parallel !!!
    #
    stdout_write("\nStarting work, results in %s ...\n\n" % output_filename)
    results = open(output_filename, "w")
    for ext in range(1, len(hdulist)):

        #
        # Make sure this is one of the OTAs covered by the specified filter
        #
        fppos = int(hdulist[ext].header['FPPOS'][2:4])
        print "fppos=",fppos
        print filter, "-->", which_otas_to_use[filter]
        if (not fppos in otas_for_photometry[filter]):
            stdout_write("Excluding OTA %02d from photometric calibration for filter %s\n" % (fppos, filter))
            continue
        
        #
        # Select a sub-catalog with stars in this OTA
        #
        ota_cat = sex_cat[sex_cat[:,1] == ext]
        stdout_write("%d stars on OTA %02d (extension %d)\n" % (ota_cat.shape[0],fppos,ext))

        # Determine the approximate center position of this OTA
        c_ra  = numpy.median(ota_cat[:,6])
        c_dec = numpy.median(ota_cat[:,7])

        # Compute the F.o.v. of this OTA
        w_ra  = (9./60.) / math.cos(math.radians(c_dec))
        w_dec = (9./60.)
        #print w_ra, w_dec
        
        # Pre-select standard stars in the vicinity
        nearby = (std_stars[:,0] > c_ra-w_ra) & (std_stars[:,0] < c_ra+w_ra) & (std_stars[:,1] > c_dec-w_dec) & (std_stars[:,1] < c_dec+w_dec)
        ref_cat = std_stars[nearby]

        stdout_write("OTA %d: %.4f %.4f -> %d stds nearby\n" % (ext, c_ra, c_dec, ref_cat.shape[0]))

        #
        # Now for each star, find the closest match in the reference catalog
        #
        matching_radius = 2 * arcsec
        matching_radius_squared = matching_radius * matching_radius
        for star in range(ota_cat.shape[0]):

            d_ra  = (std_stars[:,0] - ota_cat[star,6]) * math.cos(math.radians(ota_cat[star,7]))
            d_dec = (std_stars[:,1] - ota_cat[star,7])
            sq_d = d_ra * d_ra + d_dec * d_dec

            si = numpy.argsort(sq_d)
            if (sq_d[si[0]] < matching_radius_squared):
                # The closest match is a valid one, with distance < matching_radius
                print >>results, ota_cat[star,6], ota_cat[star,7], \
                                    std_stars[si[0],0], std_stars[si[0],1], \
                                    ota_cat[star,2], ota_cat[star,3], \
                                    std_stars[si[0],2], std_stars[si[0],4], \
                                    ext, \
                                    ota_cat[star,4], ota_cat[star,5]
                
    results.close()

    
if __name__ == "__main__":

    
    if (cmdline_arg_isset("-multi")):
        calibdir = get_cmdline_arg("-calib")
        for infile in get_clean_cmdline()[1:]:
            output_filename = infile[0:-5]+".photcalib.dat"
            photcalib(infile, output_filename, calibdir)
    else:
        fitsfile = get_clean_cmdline()[1]
        output_filename = get_clean_cmdline()[2]
        calibdir = get_clean_cmdline()[3]

        photcalib(fitsfile, output_filename, calibdir)
