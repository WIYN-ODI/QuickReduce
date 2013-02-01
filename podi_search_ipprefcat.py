#!/usr/bin/env python


import sys
import numpy
import os
import math

from podi_definitions import *
import pyfits

IPP_DIR = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

def get_reference_catalog(ra, dec, radius, basedir):

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (basedir)
    skytable_hdu = pyfits.open(skytable_filename)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data
    #print skytable[:3]
    
    # Select entries that match our list
    print "# Searching for stars within %.1f degress around %f , %f ..." % (radius, ra, dec)

    min_dec = dec - radius
    max_dec = dec + radius
    min_ra = ra - radius/math.cos(math.radians(dec))
    max_ra = ra + radius/math.cos(math.radians(dec))

    if (max_ra > 360.):
        # This wraps around the high end, shift all ra values by -180
        # Now all search RAs are ok and around the 180, next also move the catalog values
        selected = skytable['R_MIN'] < 180
        skytable['R_MAX'][selected] += 360
        skytable['R_MIN'][selected] += 360
    if (min_ra < 0):
        # Wrap around at the low end
        selected = skytable['R_MAX'] > 180
        skytable['R_MAX'][selected] -= 360
        skytable['R_MIN'][selected] -= 360

    print "# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec)
    
    needed_catalogs = (skytable['PARENT'] > 0) & (skytable['PARENT'] < 25) & \
         (skytable['R_MAX'] > min_ra)  & (skytable['R_MIN'] < max_ra) & \
         (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)

    #print skytable[needed_catalogs]
    
    files_to_read = skytable['NAME'][needed_catalogs]

    # Now quickly go over the list and take care of all filenames that still have a 0x00 in them
    for i in range(len(files_to_read)):
        found_at = files_to_read[i].find('\0')
        if (found_at > 0):
            files_to_read[i] = files_to_read[i][:found_at]
        
    # Now we are with the skytable catalog, so close it
    skytable_hdu.close()
    del skytable

    #print files_to_read

    # Load all frames, one by one, and select all stars in the valid range.
    # Then add them to the catalog with RAs and DECs
    full_catalog = numpy.zeros(shape=(0,6))
    for catalogname in files_to_read:
        catalogfile = "%s/%s.cpt" % (basedir, catalogname)

        hdu_cat = pyfits.open(catalogfile)
        # hdu_cat.info()

        # Read the RA and DEC values
        cat_ra  = hdu_cat['DVO_AVERAGE_ELIXIR'].data.field('RA')
        cat_dec = hdu_cat['DVO_AVERAGE_ELIXIR'].data.field('DEC')

        photom_raw = hdu_cat['DVO_MEASURE_ELIXIR'].data.field('MAG')
        photom_grizy = photom_raw.reshape((cat_ra.shape[0], 5))

        # To slect the right region, shift a temporary catalog
        cat_ra_shifted = cat_ra
        if (max_ra > 360.):
            cat_ra_shifted[cat_ra < 180] += 360
        elif (min_ra < 0):
            cat_ra_shifted[cat_ra > 180] -= 360

        select_from_cat = (cat_ra_shifted > min_ra) & (cat_ra_shifted < max_ra ) & (cat_dec > min_dec) & (cat_dec < max_dec)

        array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),6))
        array_to_add[:,0] = cat_ra[select_from_cat]
        array_to_add[:,1] = cat_dec[select_from_cat]

        array_to_add[:,2:6] = photom_grizy[:,0:4][select_from_cat] / 1e3
        
        #print cat_ra.shape, cat_dec.shape, photom_grizy.shape
        #print array_to_add.shape

        print "# Read %d stars from catalog %s ..." % (array_to_add.shape[0], catalogfile)
        
        #full_catalog = numpy.concatenate((full_catalog, array_to_add), axis=0)
        full_catalog = numpy.append(full_catalog, array_to_add, axis=0)
        #print photom_grizy[:3,:]
        
    print "# Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read))
    return full_catalog

    
if __name__ == "__main__":
    
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    radius = float(sys.argv[3])

    basedir = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

    catalog = get_reference_catalog(ra, dec, radius, basedir)

    if (True): #False):
        for i in range(catalog.shape[0]):
            for j in range(catalog.shape[1]):
                print catalog[i,j],
            print

