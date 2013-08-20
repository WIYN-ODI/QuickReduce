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


import sys
import numpy
import os
import math

from podi_definitions import *
import pyfits

IPP_DIR = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

def twomass_from_cds(ra, dec, radius, verbose):

    # Download the 2MASS PSC from Vizier
    print "   Downloading 2MASS PSC catalog from Vizier ... (%.4f/%.4f +/- %.1f)" % (ra, dec, radius)
    twomass_cat = "/tmp/2mass.cat"
    cds = "find2mass -c %f %f -rd %f -e b -smJ > %s" % (ra, dec, radius, twomass_cat)
    os.system(cds)

    file = open(twomass_cat)
    catalog = []
    for line in file:
        if (line[0] == "#"):
            continue

        items = line.split("|")

        ra_dec = items[0].split()
        ra, dec = float(ra_dec[0]), float(ra_dec[1])

        jmag = items[2].split()
        hmag = items[3].split()
        kmag = items[4].split()

        mags = [float(jmag[0]), float(hmag[0]), float(kmag[0])]
        
        flags = items[5].split()
        for i in range(3):
            if (not flags[0][i] in ('A', 'B', 'C', 'D')):
                mags[i] = numpy.NaN

        entry = [ra, dec, mags[0], mags[1], mags[2]]
        catalog.append(entry)

    return numpy.array(catalog)


def get_reference_catalog(ra, dec, radius, basedir, cattype="2mass_opt", verbose=False):

    if (cattype == "2mass_web"):
        return twomass_from_cds(ra, dec, radius, verbose)

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (basedir)
    skytable_hdu = pyfits.open(skytable_filename)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data
    #print skytable[:3]
    
    # Select entries that match our list
    if (verbose): print "# Searching for stars within %.1f degress around %f , %f ..." % (radius, ra, dec)

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

    if (verbose): print "# Search radius: RA=%.1f ... %.1f   DEC=%.1f ... %.1f" % (min_ra, max_ra, min_dec, max_dec)
    
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
    full_catalog = None #numpy.zeros(shape=(0,6))
    for catalogname in files_to_read:

        if (cattype == "IPPRef"):
            if (verbose): print "Reading IPPRef catalog"
            catalogfile = "%s/%s.cpt" % (basedir, catalogname)

            hdu_cat = pyfits.open(catalogfile)

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

            if (verbose): print "# Read %d stars from catalog %s ..." % (array_to_add.shape[0], catalogfile)

        elif (cattype in ('2mass_opt', '2mass_nir')):
            if (verbose): print "Reading 2mass catalog"
            catalogfile = "%s/%s.fits" % (basedir, catalogname)

            hdu_cat = pyfits.open(catalogfile)

            # Read the RA and DEC values
            cat_ra  = hdu_cat[1].data.field('ra')
            cat_dec = hdu_cat[1].data.field('dec')

            # To slect the right region, shift a temporary catalog
            cat_ra_shifted = cat_ra
            if (max_ra > 360.):
                cat_ra_shifted[cat_ra < 180] += 360
            elif (min_ra < 0):
                cat_ra_shifted[cat_ra > 180] -= 360

            select_from_cat = (cat_ra_shifted > min_ra) & (cat_ra_shifted < max_ra ) & (cat_dec > min_dec) & (cat_dec < max_dec)


            if (cattype == ('2mass_opt')):
                # select only sources with optical counterparts, i.e. valid b or vr magnitudes
                optical_counterparts = numpy.isfinite(hdu_cat[1].data.field('mag_b')) | numpy.isfinite(hdu_cat[1].data.field('mag_vr'))
                selected = select_from_cat & optical_counterparts

                array_to_add = numpy.zeros(shape=(numpy.sum(selected),4))
                array_to_add[:,0] = cat_ra[selected]
                array_to_add[:,1] = cat_dec[selected]
                array_to_add[:,2] = hdu_cat[1].data.field('mag_b')[selected]
                array_to_add[:,3] = hdu_cat[1].data.field('mag_vr')[selected]

                # make sure that column 3 (which holds the magnitudes used for 
                # alignment in fixwcs) is a valid entry
                array_to_add[:,3][numpy.isnan(array_to_add[:,3])] = array_to_add[:,2]

            else:

                array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),5))
                array_to_add[:,0] = cat_ra[select_from_cat]
                array_to_add[:,1] = cat_dec[select_from_cat]
                array_to_add[:,2] = hdu_cat[1].data.field('mag_j')[select_from_cat]
                array_to_add[:,3] = hdu_cat[1].data.field('mag_h')[select_from_cat]
                array_to_add[:,4] = hdu_cat[1].data.field('mag_k')[select_from_cat]

            if (verbose): print "# Read %d stars from catalog %s ..." % (array_to_add.shape[0], catalogfile)

        if (full_catalog == None):
            full_catalog = array_to_add
        else:
            full_catalog = numpy.append(full_catalog, array_to_add, axis=0)
        #print photom_grizy[:3,:]
        
    if (verbose): print "# Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read))
    return full_catalog





if __name__ == "__main__":
    
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    radius = float(sys.argv[3])

    basedir = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

    basedir = cmdline_arg_set_or_default("-basedir", "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/")
    catalog_type = cmdline_arg_set_or_default("-cattype", "IPPRef")


    catalog = get_reference_catalog(ra, dec, radius, basedir=basedir, cattype=catalog_type)

    if (True): #False):
        for i in range(catalog.shape[0]):
            for j in range(catalog.shape[1]):
                print catalog[i,j],
            print

