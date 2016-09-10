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

import podi_logging
import logging

from podi_definitions import *
from podi_commandline import *
import pyfits

IPP_DIR = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

def twomass_from_cds(ra, dec, radius, verbose):

    # Download the 2MASS PSC from Vizier
    print "   Downloading 2MASS PSC catalog from Vizier ... (%.4f/%.4f +/- %.1f)" % (ra, dec, radius)
    twomass_cat = "/tmp/2mass.cat"
    cds = "find2mass -c %f %f -bd %f -e b -smJ > %s" % (ra, dec, radius, twomass_cat)
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




def load_catalog_from_sdss(ra, dec, sdss_filter, verbose=False, return_query=False, max_catsize=-1):
    """
    Query the SDSS online and return a catalog of sources within the specified 
    area

    Parameters
    ----------
    ra : float[2] (e.g. [165.5, 165.9]

        Min and max values for RA

    dec : float[2]

        Same as ra, just for declination

    sdss_filter : string (allowed are u,g,r,i,z)

        Name of the SDSS filter for which to return magnitudes

    verbose : Bool

        Add debugging output

    return_query : Bool

        If set to false, the return value also contains the SQL query used to 
        query the SDSS catalog

    max_catsize : int

        Maximum number of SDSS stars to be returned

    Returns
    -------
        The SDSS catalog, in the format

        Ra, Dec, mag_u, magerr_u, mag_g, magerr_g, ..r, ..i, ..z

    """


    #import sqlcl
    logger = logging.getLogger("GetCatalogFromSDSS")

    #ra = 0
    #print "# Loading catalog from SDSS online..."
    
    if (numpy.array(ra).ndim > 0):
        min_ra = ra[0]
        max_ra = ra[1]
    else:
        min_ra = ra - 0.6/math.cos(math.radians(dec))
        max_ra = ra + 0.6/math.cos(math.radians(dec))
        
    if (min_ra < 0):
        ra_query = "( ra > %(min_ra)f or ra < %(max_ra)f )" % {"min_ra": min_ra+360, "max_ra": max_ra,} 
    elif (max_ra > 360):
        ra_query = "( ra > %(min_ra)f or ra < %(max_ra)f )" % {"min_ra": min_ra, "max_ra": max_ra-360.,} 
    else:
        ra_query = "ra BETWEEN %(min_ra)f and %(max_ra)f" % {"min_ra": min_ra, "max_ra": max_ra,} 
        
    if (numpy.array(dec).ndim > 0):
        min_dec = dec[0]
        max_dec = dec[1]
    else:
        min_dec = dec - 0.6
        max_dec = dec + 0.6

    #
    # This query is taken from the SDSS website and selects stars with clean photometry
    # --> http://skyserver.sdss3.org/dr8/en/help/docs/realquery.asp#cleanStars
    #
    sql_query = """\
SELECT ra,dec, psfMag_u, psfMagErr_u, psfMag_g, psfMagErr_g, psfMag_r, psfMagErr_r, psfMag_i, psfMagErr_i, psfMag_z, psfMagErr_z
FROM Star 
WHERE 
%(ra_query)s AND dec BETWEEN %(min_dec)f and %(max_dec)f
AND ((flags_r & 0x10000000) != 0)
-- detected in BINNED1
AND ((flags_r & 0x8100000c00a4) = 0)
-- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
-- SATURATED, or BAD_COUNTS_ERROR
AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2))
-- not DEBLEND_NOPEAK or small PSF error
-- (substitute psfmagerr in other band as appropriate)
AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)
-- not INTERP_CENTER or not COSMIC_RAY
""" % {"filter": sdss_filter,
       "min_ra": min_ra, "max_ra": max_ra,
       "min_dec": min_dec, "max_dec": max_dec,
       "ra_query": ra_query,
       }

    if (verbose): print sql_query
    logger.debug("Downloading catalog from SDSS ...")
    logger.debug(sql_query)

    # stdout_write("Downloading catalog from SDSS ...")

    # Taken from Tomas Budavari's sqlcl script
    # see http://skyserver.sdss3.org/dr8/en/help/download/sqlcl/default.asp 
    import urllib
    # Filter out comments starting with "--"
    fsql = ""
    for line in sql_query.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep;
    params = urllib.urlencode({'cmd': fsql, 'format': 'csv'})
    # SDSS-DR 8
    # url = 'http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'

    # SDSS DR13
    url = 'http://skyserver.sdss.org/dr13/en/tools/search/sql.aspx'
    sdss = urllib.urlopen(url+'?%s' % params)
    # Budavari end


    answer = []
    for line in sdss:
        if (max_catsize > 0 and len(answer) >= max_catsize):
            break
        answer.append(line)
        if (((len(answer)-1)%10) == 0):
            if (verbose): stdout_write("\rFound %d stars so far ..." % (len(answer)-1))
    #answer = sdss.readlines()
    if (answer[0].strip() == "No objects have been found"):
        stdout_write(" nothing found\n")
        if (return_query):
            return numpy.zeros(shape=(0,12)), fsql #sql_query
        return numpy.zeros(shape=(0,12))

    # stdout_write(" found %d stars!\n" % (len(answer)-1))
    logger.debug(" found %d stars!\n" % (len(answer)-1))

    if (verbose):
        print "Returned from SDSS:"
        print "####################################################"
        print ''.join(answer)
        print "####################################################"

    # If we are here, then the query returned at least some results.
    # Dump the first line just repeating what the output was
    del answer[0]

    
    # if (verbose): print "Found %d results" % (len(answer))
    results = numpy.zeros(shape=(len(answer),12))
    # Results are comma-separated, so split them up and save as numpy array
    for i in range(len(answer)):
        items = answer[i].split(",")
        for col in range(len(items)):
            results[i, col] = float(items[col])
        #ra, dec = float(items[0]), float(items[1])
        #mag, mag_err =  float(items[2]), float(items[3])
        #results[i, :] = [ra, dec, mag, mag, mag_err, mag_err]
    
    if (return_query):
        return results, fsql #sql_query
    return results




def get_reference_catalog(ra, dec, radius, basedir, cattype="2mass_opt", verbose=False):

    logger = logging.getLogger("ReadCatalog")

    # Handle the separate case of web-based catalogs
    if (cattype.startswith("web:")):

        if (cattype == "web:2mass"):
            return twomass_from_cds(ra, dec, radius, verbose)
        elif (cattype == "web:sdss"):
            return load_catalog_from_sdss(ra, dec, 
                                          verbose=verbose, 
                                          return_query=False, 
                                          max_catsize=-1)

    # print "In get_ref_catalog, cattype=%s, dir=%s" % (cattype, basedir)

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (basedir)
    if (not os.path.isfile(skytable_filename)):
        logger.error("Unable to find catalog index file in %s!" % (basedir))
        return None

    skytable_hdu = pyfits.open(skytable_filename)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data
    #print skytable[:3]
    
    # Select entries that match our list
    # print ra, dec, radius, type(ra), type(dec), type(radius)
    #logger.debug("# Searching for stars within %.1f degress around %f , %f ..." % (radius, ra, dec))

    if (not radius == None and radius > 0):
        min_dec = dec - radius
        max_dec = dec + radius
        min_ra = ra - radius/math.cos(math.radians(dec))
        max_ra = ra + radius/math.cos(math.radians(dec))
    else:
        min_dec, max_dec = dec[0], dec[1]
        min_ra, max_ra = ra[0], ra[1]

    logger.debug("Querying catalog (%s): Ra=%f...%f Dec=%f...%f" % (
        cattype, min_ra, max_ra, min_dec, max_dec))

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
    files_to_read = [f.strip() for f in files_to_read]
    logger.debug(files_to_read)

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
            logger.debug("Reading from IPPRef catalog")
            if (verbose): print "Reading IPPRef catalog"
            catalogfile = "%s/%s.cpt" % (basedir, catalogname)

            try:
                hdu_cat = pyfits.open(catalogfile)
            except:
                logger.warning("Unable to open IPPRef file %s" % (catalogfile))
                continue

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

            array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),2+4*2))  # Ra, dec + 4*(mag/error)
            array_to_add[:,0] = cat_ra[select_from_cat]
            array_to_add[:,1] = cat_dec[select_from_cat]

            for i_mag in range(4): # this only takes the griz data, excluding the Y photometry
                # convert data from milli-mags
                # no errors, so assume all errors to be 0.000 mag
                array_to_add[:,i_mag*2+2] = photom_grizy[:,i_mag][select_from_cat] / 1e3

            if (verbose): print "# Read %d stars from catalog %s ..." % (array_to_add.shape[0], catalogfile)

        elif (cattype in ('2mass_opt', 
                          '2mass_nir', 
                          'ucac4',
                          'sdss')):
            if (verbose): print "Reading 2mass catalog"
            catalogfile = "%s/%s.fits" % (basedir, catalogname)

            try:
                hdu_cat = pyfits.open(catalogfile)
            except:
                logger.warning("Unable to open catalog (%s) file %s" % (
                    cattype, catalogfile))
                continue

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

            # sys.stderr.write("cattype=%s\n" % (cattype))
            if (cattype == '2mass_opt'):
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

            elif (cattype == '2mass_nir'):
                catalog_columns = ['RA', 'DEC',
                                   'MAG_J', 'ERR_J',
                                   'MAG_H', 'ERR_H',
                                   'MAG_K', 'ERR_K',
                               ]
                array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),len(catalog_columns)))
                for i in range(len(catalog_columns)):
                    if (catalog_columns[i] in [x.upper() for x in hdu_cat[1].data.columns.names]):
                        array_to_add[:,i] = hdu_cat[1].data.field(catalog_columns[i])[select_from_cat]
             
                # array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),5))
                # array_to_add[:,0] = cat_ra[select_from_cat]
                # array_to_add[:,1] = cat_dec[select_from_cat]
                # array_to_add[:,2] = hdu_cat[1].data.field('mag_j')[select_from_cat]
                # array_to_add[:,3] = hdu_cat[1].data.field('mag_h')[select_from_cat]
                # array_to_add[:,4] = hdu_cat[1].data.field('mag_k')[select_from_cat]

            elif (cattype == 'ucac4'):
                # sys.stderr.write("Using ucac4 catalog\n")
                catalog_columns = ['RA', 'DEC',
                                   'MAG_UCAC', 'ERR_UCAC',
                                   'MAG_B', 'ERR_B',
                                   'MAG_V', 'ERR_V',
                                   'MAG_G', 'ERR_G',
                                   'MAG_R', 'ERR_R',
                                   'MAG_I', 'ERR_I',
                               ]
                array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),len(catalog_columns)))
                for i in range(len(catalog_columns)):
                    if (catalog_columns[i] in [x.upper() for x in hdu_cat[1].data.columns.names]):
                        array_to_add[:,i] = hdu_cat[1].data.field(catalog_columns[i])[select_from_cat]
                
            elif (cattype == 'sdss'):
                catalog_columns = ['RA', 'DEC',
                                   'MAG_U', 'MAGERR_U',
                                   'MAG_G', 'MAGERR_G',
                                   'MAG_R', 'MAGERR_R',
                                   'MAG_I', 'MAGERR_I',
                                   'MAG_Z', 'MAGERR_Z',
                ]
                array_to_add = numpy.zeros(shape=(numpy.sum(select_from_cat),len(catalog_columns)))
                for col in range(len(catalog_columns)):
                    if (catalog_columns[i] in [x.upper() for x in hdu_cat[1].data.columns.names]):
                        array_to_add[:,col] = hdu_cat[1].data.field(catalog_columns[col])[select_from_cat]

            if (verbose): print "# Read %d stars from catalog %s ..." % (array_to_add.shape[0], catalogfile)

        else:
            print "This catalog name is not known"
            return None

        if (type(full_catalog) == type(None)):
            full_catalog = array_to_add
        else:
            full_catalog = numpy.append(full_catalog, array_to_add, axis=0)
        #print photom_grizy[:3,:]
        
    if (verbose): print "# Read a total of %d stars from %d catalogs!" % (full_catalog.shape[0], len(files_to_read))
    return full_catalog





if __name__ == "__main__":
    
    ra = float(get_clean_cmdline()[1])
    dec = float(get_clean_cmdline()[2])
    radius = float(get_clean_cmdline()[3])

    dec_min, dec_max = dec-radius, dec+radius
    dec_range = numpy.array([dec_min, dec_max])
    cos_dec = numpy.min( numpy.cos(numpy.radians(dec_range)))
    ra_range = numpy.array([ ra-radius/cos_dec, ra+radius/cos_dec])

    basedir = "/Volumes/odifile/Catalogs/IPPRefCat/catdir.synth.grizy/"

    import podi_sitesetup as sitesetup
    basedir = sitesetup.wcs_ref_dir
    catalog_type = sitesetup.wcs_ref_type

    basedir = cmdline_arg_set_or_default("-basedir", sitesetup.wcs_ref_dir)
    catalog_type = cmdline_arg_set_or_default("-cattype", sitesetup.wcs_ref_type)
    verbose = cmdline_arg_isset("-v")


    # catalog = get_reference_catalog(ra, dec, radius, basedir=basedir, cattype=catalog_type, verbose=verbose)
    catalog = get_reference_catalog(ra_range, dec_range, radius=None, basedir=basedir, cattype=catalog_type, verbose=verbose)

    # if (True): #False):
    #     for i in range(catalog.shape[0]):
    #         for j in range(catalog.shape[1]):
    #             print catalog[i,j],
    #         print
    if (not catalog == None):
        numpy.savetxt(sys.stdout, catalog)

