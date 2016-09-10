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

import podi_photcalib
import podi_search_ipprefcat

def run_query(sql_query):

    # Taken from Tomas Budavari's sqlcl script
    # see http://skyserver.sdss3.org/dr8/en/help/download/sqlcl/default.asp 
    import urllib
    # Filter out comments starting with "--"
    fsql = ""
    for line in sql_query.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep;
    #print "_______",fsql,"_________"
    params = urllib.urlencode({'cmd': fsql, 'format': 'csv'})
    #print params

    # SDSS-DR 8
    # url = 'http://skyserver.sdss3.org/dr8/en/tools/search/x_sql.asp'

    # SDSS DR13
    url = 'http://skyserver.sdss.org/dr13/en/tools/search/sql.aspx'

    sdss = urllib.urlopen(url+'?%s' % params)
    # Budavari end

    answer = []
    for line in sdss:
        answer.append(line)
        if (((len(answer)-1)%10) == 0):
            stdout_write("\rFound %d stars so far ..." % (len(answer)-1))
    #answer = sdss.readlines()
    if (answer[0].strip() == "No objects have been found"):
        return numpy.zeros(shape=(0,0))
    else:
        stdout_write("\rFound a total of %d stars in SDSS catalog\n" % (len(answer)-1))

    return answer


if __name__ == "__main__":
    
    basedir = sys.argv[1]
    startat  = int(cmdline_arg_set_or_default("-startat", 0))

    # Load the SkyTable so we know in what files to look for the catalog"
    skytable_filename = "%s/SkyTable.fits" % (basedir)
    skytable_hdu = pyfits.open(skytable_filename)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data
    #print skytable[:3]
    

    r_min = skytable.field('R_MIN')
    r_max = skytable.field('R_MAX')
    
    d_min = skytable.field('D_MIN')
    d_max = skytable.field('D_MAX')

    outputfile = skytable.field("name")

    print r_min
    print outputfile[0:5]

    for i in range(startat, len(outputfile)):
        directory, name = os.path.split(outputfile[i])
        stdout_write("\nNext: %s/%s ---> RA=%.3f-%.3f, DEC=%.3f...%.3f\n" % (
                directory, name, r_min[i], r_max[i], d_min[i], d_max[i], ))

        if (not os.path.isdir(directory)):
            os.mkdir(directory)
        
        cat_fitsfile = outputfile[i]+".fits"
        if (os.path.isfile(cat_fitsfile)):
            hdu = pyfits.open(cat_fitsfile)
            expected = hdu[0].header['EXPCOUNT']
            read = hdu[0].header['CNTSRECV']
            print "%s: exp. %d, completed %d" % (cat_fitsfile, expected, read)
            hdu.close()
            if (read >= expected):
                print "This file already is complete, continuing with the next one"
                continue
            else:
                print "#########################"
                print "This file was started, but not complete, trying again:"
                print "@@@", cat_fitsfile
                print "#########################"
            
            
        count_entries = """\
SELECT COUNT(*) from Star
where
ra>%(ramin)f and ra<%(ramax)f and
dec between %(min_dec)f and %(max_dec)f
AND ((flags_r & 0x10000000) != 0)
-- detected in BINNED1
AND ((flags_r & 0x8100000c00a4) = 0)
-- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
-- SATURATED, or BAD_COUNTS_ERROR
AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2))
-- not DEBLEND_NOPEAK or small PSF error
-- (substitute psfmagerr in other band as appropriate)
AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)
-- not INTERP_CENTER or not COSMIC_RAY""" % \
            {"ramin": r_min[i], 
             "ramax": r_max[i],
#             "min_dec": -30, 
#             "max_dec": -32,}
             "min_dec": d_min[i], 
             "max_dec": d_max[i],}

        #print count_entries

        answer = run_query(count_entries)
        #print "___",answer,"___"

        if (len(answer) > 2):
            print ''.join(answer)
            continue

        count = int(answer[1].strip())
        print "---> expecting %s stars as result" % (count)
        #print ''.join(answer)

        if (count > 0):
            
            catalog, sqlquery = podi_search_ipprefcat.load_catalog_from_sdss(
                [r_min[i], r_max[i]], [d_min[i], d_max[i]], "r", return_query=True)
            )
            # podi_photcalib.load_catalog_from_sdss([r_min[i], r_max[i]], [d_min[i], d_max[i]], "r", return_query=True)
            #print catalog.shape

            #print catalog[:10]
            sdss_cat = catalog
        else:
            sqlquery = "No results expected, did not query"
            sdss_cat = numpy.zeros(shape=(0,12))
            print "nothing expected..."

        # Create a catalog to hold the results:
        
        #sys.exit(0)

        columns = [\
            pyfits.Column(name='RA',       format='E', unit='degrees',    array=sdss_cat[:, 0], disp='right ascension'),
            pyfits.Column(name='DEC',      format='E', unit='degrees',    array=sdss_cat[:, 1], disp='declination'),
            pyfits.Column(name='MAG_U',    format='E', unit='magnitudes', array=sdss_cat[:, 2], disp='photometry in u'),
            pyfits.Column(name='MAGERR_U', format='E', unit='magnitudes', array=sdss_cat[:, 3], disp='phot. error in u'),
            pyfits.Column(name='MAG_G',    format='E', unit='magnitudes', array=sdss_cat[:, 4], disp='photometry in g'),
            pyfits.Column(name='MAGERR_G', format='E', unit='magnitudes', array=sdss_cat[:, 5], disp='phot. error in g'),
            pyfits.Column(name='MAG_R',    format='E', unit='magnitudes', array=sdss_cat[:, 6], disp='photometry in r'),
            pyfits.Column(name='MAGERR_R', format='E', unit='magnitudes', array=sdss_cat[:, 7], disp='phot. error in r'),
            pyfits.Column(name='MAG_I',    format='E', unit='magnitudes', array=sdss_cat[:, 8], disp='photometry in i'),
            pyfits.Column(name='MAGERR_I', format='E', unit='magnitudes', array=sdss_cat[:, 9], disp='phot. error in i'),
            pyfits.Column(name='MAG_Z',    format='E', unit='magnitudes', array=sdss_cat[:,10], disp='photometry in z'),
            pyfits.Column(name='MAGERR_Z', format='E', unit='magnitudes', array=sdss_cat[:,11], disp='phot. error in z'),
            ]

        stdout_write("writing catalog ...")
        primhdu = pyfits.PrimaryHDU()
        #print sqlquery.replace("\n", " ")
        primhdu.header.add_history("SQL-Query sent to SDSS:")
        primhdu.header.add_history(sqlquery.replace("\n", " "))
        primhdu.header["EXPCOUNT"] = (count, "expected return count")
        primhdu.header["CNTSRECV"] = (sdss_cat.shape[0], "returned count")
        primhdu.header["R_MIN"] = (r_min[i], "min RA")
        primhdu.header["R_MAX"] = (r_max[i], "max RA")
        primhdu.header["D_MIN"] = (d_min[i], "min DEC")
        primhdu.header["D_MAX"] = (d_max[i], "max DEC")

        coldefs = pyfits.ColDefs(columns)
        tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

        hdulist = pyfits.HDUList([primhdu, tbhdu])
        hdulist.writeto(cat_fitsfile, clobber=True)

        stdout_write(" done!\n")
