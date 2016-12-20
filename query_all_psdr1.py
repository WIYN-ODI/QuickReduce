#!/usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
#  Modified 2016 by Daniel Harbeck for PS1DR1 downloaid
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

"""
Usage:

  query_all_psdr1 download [basedir]

  Start to download 1x1 degree chunks of PS1 DR1 from stsci mast via http requests.

  [param] points to the base directory in which the catalog will be generated.

  If not existsnt,a new directory basedir/raw will be created and contain the csv files with stellar photometry.


"""


import sys
import numpy
import os
import os.path
import math

from podi_definitions import *
import pyfits

import podi_photcalib
import urllib2
import concurrent.futures



def run_query_ps1dr1(minra,maxra,mindec,maxdec,  outfile, maxobj=1000000, mindet=5):
    """
    Query STSCI for PS1 DR1 data and write to file if defined.

    :param minra:
    :param maxra:
    :param mindec:
    :param maxdec:
    :param outfile:
    :param maxobj:
    :param mindet:
    :return:
    """


    # Do not duplicate efforts
    if os.path.isfile(cat_csvfile):
        print ("\t %s  already exists. skipping")
        return None

    # But tell what we are up to
    stdout_write("\nNext: %s ---> RA=%.3f...%.3f, DEC=%.3f...%.3f\n" % (
        outfile, minra, maxra, mindec, maxdec ))

    # Define the url to query catalog
    urltemplate = "http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?BBOX=%.1f,%.1f,%.1f,%.1f&FORMAT=csv&CAT=PS1V3OBJECTS&MINDET=%s&MAXOBJ=%s"


    query = urltemplate % (minra,mindec,maxra,maxdec,mindet,maxobj)
    print "\tqueryurl is: %s" % (query)

    # Get ready and query the url
    answer = []
    try:
        psdr1 = urllib2.urlopen(query, timeout=3600)

        for line in psdr1:
            answer.append(line)
            if (((len(answer)-1)%100) == 0):
                stdout_write("\r\tFound %d stars so far ..." % (len(answer)-1))

        psdr1.close()

    except Exception as e:
        print ("\tError while retrieving catlog data: %s" %(e))

    # We do not want errors. If we hae some, tell about it.
    if "Error" in answer[0]:
        print "\tError while querying catalog: %s" % (answer[0])
        return None
    else:
        stdout_write("\r\tFound a total of %d stars in PanSTARRS catalog\n" % (len(answer)-2))

    # Now we are ready to write the output file. Note that on identified query error conditions, we will not write a file.
    if (outfile != None) and (len (answer) > 1):
        try:
            thefile = open (outfile, 'w')
            for line in answer:
                thefile.write (line)
            thefile.close()
        except Exception as e:
            print ("Error opening file: %s:\n %s" % (outfile, e) )


    return answer


if __name__ == "__main__":

    #
    # Step one: download data from STSCI and stage locally.
    #

    if sys.argv[1] == 'download':
        basedir = sys.argv[2]

        rawDirectory = basedir + "/raw"
        # Increase this to make new friends at STSCI
        nsimuquery = 25


        if  not os.path.isdir(rawDirectory):
           print  "creating missing  %s " % (rawDirectory)
           os.makedirs (rawDirectory)

        pool = concurrent.futures.ThreadPoolExecutor (max_workers=nsimuquery)
        futures=[]
        increment = 0.1
        for RA in numpy.arange (0,360, increment):
            for DEC in numpy.arange (-30,90, increment):
                cat_csvfile = rawDirectory + "/" + "ps1dr%05d_%05d.csv" % (RA*10,DEC*10)

                futures.append ( pool.submit (run_query_ps1dr1, RA,RA+increment, DEC, DEC+increment, cat_csvfile) )
                #break
            break
        print ("Now let's wait for the plan to unfold")
        concurrent.futures.wait(futures)
        print ("Download done.")
        sys.exit (0)


    if sys.argv[1] == 'process':
        sys.exit (0)

    print ("No suitable action indicated on command line.")


    # startat  = int(cmdline_arg_set_or_default("-startat", 0))
#
# # Load the SkyTable so we know in what files to look for the catalog"
# skytable_filename = "%s/SkyTable.fits" % (basedir)
# skytable_hdu = pyfits.open(skytable_filename)
#
# #print skytable_hdu.info()
#
# skytable = skytable_hdu['SKY_REGION'].data
# #print skytable[:3]
#
#
# r_min = skytable.field('R_MIN')
# r_max = skytable.field('R_MAX')
#
# d_min = skytable.field('D_MIN')
# d_max = skytable.field('D_MAX')

        #
        # if (os.path.isfile(cat_csvfile)):
        #     hdu = pyfits.open(cat_fitsfile)
        #     expected = hdu[0].header['EXPCOUNT']
        #     read = hdu[0].header['CNTSRECV']
        #     print "%s: exp. %d, completed %d" % (cat_fitsfile, expected, read)
        #     hdu.close()
        #     if (read >= expected):
        #         print "This file already is complete, continuing with the next one"
        #         continue
        #     else:
        #         print "#########################"
        #         print "This file was started, but not complete, trying again:"
        #         print "@@@", cat_fitsfile
        #         print "#########################"
            
            






        # Create a catalog to hold the results:
        
        #sys.exit(0)

        # columns = [\
        #     pyfits.Column(name='RA',       format='E', unit='degrees',    array=sdss_cat[:, 0], disp='right ascension'),
        #     pyfits.Column(name='DEC',      format='E', unit='degrees',    array=sdss_cat[:, 1], disp='declination'),
        #     pyfits.Column(name='MAG_U',    format='E', unit='magnitudes', array=sdss_cat[:, 2], disp='photometry in u'),
        #     pyfits.Column(name='MAGERR_U', format='E', unit='magnitudes', array=sdss_cat[:, 3], disp='phot. error in u'),
        #     pyfits.Column(name='MAG_G',    format='E', unit='magnitudes', array=sdss_cat[:, 4], disp='photometry in g'),
        #     pyfits.Column(name='MAGERR_G', format='E', unit='magnitudes', array=sdss_cat[:, 5], disp='phot. error in g'),
        #     pyfits.Column(name='MAG_R',    format='E', unit='magnitudes', array=sdss_cat[:, 6], disp='photometry in r'),
        #     pyfits.Column(name='MAGERR_R', format='E', unit='magnitudes', array=sdss_cat[:, 7], disp='phot. error in r'),
        #     pyfits.Column(name='MAG_I',    format='E', unit='magnitudes', array=sdss_cat[:, 8], disp='photometry in i'),
        #     pyfits.Column(name='MAGERR_I', format='E', unit='magnitudes', array=sdss_cat[:, 9], disp='phot. error in i'),
        #     pyfits.Column(name='MAG_Z',    format='E', unit='magnitudes', array=sdss_cat[:,10], disp='photometry in z'),
        #     pyfits.Column(name='MAGERR_Z', format='E', unit='magnitudes', array=sdss_cat[:,11], disp='phot. error in z'),
        #     ]
        #
        # stdout_write("writing catalog ...")
        # primhdu = pyfits.PrimaryHDU()
        # #print sqlquery.replace("\n", " ")
        # primhdu.header.add_history("SQL-Query sent to SDSS:")
        # primhdu.header.add_history(sqlquery.replace("\n", " "))
        # primhdu.header["EXPCOUNT"] = (count, "expected return count")
        # primhdu.header["CNTSRECV"] = (sdss_cat.shape[0], "returned count")
        # primhdu.header["R_MIN"] = (r_min[i], "min RA")
        # primhdu.header["R_MAX"] = (r_max[i], "max RA")
        # primhdu.header["D_MIN"] = (d_min[i], "min DEC")
        # primhdu.header["D_MAX"] = (d_max[i], "max DEC")
        #
        # coldefs = pyfits.ColDefs(columns)
        # tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
        #
        # hdulist = pyfits.HDUList([primhdu, tbhdu])
        # hdulist.writeto(cat_fitsfile, clobber=True)
        #
        # stdout_write(" done!\n")
