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
import time

sys.path.append(os.path.dirname(__file__)+"/../")
from podi_definitions import *
import podi_logging
import logging
import pyfits

import podi_photcalib
import urllib2
# import concurrent.futures

import multiprocessing
import itertools
import subprocess
import time

def parallel_handler(queue, maxobj=1000000, mindet=5, id=-1):

    while (True):
        cmd = queue.get()
        if (cmd is None):
            return

        (minra,maxra,mindec,maxdec,  outfile)  = cmd
        # print cmd

        run_query_ps1dr1(minra,maxra,mindec,maxdec,  outfile, maxobj=maxobj,
                         mindet=mindet, id=id)


def wait_for_free_space(dirname, min_free, sleep=5):

    first_time = True
    logger = logging.getLogger("DiskSpace")

    if (not os.path.isdir(dirname)):
        dirname, _ = os.path.split(dirname)

    # print dirname
    while (True):
        p = subprocess.Popen(["df", dirname], stdout=subprocess.PIPE)
        output = p.communicate()[0]
        _, _, _, available, _, _ = output.split("\n")[1].split()
        free_space = int(available)

        if (free_space > min_free):
            break

        if (first_time):
            logger.warning("Waiting for more free diskspace (has %d, needs %d kilobytes)" % (
                free_space, min_free)
            )
            first_time = False

        time.sleep(sleep)



def run_query_ps1dr1(minra,maxra,mindec,maxdec,  outfile, maxobj=1000000,
                     mindet=5, id=None):
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

    worker_id = ""
    if (id is not None):
        worker_id = "Worker %3d: " % (id)

    # Do not duplicate efforts
    if os.path.isfile(outfile):
        # print ("\t %s  already exists. skipping" % (outfile))
        return None

    # if (numpy.random.random() < 0.9):
    #     return None

    # But tell what we are up to
    # stdout_write("\nNext: %s ---> RA=%.3f...%.3f, DEC=%.3f...%.3f\n" % (
    #     outfile, minra, maxra, mindec, maxdec ))

    # Define the url to query catalog
    urltemplate = "http://gsss.stsci.edu/webservices/vo/CatalogSearch.aspx?BBOX=%.1f,%.1f,%.1f,%.1f&FORMAT=csv&CAT=PS1V3OBJECTS&MINDET=%s&MAXOBJ=%s"


    query = urltemplate % (minra,mindec,maxra,maxdec,mindet,maxobj)
    #print "\tqueryurl is: %s" % (query)

    # Get ready and query the url
    answer = []
    try:
        psdr1 = urllib2.urlopen(query, timeout=3600)

        for line in psdr1:
            answer.append(line)
            #if (((len(answer)-1)%100) == 0):
            #    stdout_write("\r\tFound %d stars so far ..." % (len(
            # answer)-1))

        psdr1.close()

    except Exception as e:
        print ("\tError while retrieving catlog data: %s" %(e))

    # We do not want errors. If we hae some, tell about it.
    if "Error" in answer[0]:
        print("\tError while querying catalog: %s" % (answer[0]))
        return None
    else:
        print("%sRA=%5.1f ... %5.1f, DEC= %+5.1f ... %+5.1f ==> %5d stars" % (
            worker_id, minra, maxra, mindec, maxdec, len(answer)-2))

    # Now we are ready to write the output file. Note that on identified query error conditions, we will not write a file.
    if (outfile is not None) and (len (answer) > 0):
        try:

            # check for available disk-space (require >100MB)
            # using are kilobytes
            wait_for_free_space(outfile, 1e5)

            with  open (outfile, 'w') as thefile:
                thefile.write(os.linesep.join(answer))
            # thefile =
            # for line in answer:
            #     thefile.write (line)
            # thefile.close()
        except Exception as e:
            print ("Error opening file: %s:\n %s" % (outfile, e) )


    return answer


def convert_ascii_to_fits(raw_dir, out_dir, indexfile):

    indexhdu = pyfits.open(indexfile)
    table = indexhdu['SKY_REGION'].data

    # read data
    data = {}
    for fieldname in ['NAME',
                      'R_MIN', 'R_MAX',
                      'D_MIN', 'D_MAX']:
        data[fieldname] = table.field(fieldname)

    print data['NAME'].shape

    boxsize = 0.1
    error_log = open("%s/conversion.errors" % (out_dir), "a")
    for field in range(data['NAME'].shape[0]):

        print
        name = data['NAME'][field].strip()
        dir_name,_ = os.path.split(os.path.join(out_dir,name))
        print dir_name
        if (not os.path.isdir(dir_name)):
            os.mkdir(dir_name)
        cat_fitsfile = "%s.fits" % (os.path.join(out_dir,name))
        print cat_fitsfile

        if (os.path.isfile(cat_fitsfile)):
            print "FITS catalog %s already exists, skipping" % (cat_fitsfile)
            continue

        r_min = data['R_MIN'][field]
        r_max = data['R_MAX'][field]
        d_min = data['D_MIN'][field]
        d_max = data['D_MAX'][field]

        d_lower = math.floor(d_min / boxsize) * boxsize
        d_upper = math.ceil(d_max / boxsize) * boxsize
        r_lower = math.floor(data['R_MIN'][field] / boxsize) * boxsize
        r_upper = math.ceil(data['R_MAX'][field] / boxsize) * boxsize
        print r_min, r_max, d_min, d_max, "   ==>  ", r_lower, r_upper, d_lower, d_upper


        n_ra = (r_upper - r_lower) / boxsize
        n_dec = (d_upper - d_lower) / boxsize

        list_ra = numpy.arange(int(n_ra), dtype=numpy.float) * boxsize + r_lower
        list_dec = numpy.arange(int(n_ra), dtype=numpy.float) * boxsize + d_lower
        print list_ra[0], list_ra[-1], list_dec[0], list_dec[-1]

        collected_catalog = None
        error_abort = False
        for ra, dec in itertools.product(list_ra, list_dec):
            fn = "%s/ps1dr%05d_%05d.csv" % (raw_dir, ra*10, dec*10)
            # print ra, dec, fn, os.path.isfile(fn)

            # load catalog
            # only extract the relevant columns to speed things up
            # columns to extract: 0-indexed
            # 4,5,6,7,8: nDetections in grizy
            # 11,12: ra/dec mean
            # 13,14 ra/dec error
            # 25,26: g-psf mag/err
            # 31,32: r
            # 37,38: i
            # 43,44: z
            # 49,50: y
            # 79: r mean kron magnitude - to isolate stars
            if (not os.path.isfile(fn)):
                print "missing file %s for %s" % (fn, cat_fitsfile)
                print >> error_log, "missing file %s for %s" % (
                fn, cat_fitsfile)
                error_abort=True
                break
            try:
                skiprows=None
                with open(fn, "r") as cf:
                    lines = cf.readlines()
                    for i,l in enumerate(lines):
                        first_item = l.split(" ")[0]
                        if (first_item == "PSO"):
                            skiprows = i
                            break
                if (skiprows is None):
                    print("No data found in %s for %s" % (fn, cat_fitsfile))
                    print >>error_log, "No data found in %s for %s" % (fn, cat_fitsfile)
                    continue

                filedata = numpy.loadtxt(
                    fn, delimiter=',', skiprows=skiprows,
                    usecols=(4,5,6,7,8,  # 4,5,6,7,8: nDetections in grizy --> 0,1,2,3,4
                             11,12,      # 11,12: ra/dec mean              --> 5,6
                             13,14,      # 13,14 ra/dec error              --> 7,8
                             25,26,      # 25,26: g-psf mag/err            --> 9,10
                             31,32,      # 31,32: r                        --> 11,12
                             37,38,      # 37,38: i                        --> 13,14
                             43,44,      # 43,44: z                        --> 15,16
                             49,50,      # 49,50: y                        --> 17,18
                             79),        # 79: r mean kron magnitude - to isolate stars --> 19
                )
            except (IndexError, ValueError, StopIteration, TypeError) as e:
                print("Unable to load %s for %s: %s" % (fn, cat_fitsfile, type(e)))
                print >>error_log, "Unable to load %s for %s: %s" % (fn, cat_fitsfile, type(e))
                continue


            #
            # Now select only stars with enough detections
            #
            try:
                n_min_detections = 5
                max_phot_error = 0.2
                is_star = (filedata[:,11] - filedata[:,19]) < 0.05
                in_field = (filedata[:,5] >= r_min) & (filedata[:,5] < r_max) & \
                           (filedata[:,6] >= d_min) & (filedata[:,6] < d_max)
                good_detection = (filedata[:, 0] > n_min_detections) & \
                                 (filedata[:, 1] > n_min_detections) & \
                                 (filedata[:, 2] > n_min_detections) & \
                                 (filedata[:, 3] > n_min_detections) & \
                                 (filedata[:, 4] > n_min_detections)
                small_errors = (filedata[:, 10] > 0.) & (filedata[:, 10] < max_phot_error) & \
                               (filedata[:, 12] > 0.) & (filedata[:, 12] < max_phot_error) & \
                               (filedata[:, 14] > 0.) & (filedata[:, 14] < max_phot_error) & \
                               (filedata[:, 16] > 0.) & (filedata[:, 16] < max_phot_error) & \
                               (filedata[:, 18] > 0.) & (filedata[:, 18] < max_phot_error)
            except IndexError:
                print("Error down-selecting from %s (%s)" % (fn, cat_fitsfile))
                print >>error_log, "Error down-selecting from %s (%s)" % (fn, cat_fitsfile)
                continue

            good = is_star & in_field & good_detection & small_errors

            # trim off the n_detection fields
            cat_data = filedata[:, 5:][good]

            if (collected_catalog is None):
                collected_catalog = cat_data
            else:
                collected_catalog = numpy.append(
                    collected_catalog, cat_data, axis=0)

        if (error_abort):
            continue

        #
        # Now we have the entire catalog, prepare to write the
        # catalog as FITS file
        #
        catalog_fields = [
            'RA', 'DEC',
            'RA_err', 'DEC_err',
            'g', 'g_err',
            'r', 'r_err',
            'i', 'i_err',
            'z', 'z_err',
            'y', 'y_err',
        ]
        columns = [None] * len(catalog_fields)
        for i, fn in enumerate(catalog_fields):
            columns[i] = pyfits.Column(
                name=fn, format='D', array=collected_catalog[:,i]
            )

        primhdu = pyfits.PrimaryHDU()
        # print sqlquery.replace("\n", " ")
        primhdu.header["R_MIN"] = (r_min, "min RA")
        primhdu.header["R_MAX"] = (r_max, "max RA")
        primhdu.header["D_MIN"] = (d_min, "min DEC")
        primhdu.header["D_MAX"] = (d_max, "max DEC")

        coldefs = pyfits.ColDefs(columns)
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)

        hdulist = pyfits.HDUList([primhdu, tbhdu])
        hdulist.writeto(cat_fitsfile, clobber=True)

        print(" done!\n")

        # if (field > 10):
        #    break



if __name__ == "__main__":

    #
    # Step one: download data from STSCI and stage locally.
    #
    options = set_default_options()
    podi_logging.setup_logging(options)

    if sys.argv[1] == 'download':
        basedir = sys.argv[2]

        rawDirectory = basedir + "/raw"
        # Increase this to make new friends at STSCI
        nsimuquery = 75

        try:
            thres = float(sys.argv[3])
        except:
            thres = 0.9

        start_RA = 0
        end_RA = 360
        try:
            start_RA = float(sys.argv[4])
            end_RA = float(sys.argv[5])
        except:
            pass

        if  not os.path.isdir(rawDirectory):
           print ("creating missing  %s " % (rawDirectory))
           os.makedirs (rawDirectory)

        queue = multiprocessing.JoinableQueue()
        fields_to_download = 0
        increment = 0.1
        print("Filling up work queue")
        for RA in numpy.arange (start_RA, end_RA, increment):
            sys.stdout.write("\rRA: %6.1f" % (RA))
            sys.stdout.flush()
            for DEC in numpy.arange (-30,90, increment):

                cat_csvfile = rawDirectory + "/" + "ps1dr%05d_%05d.csv" % (RA*10,DEC*10)

                # if (numpy.random.random() < thres):
                #     continue
                if (os.path.isfile(cat_csvfile)):
                    continue

                print RA, DEC, cat_csvfile

                queue.put( (RA, RA+increment, DEC, DEC+increment, cat_csvfile))
                fields_to_download += 1
                
                #futures.append ( pool.submit (run_query_ps1dr1, RA,
                # RA+increment, DEC, DEC+increment, cat_csvfile) )

                #break
        print("Done filling up work queue!")

        print("%d fields left to be queried!" % (fields_to_download))
        time.sleep(2)

        print("Starting workers")
        processes = []
        for i in range(nsimuquery):
            kwargs = dict(
                queue=queue,
                maxobj=1000000,
                mindet=5,
                id=i+1
            )
            p = multiprocessing.Process(
                target=parallel_handler,
                kwargs=kwargs,
            )
            p.start()
            processes.append(p)
            queue.put(None)


        for p in processes:
            p.join()

        #pool = concurrent.futures.ThreadPoolExecutor (max_workers=nsimuquery)
        #futures=[]
            #break
        #print ("Now let's wait for the plan to unfold")
        #concurrent.futures.wait(futures)
        print ("Download done.")

    elif sys.argv[1] == 'process':

        raw_dir = sys.argv[2]

        fits_dir = sys.argv[3]
        index_filename = "%s/SkyTable.fits" % (fits_dir)

        convert_ascii_to_fits(raw_dir, fits_dir, index_filename)

    else:
        print ("No suitable action indicated on command line.")

    podi_logging.shutdown_logging(options)


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
