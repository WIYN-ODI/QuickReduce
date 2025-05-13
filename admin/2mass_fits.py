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
import os
import astropy.io.fits as pyfits
import numpy
import math



def ascii2fits(asciifile, fitsfile):

    file = open(asciifile, "rb")
    f = file.readlines()


    list_of_entries = []

    linecounter = 0
    for line in f:

        linecounter += 1
        sys.stdout.write("\rReading line %d ..." % (linecounter))
        sys.stdout.flush()

        # Split line into items - separator is '|'

        items = line.split('|')

        ra = float(items[0])
        dec = float(items[1])

        #print items[0:10]

        err_maj, err_min, err_ang = float(items[2]), float(items[3]), float(items[4])
        err_ra =  math.sqrt((err_maj * math.sin(math.radians(err_ang)))**2 + 
                            (err_min * math.cos(math.radians(err_ang)))**2)
        err_dec = math.sqrt((err_maj * math.cos(math.radians(err_ang)))**2 + 
                            (err_min * math.sin(math.radians(err_ang)))**2)

        mag_jhk = numpy.zeros(shape=(3))
        err_jhk = numpy.zeros(shape=(3))

        mag_jhk[:] = numpy.nan
        err_jhk[:] = numpy.nan

        try:
            mag_jhk[0], err_jhk[0] = float(items[6]), float(items[8])
        except:
            pass

        try:
            mag_jhk[1], err_jhk[1] = float(items[10]), float(items[12])
        except:
            pass

        try:
            mag_jhk[2], err_jhk[2] = float(items[14]), float(items[16])
        except:
            pass

        #Here: Add some logic to take care of photflat and rdflag
        phqual, rdflag, blflag, ccflag = items[18], items[19], items[20], items[21]
        for i in range(3):
            if (numpy.isnan(mag_jhk[i])):
                continue

            if (ccflag[i] != '0'):
                # Photometry affected by persistence, spikes, etc.
                err_jhk[i] = numpy.nan
            elif (rdflag[i] == 0):
                # Undetected
                mag_jhk[i] = numpy.nan
                err_jhk[i] = numpy.nan
            elif (rdflag[i] in ('6', '9')):
                # Upper limit
                err_jhk[i] = -99.9



        prox = float(items[23])


        try:
            mag_b = float(items[53])
        except:
            mag_b = numpy.nan
            pass

        try:
            mag_vr = float(items[54])
        except:
            mag_vr = numpy.nan
            pass


        this_star = [ra, dec, err_ra, err_dec,  
                     mag_jhk[0], err_jhk[0], mag_jhk[1], err_jhk[1], mag_jhk[2], err_jhk[2],
                     prox, mag_b, mag_vr
                     ]

        #print items[0], mag_jhk, phqual, mag_b, mag_vr
        #print this_star

        list_of_entries.append(this_star)

    full_data = numpy.array(list_of_entries)

    #Create fits table
    columns = [
        pyfits.Column(name='ra',      format='E', array=full_data[:,0]),
        pyfits.Column(name='dec',     format='E', array=full_data[:,1]),
        pyfits.Column(name='ra_err',  format='E', array=full_data[:,2]),
        pyfits.Column(name='dec_err', format='E', array=full_data[:,3]),
        pyfits.Column(name='mag_j',   format='E', array=full_data[:,4]),
        pyfits.Column(name='err_j',   format='E', array=full_data[:,5]),
        pyfits.Column(name='mag_h',   format='E', array=full_data[:,6]),
        pyfits.Column(name='err_h',   format='E', array=full_data[:,7]),
        pyfits.Column(name='mag_k',   format='E', array=full_data[:,8]),
        pyfits.Column(name='err_k',   format='E', array=full_data[:,9]),
        pyfits.Column(name='prox',    format='E', array=full_data[:,10]),
        pyfits.Column(name='mag_b',   format='E', array=full_data[:,11]),
        pyfits.Column(name='mag_vr',  format='E', array=full_data[:,12]),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    #coldefs = pyfits.ColDefs([c1, c2, c3, c4, c5])
    #tbhdu = pyfits.new_table(coldefs)



    #tbhdu.data = full_data

    #binhdu = pyfits.BinTableHDU(header=tbhdu.header, data=tbhdu.data)

    primhdu = pyfits.PrimaryHDU()

    hdulist = pyfits.HDUList([primhdu, tbhdu])

    hdulist.writeto(fitsfile, overwrite=True)





if __name__ == "__main__":

    for gzip in sys.argv[1:]:
        #print gzip

        asciifile = "/scratch/temp.ascii"
        fitsfile = gzip[:-3] + ".fits"

        cmd = "gunzip -c %s > %s" % (gzip, asciifile)
        os.system(cmd)

        print(cmd, fitsfile)
        ascii2fits(asciifile, fitsfile)
