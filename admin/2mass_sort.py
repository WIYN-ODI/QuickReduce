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
import pyfits
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

        mag_jhk[:] = numpy.NaN
        err_jhk[:] = numpy.NaN

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
                err_jhk[i] = numpy.NaN
            elif (rdflag[i] == 0):
                # Undetected
                mag_jhk[i] = numpy.NaN
                err_jhk[i] = numpy.NaN
            elif (rdflag[i] in ('6', '9')):
                # Upper limit
                err_jhk[i] = -99.9



        prox = float(items[23])


        try:
            mag_b = float(items[53])
        except:
            mag_b = numpy.NaN
            pass

        try:
            mag_vr = float(items[54])
        except:
            mag_vr = numpy.NaN
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
    columns = [\
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

    hdulist.writeto(fitsfile, clobber=True)



def sort_fits_table(idx, fits):

    # Load the SkyTable so we know in what files to look for the catalog"
    print("Loading lookup table")
    skytable_hdu = pyfits.open(idx)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data
    #print skytable[:3]

    print("\n\nloading data file", fits)
    fits_hdu = pyfits.open(fits)

    decs = fits_hdu[1].data.field('dec')
    ras  = fits_hdu[1].data.field('ra')
    min_dec = numpy.min(decs)
    max_dec = numpy.max(decs)

    print(min_dec, max_dec)

    # Now select all lookup entries in the valid range
    
    needed_catalogs = (skytable['D_MAX'] > min_dec) & (skytable['D_MIN'] < max_dec)
    
    matched_files = skytable[needed_catalogs]
    #print matched_files

    for i in range(len(matched_files)):
        #print matched_files[i]
        min_ra, max_ra = matched_files[i].field('R_MIN'), matched_files[i].field('R_MAX')
        min_dec, max_dec = matched_files[i].field('D_MIN'), matched_files[i].field('D_MAX')
    
        outputfits = matched_files[i].field('NAME')+".fits"
        #print outputfits

        # Select all the data that goes into this file
        data_for_this_file = (fits_hdu[1].data.field('dec') >= min_dec) \
                           & (fits_hdu[1].data.field('dec')  < max_dec) \
                           & (fits_hdu[1].data.field('ra')  >= min_ra) \
                           & (fits_hdu[1].data.field('ra')   < max_ra) 

        copy_data = fits_hdu[1].data[data_for_this_file]

        #print copy_data.shape

        if (not os.path.isfile(outputfits)):
            print("Preparing to copy", len(copy_data), "entries to new file", outputfits)
            # This is a new file
            primhdu = pyfits.PrimaryHDU()
            coldefs = pyfits.ColDefs(fits_hdu[1].columns)
            tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')
            tbhdu.data = copy_data
            table_hdu_out = pyfits.HDUList([primhdu, tbhdu])
            table_hdu_out.writeto(outputfits)
        else:
            print("Preparing to copy", len(copy_data), "entries to existing file", outputfits)
            # This file already exists
            table_hdu_out = pyfits.open(outputfits, mode="update")
            nrows1 = table_hdu_out[1].data.shape[0]
            nrows2 = copy_data.shape[0]
            nrows = nrows1 + nrows2
            #print "nrows=",nrows

            tbhdu = pyfits.new_table(fits_hdu[1].columns, nrows=nrows)
            for colname in range(len(fits_hdu[1].columns.names)):
                #p1 = table_hdu_out[1].data.field(colname)
                #p2 = copy_data.field(colname)
                #print p1.shape, p2.shape
                #p12 = numpy.append(p1, p2)
                #print p12.shape

                tbhdu.data.field(colname)[:nrows1] = table_hdu_out[1].data.field(colname)[:]
                tbhdu.data.field(colname)[nrows1:] = copy_data.field(colname)[:]

            table_hdu_out[1] = tbhdu
            table_hdu_out.flush()

        table_hdu_out.close()
        del table_hdu_out
            


if __name__ == "__main__":

    index = sys.argv[1]

    for fits in sys.argv[2:]:
        #print gzip

        sort_fits_table(index, fits)
