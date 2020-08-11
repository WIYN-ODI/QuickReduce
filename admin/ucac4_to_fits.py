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


# col byte item   fmt unit       explanation                            notes
# ---------------------------------------------------------------------------
#  1  1- 3 ra     I*4 mas        right ascension at  epoch J2000.0 (ICRS) (1)
#  2  5- 8 spd    I*4 mas        south pole distance epoch J2000.0 (ICRS) (1)
#  3  9-10 magm   I*2 millimag   UCAC fit model magnitude                 (2)
#  4 11-12 maga   I*2 millimag   UCAC aperture  magnitude                 (2)
#  5 13    sigmag I*1 1/100 mag  error of UCAC magnitude                  (3)
#  6 14    objt   I*1            object type                              (4)
#  7 15    cdf    I*1            combined double star flag                (5)
#          15 bytes
#  8 16    sigra  I*1 mas        s.e. at central epoch in RA (*cos Dec)   (6)
#  9 17    sigdc  I*1 mas        s.e. at central epoch in Dec             (6)
# 10 18    na1    I*1            total # of CCD images of this star
# 11 19    nu1    I*1            # of CCD images used for this star       (7)
# 12 20    cu1    I*1            # catalogs (epochs) used for proper motions
#           5 bytes
# 13 21-22 cepra  I*2 0.01 yr    central epoch for mean RA, minus 1900
# 14 23-24 cepdc  I*2 0.01 yr    central epoch for mean Dec,minus 1900
# 15 25-26 pmrac  I*2 0.1 mas/yr proper motion in RA*cos(Dec)             (8)
# 16 27-28 pmdc   I*2 0.1 mas/yr proper motion in Dec
# 17 29    sigpmr I*1 0.1 mas/yr s.e. of pmRA * cos Dec                   (9)
# 18 30    sigpmd I*1 0.1 mas/yr s.e. of pmDec                            (9)
#          10 bytes
# 19 31-34 pts_key I*4           2MASS unique star identifier            (10)
# 20 35-36 j_m    I*2 millimag   2MASS J  magnitude
# 21 37-38 h_m    I*2 millimag   2MASS H  magnitude
# 22 39-40 k_m    I*2 millimag   2MASS K_s magnitude
# 23 41    icqflg I*1            2MASS cc_flg*10 + ph_qual flag for J    (11)
# 24 42     (2)   I*1            2MASS cc_flg*10 + ph_qual flag for H    (11)
# 25 43     (3)   I*1            2MASS cc_flg*10 + ph_qual flag for K_s  (11)
# 26 44    e2mpho I*1 1/100 mag  error 2MASS J   magnitude               (12)
# 27 45     (2)   I*1 1/100 mag  error 2MASS H   magnitude               (12)
# 28 46     (3)   I*1 1/100 mag  error 2MASS K_s magnitude               (12)
#          16 bytes
# 29 47-48 apasm  I*2 millimag   B magnitude from APASS                  (13)
# 30 49-50  (2)   I*2 millimag   V magnitude from APASS                  (13)
# 31 51-52  (3)   I*2 millimag   g magnitude from APASS                  (13)
# 32 53-54  (4)   I*2 millimag   r magnitude from APASS                  (13)
# 33 55-56  (5)   I*2 millimag   i magnitude from APASS                  (13)
# 34 57    apase  I*1 1/100 mag  error of B magnitude from APASS         (14)
# 35 58     (2)   I*1 1/100 mag  error of V magnitude from APASS         (14)
# 36 59     (3)   I*1 1/100 mag  error of g magnitude from APASS         (14)
# 37 60     (4)   I*1 1/100 mag  error of r magnitude from APASS         (14)
# 38 61     (5)   I*1 1/100 mag  error of i magnitude from APASS         (14)
# 39 62    gcflg  I*1            Yale SPM g-flag*10  c-flag              (15)
#          16 bytes
# 40 63-66 icf(1) I*4            FK6-Hipparcos-Tycho source flag         (16)
# 41       icf(2) ..             AC2000       catalog match flag         (17)
# 42       icf(3) ..             AGK2 Bonn    catalog match flag         (17)
# 43       icf(4) ..             AKG2 Hamburg catalog match flag         (17)
# 44       icf(5) ..             Zone Astrog. catalog match flag         (17)
# 45       icf(6) ..             Black Birch  catalog match flag         (17)
# 46       icf(7) ..             Lick Astrog. catalog match flag         (17)
# 47       icf(8) ..             NPM  Lick    catalog match flag         (17)
# 48       icf(9) ..             SPM  YSJ1    catalog match flag         (17)
#           4 bytes
# 49 67    leda   I*1            LEDA galaxy match flag                  (18)
# 50 68    x2m    I*1            2MASS extend.source flag                (19)
# 51 69-72 rnm    I*4            unique star identification number       (20)
# 52 73-74 zn2    I*2            zone number of UCAC2 (0 = no match)     (21)
# 53 75-78 rn2    I*4            running record number along UCAC2 zone  (21)
#          12 bytes
# ---------------------------------------------------------------------------
#          78 = total number of bytes per star record



def import_ucac(catalog_dir, ucac_ascii):

    # Load the SkyTable so we know in what files to look for the catalog"
    # print "Loading lookup table"
    idx = "%s/SkyTable.fits" % (catalog_dir)
    skytable_hdu = pyfits.open(idx)

    #print skytable_hdu.info()

    skytable = skytable_hdu['SKY_REGION'].data



    # Load the UCAC file
    print("\nWorking on", ucac_ascii)

    #ucac = numpy.loadtxt(ucac_ascii)
    x = open(ucac_ascii)
    lines = x.readlines()
    all_data = []
    for line in lines:
        linevalues = [0.] * 53
        items = line.split()
        for i in range(53):
            try:
                linevalues[i] = float(items[i])
            except:
                pass
        all_data.append(linevalues)
    ucac = numpy.array(all_data)

    # Do some coordinate conversion from milli-arcsec and milli-mag to mags
    # Ra/dec are given in mas
    for i in [1,2]:
        ucac[:,i-1] /= 3600000.

    # Fix declination (given is distance from south pole)
    ucac[:,1] -= 90.

    # magnitudes are given in millimag
    for i in [3,4,      # UCAC fit/aperture magnitude
              20,21,22, # 2MASS JHK
              29,30,31,32,33 # APASS B,V,g,r,i
              ]:
        ucac[:,i-1] /= 1000.

    # magnitude errors are given in 1/100 mag
    for i in [5,
              26,27,28,
              34,35,36,37,38
          ]:
        ucac[:,i-1] *= 0.01

    # convert the central epochs for proper motion
    for i in [13,14]: ucac[:,i-1] = ucac[:,i-1] * 0.01 + 1900
    # print ucac.shape

    # Find minimum and maximum declination - each catalog has all stars at all RAs
    dec_min = numpy.min(ucac[:,1])
    dec_max = numpy.max(ucac[:,1])
    print("%d entries - declination %.2f -- %.2f" % (ucac.shape[0], dec_min, dec_max))

    # Select all catalog files in the declination range
    
    matched_declination_range = (skytable.field('D_MAX') > dec_min) & \
                                (skytable.field('D_MIN') < dec_max)

    # print dec_min, dec_max
    dmin = skytable.field('D_MIN')[matched_declination_range]
    dmax = skytable.field('D_MAX')[matched_declination_range]
    rmin = skytable.field('R_MIN')[matched_declination_range]
    rmax = skytable.field('R_MAX')[matched_declination_range]
    catname = skytable.field('NAME')[matched_declination_range]
    # print catname

    #
    # Now loop over each output catalog
    # extract the part of the UCAC catalog for this catalog file
    # and add it to the catalog. If there is no catalog yet, create it
    # 
    total = 0
    n_new = 0
    n_added = 0
    for cur_cat in range(len(catname)):
        
        full_catalog_filename = "%s/%s.fits" % (catalog_dir, catname[cur_cat])
        # print full_catalog_filename
        _dir, _file = os.path.split(full_catalog_filename)
        if (not os.path.isdir(_dir)):
            os.mkdir(_dir)

        in_this_catalog_region = (ucac[:,0] >= rmin[cur_cat]) & \
                                 (ucac[:,0]  < rmax[cur_cat]) & \
                                 (ucac[:,1] >= dmin[cur_cat]) & \
                                 (ucac[:,1]  < dmax[cur_cat])
        ucac_portion = ucac[in_this_catalog_region]
        # print ucac_portion.shape
        total += ucac_portion.shape[0]

        # Create the data
        columns = [
            pyfits.Column(name='ra',      format='E', array=ucac_portion[:, 0]),
                   pyfits.Column(name='dec',     format='E', array=ucac_portion[:, 1]),
                   pyfits.Column(name='ra_err',  format='E', array=ucac_portion[:, 7]),
                   pyfits.Column(name='dec_err', format='E', array=ucac_portion[:, 8]),

                   pyfits.Column(name='mag_ucac', format='E', array=ucac_portion[:, 3]),
                   pyfits.Column(name='err_ucac', format='E', array=ucac_portion[:, 4]),

                   pyfits.Column(name='mag_B',   format='E', array=ucac_portion[:,28]),
                   pyfits.Column(name='err_B',   format='E', array=ucac_portion[:,33]),
                   pyfits.Column(name='mag_V',   format='E', array=ucac_portion[:,29]),
                   pyfits.Column(name='err_V',   format='E', array=ucac_portion[:,34]),
                   pyfits.Column(name='mag_g',   format='E', array=ucac_portion[:,30]),
                   pyfits.Column(name='err_g',   format='E', array=ucac_portion[:,35]),
                   pyfits.Column(name='mag_r',   format='E', array=ucac_portion[:,31]),
                   pyfits.Column(name='err_r',   format='E', array=ucac_portion[:,36]),
                   pyfits.Column(name='mag_i',   format='E', array=ucac_portion[:,32]),
                   pyfits.Column(name='err_i',   format='E', array=ucac_portion[:,37]),

                   pyfits.Column(name='mag_j',   format='E', array=ucac_portion[:,19]),
                   pyfits.Column(name='err_j',   format='E', array=ucac_portion[:,26]),
                   pyfits.Column(name='mag_h',   format='E', array=ucac_portion[:,20]),
                   pyfits.Column(name='err_h',   format='E', array=ucac_portion[:,27]),
                   pyfits.Column(name='mag_k',   format='E', array=ucac_portion[:,21]),
                   pyfits.Column(name='err_k',   format='E', array=ucac_portion[:,28]),

                   pyfits.Column(name='cep_ra',    format='E', array=ucac_portion[:,12]),
                   pyfits.Column(name='cep_dec',   format='E', array=ucac_portion[:,13]),
                   pyfits.Column(name='pm_ra',     format='E', array=ucac_portion[:,14]),
                   pyfits.Column(name='pm_dec',    format='E', array=ucac_portion[:,15]),
                   pyfits.Column(name='pmerr_ra',  format='E', array=ucac_portion[:,16]),
                   pyfits.Column(name='pmerr_dec', format='E', array=ucac_portion[:,17]),

               ]

        # print columns
        coldefs = pyfits.ColDefs(columns)
        tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

        # Check if this catalog file already exist
        if (os.path.isfile(full_catalog_filename)):
            # File exist, append the data

            hdulist = pyfits.open(full_catalog_filename, mode='update')
            nrows = hdulist[1].header['NAXIS2']
            
            # Create a bigger table that holds both the existing and the new table
            nrows_total = nrows + ucac_portion.shape[0]
            tbhdu_both = pyfits.new_table(tbhdu.columns, 
                                          tbtype='BinTableHDU',
                                          nrows=nrows_total)
            # print "data.shape=",hdulist[1].data.shape
            # print "both.shape",tbhdu_both.data.shape

            # Now merge both tables
            for col in range(len(columns)):
                tbhdu_both.data.field(col)[:nrows] = hdulist[1].data.field(col)[:]
                tbhdu_both.data.field(col)[nrows:] = tbhdu.data.field(col)

            # print "Adding %s columns to file %s, new total %d" % (
            #     ucac_portion.shape[0], full_catalog_filename, nrows_total)
            n_added += 1

            hdulist[1] = tbhdu_both
            hdulist.flush()
            hdulist.close()
            pass
        else:
            # no file yet, create a new one
            # print "Creating catalog",full_catalog_filename
            primhdu = pyfits.PrimaryHDU()
            hdulist = pyfits.HDUList([primhdu, tbhdu])
            hdulist.writeto(full_catalog_filename, overwrite=True)
            n_new += 1
            pass

    # print ucac.shape, total
    print("% 3d new files, added to % 3d files." % (n_new, n_added))

    return

    print("\n\nloading data file", fits)
    fits_hdu = pyfits.open(fits)

    decs = fits_hdu[1].data.field('dec')
    ras  = fits_hdu[1].data.field('ra')
    min_dec = numpy.min(decs)
    max_dec = numpy.max(decs)

    print()
    var = min_dec, max_dec

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

    catalog_directory = sys.argv[1]

    for ucac_file in sys.argv[2:]:
        #print gzip
        # print ucac_file

        import_ucac(catalog_directory, ucac_file)
