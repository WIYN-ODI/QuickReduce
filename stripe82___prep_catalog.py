#!/usr/bin/env python

#
# Download the necessary catalog from:
# http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82.html
# http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82calibStars_v2.6.dat.gz
#


### This is file stripe82calibStars_v2.6.dat created by ZI on Mar 8,2007
###  
### A catalog of 1,006,849 candidate standard stars from SDSS stripe 82
### 
### Details are described in Ivezic et al. 2007 (AJ, astro-ph/)
### A quick selection summary is:
###  1) unresolved in imaging, at least one band with photometric
###     error below 0.05 mag
###  2) processing flags BRIGHT, SATUR, BLENDED, or EDGE are not set
###  3) at least 4 observations in gri
###  4) non-variable (chi2 < 3 in gri)
###  5) the final standard error of the mean r band mag: <0.05 mag
###
### There is one line per star, and each line lists (in this order):
###  0) every line starts with the string CALIBSTARS
###  1) RA Dec RArms Decrms: the mean position and its rms per coordinate,
###     this is J2000, decimal degrees for RA and Dec, and arcsec for rms
###     NB: standard errors can be computed as rms/sqrt(Ntot)
###  2) Ntot: the total number of epochs
###  3) Ar: the Schlegel, Finkbeiner & Davis (1998) ISM extinction value in 
###     the r band; extinction in other bands can be computed as [Rv=3.1]: 
###     Am = Cm*Ar, with Cm=(1.873, 1.377, 0.758, 0.537) for m=(ugiz) 
###  4) and then in each band (ugriz, in this order):
###       (Nobs mmed mmu msig mrms mchi2), which are: 
###       the total number of observations in this band
###       the median magnitude 
###       the mean magnitude
###       the standard error for the mean (1.25 larger for the median)
###       the root-mean-square scatter
###       chi2 per degree of freedom (computed using the mean magnitude)
### 
### **** IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT **** 
###  1) To select sources with reliable photometry in the u and z bands
###     don't forget to require Nobs >= 4
###  2) to avoid a slight bias (~0.02 mag) at the faint end in the gri  
###     bands, require msig*sqrt(Nobs) < 0.03 
### **** IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT **** 
### 
###  For more details, in case of problems, etc., send email to Z. Ivezic
###                                         (ivezic@astro.washington.edu)
### 

#
# Columns
#

#  1 (  ): CALIBSTARS
#  2 ( 0): Ra (deg)
#  3 ( 1): Dec (deg)
#  4 ( 2): Ra_rms (arcsec)
#  5 ( 3): Dec_rms (arcsec)
#  6 ( 4): Ntot
#  7 ( 5): A_R
#
#  8 ( 6): u N_Obs 
#  9 ( 7): u mag_median 
# 10 ( 8): u mag_mean
# 11 ( 9): u std.err mean
# 12 (10): u rms
# 13 (11): u chi^2/dof
#
# 14 (12): g N_Obs 
# 14 (13): g mag_median 
# 16 (14): g mag_mean
# 17 (15): g std.err mean
# 18 (16): g rms
# 19 (17): g chi^2/dof
#
# 20 (18): r N_Obs 
# 21 (19): r mag_median 
# 22 (20): r mag_mean
# 23 (21): r std.err mean
# 24 (22): r rms
# 25 (23): r chi^2/dof
#
# 26 (24): i N_Obs 
# 27 (25): i mag_median 
# 28 (26): i mag_mean
# 29 (27): i std.err mean
# 30 (28): i rms
# 31 (29): i chi^2/dof
#
# 32 (30): z N_Obs 
# 33 (31): z mag_median 
# 34 (32): z mag_mean
# 35 (33): z std.err mean
# 36 (34): z rms
# 37 (35): z chi^2/dof


# using awk, select only stars with Nobs>=4 in u and z
# awk '{if ($8>=4 && $32>=4) print $0; }' stripe82calibStars_v2.6.dat > stripe82_ugriz.select


# in Python:
# Load catalog into numpy, split it into a number of sub-catalogs, each for one filter and save as numpy binary array

import os
import numpy

raw_file = "stripe82calibStars_v2.6.dat"
compatible_file = "stripe82calibStars_v2.6.4np"

if (not os.path.isfile(compatible_file)):
    cmd = "awk '{if (substr($0,0,1)!=\"#\") print substr($0,11,10000);}' %s > %s" % (raw_file, compatible_file)
    print cmd
    os.system(cmd)

print "Loading the catalog file"
full_cat = numpy.loadtxt(compatible_file)
print "Done loading, found",full_cat.shape[0],"objects"

# Now take apart catalog

filters = ['u', 'g', 'r', 'i', 'z']

for i in range(len(filters)):

    #
    # Apply the filters listed above as "IMPORTANT IMPORTANT"
    #
    col_NObs = i*6 + 6
    if (filters[i] in ('u', 'z')):
        select = full_cat[:,col_NObs] >= 4
    else:
        col_msig = 6*i + 9
        select = full_cat[:,col_msig]*numpy.sqrt(full_cat[:,col_NObs]) < 0.03

    #  Create the output array   
    number_good_stars = numpy.sum(select)
    catalog = numpy.zeros(shape=(number_good_stars, 6))

    # Then go ahead and fill in the values
    catalog[:,0] = full_cat[:,0][select]
    catalog[:,1] = full_cat[:,1][select]

    # Figure out in what columns the data is
    mag_median = 6*i +  7
    mag_mean   = 6*i +  8
    mag_std    = 6*i +  9
    mag_rms    = 6*i + 10

    # Finish up the catalog for this filter
    catalog[:,2] = full_cat[:,mag_median][select]
    catalog[:,3] = full_cat[:,mag_mean  ][select]
    catalog[:,4] = full_cat[:,mag_std   ][select]
    catalog[:,5] = full_cat[:,mag_rms   ][select]
    
    # Write it to file, as a numpy binary array
    catalog_filename = "stripe82__%s.cat" % filters[i]
    print "Writing catalog to file (%s)" % (catalog_filename)
    numpy.save(catalog_filename, catalog)
