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
from query_usno import query_usno
from podi_definitions import *
import pyfits
#import date
import datetime
#import pywcs  
from astLib import astWCS
import pdb
import scipy
import scipy.stats


import podi_matchcatalogs

output_debug_catalogs = False

def apply_transformation(p,x):
    x_ = x.copy()
    x_[:,0] = p[0] + x[:,0] * math.cos(p[2]) + x[:,1] * math.sin(p[2])
    x_[:,1] = p[1] - x[:,0] * math.sin(p[2]) + x[:,1] * math.cos(p[2])
    return x_

def improve_match_and_rotation(fixwcs_ref_ra, fixwcs_ref_dec,
                               fixwcs_odi_ra, fixwcs_odi_dec,
                               wcs_shift,
                               matching_radius=2, n_repeats=3,
                               verbose=False):
    
    # Apply the pre-defined correction
    fixwcs_odi_ra[:] += wcs_shift[0]
    fixwcs_odi_dec[:] += wcs_shift[1]

    # Now lets search for matching catalogs
    matching_radius = 2. / 3600. # = 2'' in degrees

    ref_full = numpy.empty(shape=(fixwcs_ref_ra.shape[0],2))
    ref_full[:,0] = fixwcs_ref_ra[:]
    ref_full[:,1] = fixwcs_ref_dec[:]
    print "Read ",ref_full.shape, "reference stars"

    odi_full = numpy.empty(shape=(fixwcs_odi_ra.shape[0],4))
    odi_full[:,0] = fixwcs_odi_ra[:]
    odi_full[:,1] = fixwcs_odi_dec[:]
    odi_full[:,2] = fixwcs_odi_ra[:]
    odi_full[:,3] = fixwcs_odi_dec[:]
    print "Read ",odi_full.shape, "ODI stars"

    odi_orig = odi_full.copy()

#    def fitfunc(p,x):
#        x_ = x.copy()
#        x_[:,0] = p[0] + x[:,0] * math.cos(p[2]) + x[:,1] * math.sin(p[2])
#        x_[:,1] = p[1] - x[:,0] * math.sin(p[2]) + x[:,1] * math.cos(p[2])
#        return x_
    def errfunc(p,ref,odi):
        return (ref - apply_transformation(p,odi)).flatten()

    p_total = [0,0,0]

    for repeat in range(n_repeats):
        return_cat = podi_matchcatalogs.match_catalogs(ref_full, odi_full)
        if (verbose): print "return_Cat=",return_cat

        # Now we should have almost matching catalogs, with the exception 
        # of the mismatch in rotation

        good_matches = return_cat[:,2] > 0
        matched_cat = return_cat[good_matches]

        if (verbose): print "left are ",matched_cat.shape,"good matches"

        if (verbose): print matched_cat[0:5,:]

        p_init = [0., 0., 0.]
        out = scipy.optimize.leastsq(errfunc, p_init,
                                     args=(matched_cat[:,0:2], matched_cat[:,2:4]),
                                     full_output=1)

        #print out
        p_final = out[0]
        covar = out[1]

        p_total += p_final

        if (verbose): print "best transformation:",p_final
        angle = math.degrees(p_final[2])
        if (verbose): print "mismatched angle=",angle*60,"arcmin"

        if (output_debug_catalogs): numpy.savetxt("matched.cat.%d" % (repeat), matched_cat)

        # Now apply the new correction to the full array of stars
        odi_full[:,0:2] = apply_transformation(p_final, odi_full[:,0:2])


    # Now we have iteratively refined the matched catalog
    # In a final step use the original catalog to get the full transformation
    # A backup of the unaltered corrdinates are at the back of the odi_full 
    # array and subsequently copied into the matched_cat array
    p_init = p_total
    out = scipy.optimize.leastsq(errfunc, p_init,
                                 args=(matched_cat[:,0:2], matched_cat[:,4:6]),
                                 full_output=1)

    p_final = out[0]
    covar = out[1]

    print "final transformation:", p_final, "angle=", 60.*math.degrees(p_final[2])
    if (verbose): print "straight sum of transformations: ", p_total

    odi_corr = apply_transformation(p_total, odi_orig)

    # and re-match the catalogs to see if now we can match more stars
    return_cat = podi_matchcatalogs.match_catalogs(ref_full, odi_corr) #, [0,360], [0,360])
    matched_cat = return_cat[return_cat[:,2]>=0]

    if (output_debug_catalogs): numpy.savetxt("matched.out", matched_cat)


    print "after some fiddling we can now match",matched_cat.shape,"stars"

    return p_final


if __name__ == "__main__":

    # First load all the data
    fixwcs_ref_ra = numpy.loadtxt("numsave.fixwcs_ref_ra.txt")
    fixwcs_ref_dec  = numpy.loadtxt("numsave.fixwcs_ref_dec.txt")
    fixwcs_odi_ra = numpy.loadtxt("numsave.fixwcs_odi_ra.txt")
    fixwcs_odi_dec = numpy.loadtxt("numsave.fixwcs_odi_dec.txt")
    wcs_shift_guess = numpy.loadtxt("numsave.wcs_shift_guess.txt")
    wcs_shift_refinement = numpy.loadtxt("numsave.wcs_shift_refinement.txt")

    wcs_shift = wcs_shift_guess + wcs_shift_refinement
       
    final_transform = improve_match_and_rotation(fixwcs_ref_ra, fixwcs_ref_dec,
                               fixwcs_odi_ra, fixwcs_odi_dec,
                               wcs_shift,
                               matching_radius=2, n_repeats=5,
                               verbose=True)

    # Use the average position of the reference catalog
    med_ra = numpy.median(fixwcs_ref_ra)
    med_dec = numpy.median(fixwcs_ref_dec)
    print med_ra, med_dec

    dx = math.cos(final_transform[2])*med_ra + math.sin(final_transform[2])*med_dec
    dy = -math.sin(final_transform[2])*med_ra + math.cos(final_transform[2])*med_dec

    print final_transform
    print dx, dy

    print med_ra-dx, med_dec-dy

