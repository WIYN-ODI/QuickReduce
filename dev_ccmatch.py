#!/usr/bin/env python


import sys
import numpy
import os
import pyfits
import datetime
import scipy
import scipy.stats
import math
import scipy.spatial
import itertools

from  podi_definitions import *
import podi_search_ipprefcat

def count_matches(src_cat, ref_cat, matching_radius=(1./60.), fine_radius=(2./3600.)):

    # 
    # Now loop over all stars in the source catalog and find nearby stars in the reference catalog
    # Use a rather large matching radius for this step
    #
    # matching_radius = 1./60. # 1 arcmin
    ref_tree = scipy.spatial.cKDTree(ref_cat)

    src_tree = scipy.spatial.cKDTree(src_cat)

    print src_cat.shape

    # find all matches
    matches = src_tree.query_ball_tree(ref_tree, matching_radius, p=2)
    # also count how many matches in total we have found
    n_matches = src_tree.count_neighbors(ref_tree, matching_radius, p=2)

    all_offsets = numpy.zeros(shape=(n_matches,2))
    cur_pair = 0
    for cur_src in range(len(matches)):
        if (len(matches[cur_src]) <= 0):
            continue

        if (verbose): print "\n",cur_src
        # print matches[cur_src]
        #
        # matches[cur_src] contains the indices of matching stars from the reference catalog
        # So extract the actual coordinates of all nearby reference stars
        #
        cur_matches = numpy.array(ref_cat[matches[cur_src]])
        if (verbose): print cur_matches

        # print cur_matches.shape

#        matches[cur_src] -= src_cat[cur_src]

        # Subtract the source position to get relative offsets
        cur_matches -= src_cat[cur_src]
        if (verbose): print cur_matches

        # And add all offsets into the global offset registry
        for cur_refstar in range(len(matches[cur_src])):
            all_offsets[cur_pair,:] = cur_matches[cur_refstar]

        cur_pair += 1

    all_offsets = all_offsets[:cur_pair,:]
    print "found", all_offsets.shape, "potential offsets", cur_pair, n_matches

    numpy.savetxt("ccmatch.dump",all_offsets)

    # 
    # At this stage, we have a catalog of all potential offsets, 
    # so we now need to figure out which one is the most likely,
    # i.e. the most frequently occuring
    #

    candidate_offset_tree = scipy.spatial.cKDTree(all_offsets)

    print all_offsets.shape

    n_coincidences = candidate_offset_tree.count_neighbors(candidate_offset_tree, fine_radius, p=2)
    coincidences = candidate_offset_tree.query_ball_tree(candidate_offset_tree, fine_radius, p=2)

    search_weights = numpy.zeros(shape=(len(coincidences),3))

    for i in range(len(coincidences)):
        search_weights[i,0] = len(coincidences[i])
        search_weights[i,1:3] = all_offsets[i,:]

    # search_weights = numpy.zeros(shape=(all_offsets.shape[0],3))
    # search_weights[:,1:3] = all_offsets
    # for i in range(all_offsets.shape[0]):
    #     neighbors = candidate_offset_tree.query_ball_point(search_weights[i,1:3], 1e-8, p=2)
    #     search_weights[i,0] = len(neighbors)
    #     #search_weights[i,0] = candidate_offset_tree.count_neighbors(
    #     #    scipy.spatial.cKDTree(search_weights[i,1:3]), (5./3600.), p=2)

    # Find which offset has the highest weight
    max_coincidence_count = numpy.argmax(search_weights[:,0])

    best_offset = all_offsets[max_coincidence_count,:]
    print "best offset", best_offset
    print "matching in",search_weights[max_coincidence_count],"fields"

    numpy.savetxt("ccmatch.offsetcount.%d" % int(round((angle*60),0)), search_weights)

    return search_weights[max_coincidence_count,0], best_offset


def rotate_shift_catalog(src_cat, center, angle, shift=None):

    center_ra, center_dec = center
    # print center_ra, center_dec

    src_rotated = numpy.zeros_like(src_cat)
    src_rel_to_center = src_cat - [center_ra, center_dec]

    # print "\n\n\nDuring rotation"
    # print "src-cat=\n",src_cat[:5,:]
    # print "src-cat rel to center=\n",src_rel_to_center[:5,:]
    # print "rotation end...\n\n"

    # angles are given in arcmin
    angle_rad = math.radians(angle)
    # print "in rot_shift: angle-rad=",angle_rad

    src_rotated[:,0] \
        = math.cos(angle_rad) * src_rel_to_center[:,0] \
        + math.sin(angle_rad) * src_rel_to_center[:,1] \
        + center_ra
    src_rotated[:,1] \
        = -math.sin(angle_rad) * src_rel_to_center[:,0] \
        +  math.cos(angle_rad) * src_rel_to_center[:,1] \
        + center_dec
    
    if (not shift == None):
        src_rotated += shift

    return src_rotated


if __name__ == "__main__":
    verbose=False

    # Load the source catalog file
    src_catfile = sys.argv[1]
    src_cat = numpy.loadtxt(src_catfile)
    # eliminate all flagged stars
    src_cat = src_cat[src_cat[:,7] == 0][:,0:2]

    print "src_cat:",src_cat.shape

    # 
    # Create the reference catalog
    #
    ref_catfile = sys.argv[2]
    ref_raw = numpy.loadtxt(ref_catfile)[:,0:2]
    print "ref. cat (raw) =",ref_raw.shape


    #
    # Reduce the reference catalog to approx. the coverage of the source catalog
    #
    ref_cat = match_catalog_areas(src_cat, ref_raw, 2./60.)
    print "area matched ref. catalog:", ref_cat.shape

    #
    # compute the center of the field
    #
    center_ra = numpy.median(src_cat[:,0])
    center_dec = numpy.median(src_cat[:,1])
    print "field center at ", center_ra, center_dec

    #
    # For testing purposes, rotate the field by a little
    #
    testing = True
    if (testing):
        angle = 10./60. # 10 arcmin
        src_cat = rotate_shift_catalog(src_cat, (center_ra, center_dec), angle)

    # 
    # Now loop over all stars in the source catalog and find nearby stars in the reference catalog
    # Use a rather large matching radius for this step
    #
    matching_radius = 1./60. # 1 arcmin
    ref_tree = scipy.spatial.cKDTree(ref_cat)

    angle_max = 2.
    d_angle = 5.
    n_angles = int(math.ceil((2 * angle_max) / (d_angle / 60.))) + 1

    all_results = numpy.zeros(shape=(n_angles, 4))

    all_results[:,0] = numpy.linspace(-angle_max, angle_max, n_angles)
    for cur_angle in range(n_angles):

        angle = all_results[cur_angle,0]
        print "\n\n\n",angle*60

        src_rotated = rotate_shift_catalog(src_cat, (center_ra, center_dec), angle, None)
        print "src-cat:",src_cat.shape
        print "src_roated:",src_rotated.shape
        print "ref-cat:",ref_cat.shape

        # angle_rad = math.radians(angle)

        # src_rotated = numpy.zeros_like(src_cat)
        # src_rel_to_center = src_cat - [center_ra, center_dec]
        # src_rotated[:,0] \
        #     = math.cos(angle_rad) * src_rel_to_center[:,0] \
        #     + math.sin(angle_rad) * src_rel_to_center[:,1] \
        #     + center_ra
        # src_rotated[:,1] \
        #     = -math.sin(angle_rad) * src_rel_to_center[:,0] \
        #     +  math.cos(angle_rad) * src_rel_to_center[:,1] \
        #     + center_dec

        n_matches, offset = count_matches(src_rotated, ref_cat, matching_radius)

        all_results[cur_angle,1] = n_matches
        all_results[cur_angle,2:4] = offset

    numpy.savetxt(sys.stdout, all_results)
    numpy.savetxt("ccmatch.allresults", all_results)


    #
    # Now find the best solution (the one with the highest matched star density)
    #
    idx_best_angle = numpy.argmax(all_results[:,1])

    best_guess = all_results[idx_best_angle]
    print best_guess, "angle=",best_guess[0]*60.,"arcmin"

    # 
    # With this best guess at hand, match each star in the source 
    # catalog to a closest match in the reference catalog
    # 

    print best_guess
    src_rotated = rotate_shift_catalog(src_cat, (center_ra, center_dec), 
                                       angle=best_guess[0], 
                                       shift=-best_guess[2:4])

    numpy.savetxt("ccmatch.roughalign", src_rotated)
