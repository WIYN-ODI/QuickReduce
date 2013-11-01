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

    # print src_cat.shape

    # find all matches
    matches = src_tree.query_ball_tree(ref_tree, matching_radius, p=2)
    # also count how many matches in total we have found
    n_matches = src_tree.count_neighbors(ref_tree, matching_radius, p=2)

    all_offsets = numpy.zeros(shape=(n_matches,2))
    cur_pair = 0
    for cur_src in range(len(matches)):
        if (len(matches[cur_src]) <= 0):
            continue

        #if (verbose): print "\n",cur_src
        # print matches[cur_src]
        #
        # matches[cur_src] contains the indices of matching stars from the reference catalog
        # So extract the actual coordinates of all nearby reference stars
        #
        cur_matches = numpy.array(ref_cat[matches[cur_src]])
        #if (verbose): print cur_matches

        # print cur_matches.shape

#        matches[cur_src] -= src_cat[cur_src]

        # Subtract the source position to get relative offsets
        cur_matches -= src_cat[cur_src]
        #if (verbose): print cur_matches

        # And add all offsets into the global offset registry
        for cur_refstar in range(len(matches[cur_src])):
            all_offsets[cur_pair,:] = cur_matches[cur_refstar]

        cur_pair += 1

    all_offsets = all_offsets[:cur_pair,:]
    # print "found", all_offsets.shape, "potential offsets", cur_pair, n_matches

    numpy.savetxt("ccmatch.dump",all_offsets)

    # 
    # At this stage, we have a catalog of all potential offsets, 
    # so we now need to figure out which one is the most likely,
    # i.e. the most frequently occuring
    #

    candidate_offset_tree = scipy.spatial.cKDTree(all_offsets)

    # print all_offsets.shape

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


def rotate_shift_catalog(src_cat, center, angle, shift=None, verbose = False):
    
    if (verbose):
        print "\n\n\nIn rotate_shift_catalog"
        print "angle =", angle
        print "shift =", shift
        print "center =", center
        print "src-cat=\n", src_cat[:3]

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
    if (verbose): print "angle radians =",angle_rad

    # print "in rot_shift: angle-rad=",angle_rad

    src_rotated[:,0] \
        = math.cos(angle_rad) * src_rel_to_center[:,0] \
        + math.sin(angle_rad) * src_rel_to_center[:,1] \
        + center_ra
    src_rotated[:,1] \
        = -math.sin(angle_rad) * src_rel_to_center[:,0] \
        +  math.cos(angle_rad) * src_rel_to_center[:,1] \
        + center_dec
    
    if (verbose): print "src_rotated=\n", src_rotated[:3]

    if (not shift == None):
        print "applying shift", shift
        src_rotated += shift

    if (verbose): print "src-final=\n", src_rotated[:3],"\n\n\n"
        
    return src_rotated

def kd_match_catalogs(src_cat, ref_cat, matching_radius, max_count=1):

    src_tree = scipy.spatial.cKDTree(src_cat[:,0:2])
    ref_tree = scipy.spatial.cKDTree(ref_cat[:,0:2])

    print src_cat[0:5]
    print ref_cat[0:5]

    # Create an array to hold the matched catalog
    output_cat = numpy.empty(shape=(src_cat.shape[0], src_cat.shape[1]+ref_cat.shape[1]))
    # and insert the source catalog
    n_src_columns = src_cat.shape[1]
    output_cat[:,0:n_src_columns] = src_cat

    print output_cat[0:5]

    # also create an array holding for which sources we found a match
    match_found = numpy.zeros(shape=(src_cat.shape[0]))

    # match the catalogs using a kD-tree
    match_indices = src_tree.query_ball_tree(ref_tree, matching_radius, p=2)

    print src_tree.count_neighbors(ref_tree, matching_radius, p=2)

    # Now loop over all matches and merge the found matches
    for cur_src in range(src_cat.shape[0]):

        # Determine how many reference stars are close to this source
        # Do not keep match if none or too many reference stars are nearby
        n_matches = len(match_indices[cur_src])
        if (n_matches <= 0 or n_matches > max_count):
            continue
    
        output_cat[cur_src, n_src_columns:] = ref_cat[match_indices[cur_src][0]]
        match_found[cur_src] = 1

    # Now eliminate all sources without matches
    final_cat = output_cat[match_found == 1]

    return final_cat



def find_best_guess(src_cat, ref_cat,
                    center_ra, center_dec,
                    matching_radius=(1./60.),
                    angle_max=2., #degrees
                    d_angle=5 # arcmin
                    ):

    # 
    # Now loop over all stars in the source catalog and find nearby stars in the reference catalog
    # Use a rather large matching radius for this step
    #
    #matching_radius = 1./60. # 1 arcmin
    ref_tree = scipy.spatial.cKDTree(ref_cat)

    #angle_max = 2.
    #d_angle = 5.
    n_angles = int(math.ceil((2 * angle_max) / (d_angle / 60.))) + 1

    all_results = numpy.zeros(shape=(n_angles, 4))

    all_results[:,0] = numpy.linspace(-angle_max, angle_max, n_angles)
    for cur_angle in range(n_angles):

        angle = all_results[cur_angle,0]
        # print "\n\n\n",angle*60

        src_rotated = rotate_shift_catalog(src_cat, (center_ra, center_dec), angle, None)
        # print "src-cat:",src_cat.shape
        # print "src_roated:",src_rotated.shape
        # print "ref-cat:",ref_cat.shape

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


    return best_guess


if __name__ == "__main__":
    verbose=False

    # Load the source catalog file
    src_catfile = sys.argv[1]
    src_raw = numpy.loadtxt(src_catfile)
    # eliminate all flagged stars
    src_cat = src_raw[src_raw[:,7] == 0][:,0:2]

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
    ref_cat = match_catalog_areas(src_cat, ref_raw, 4./60.)
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
        angle = 7./60. # 10 arcmin
        src_cat = rotate_shift_catalog(src_cat, (center_ra, center_dec), angle, [0.01, -0.005])

    #
    # Find 1st order best guess
    #
    # best_guess = find_best_guess(src_cat, ref_cat,
    #                              center_ra, center_dec,
    #                              matching_radius=(1./60.),
    #                              angle_max=2., #degrees
    #                              d_angle=5 # arcmin
    #                              )
    # print "\n\n\n\n\n found best guess:"
    # print best_guess, "\n\n\n\n\n"
    best_guess = [ -8.33333333e-02, 1.02000000e+02, -9.91732294e-03, 5.02854045e-03]

    sys.exit(0)

    # 
    # With this best guess at hand, match each star in the source 
    # catalog to a closest match in the reference catalog
    # 


    print best_guess
    src_rotated = rotate_shift_catalog(src_cat, (center_ra, center_dec), 
                                       angle=best_guess[0], 
                                       shift=best_guess[2:4])

    numpy.savetxt("ccmatch.roughalign", src_rotated)

    # Match up stars
    ref_tree = scipy.spatial.cKDTree(ref_cat)
    src_tree = scipy.spatial.cKDTree(src_rotated)
    matching_radius = 3./3600.

    matched_src_ref_idx = src_tree.query_ball_tree(ref_tree, matching_radius, p=2)

    src_ref_pairs = numpy.ones(shape=(src_rotated.shape[0],4))
    src_ref_pairs[:,0:2] = src_rotated[:,0:2]
    src_ref_pairs[:,2:4] = numpy.NaN

    for i in range(len(matched_src_ref_idx)):
        # Ignore all points with no or more than 1 closest match
        n_close_stars = len(matched_src_ref_idx[i])
        if (not n_close_stars == 1):
            continue

        src_ref_pairs[i, 2:4] = ref_cat[matched_src_ref_idx[i][0]]

    numpy.savetxt("ccmatch.srcrefmatched", src_ref_pairs)

    #
    # Further optimize the rotation angle by introducing the 
    # shift and rotation as free parameters and fitting to minimize 
    # the deviations
    #
    def difference_ref_src_catalogs(src, ref, declination):
        diff = src - ref
        diff[:,0] *= math.cos(math.radians(declination))
        distance = numpy.hypot(diff[:,0], diff[:,1])
        # return diff.ravel() #stance
        return distance

    # Optimize rotation angle by fitting offsets with rotation as a free parameter
    def optimize_rotation_angle(p, src, ref, center):
        ra, dec = center
        print p
        src_rotated = rotate_shift_catalog(src, center, 
                                          angle=p[0], 
                                          shift=p[1:2])
        return difference_ref_src_catalogs(src_rotated, ref, dec)

    r_mismatch = difference_ref_src_catalogs(src_ref_pairs[:,0:2], src_ref_pairs[:,2:4], center_dec)

    p_init = [best_guess[0], best_guess[2], best_guess[3]]
    print p_init

    x_rotated = rotate_shift_catalog(src_cat, (center_ra, center_dec), 
                                     angle=best_guess[0], shift=best_guess[2:4])

    def difference_source_reference_cat(p, src_cat, ref_cat, center, for_fitting=False):
        # (center_ra, center_dec), 
#        dummy = open("dummy.txt","w")
        src_rotated = rotate_shift_catalog(src_cat, center, 
                                           angle=p[0], 
                                           shift=p[1:3])
#        numpy.savetxt(dummy, src_cat)
#        print >>dummy, "\n\n\n\n\n"
#        numpy.savetxt(dummy, src_rotated)
#        print >>dummy, "\n\n\n\n\n"
#        numpy.savetxt(dummy, ref_cat)
#        print >>dummy, "\n\n\n\n\n"

        diff = src_rotated - ref_cat
#        numpy.savetxt(dummy, diff)
#        dummy.close()

#        print "src_cat=", src_cat.shape
#        print "src_cat_rot=", src_rotated.shape
#        print "ref_cat=", ref_cat.shape
#        print "diff_cat=", diff.shape

        if (for_fitting):
            return diff.ravel()

        return diff


    numpy.savetxt("ccmatch.x_rotated", x_rotated)

    numpy.savetxt("ccmatch.r_mismatch", r_mismatch)

    center_radec = (center_ra, center_dec)
    diff = difference_source_reference_cat(p_init, src_cat, src_ref_pairs[:,2:4], center_radec)
    numpy.savetxt("ccmatch.diff", diff)



    valid_matches = numpy.isfinite(src_ref_pairs[:,2])
    matched_src = src_cat[valid_matches]
    matched_ref = src_ref_pairs[:,2:4][valid_matches]

    diff2 = difference_source_reference_cat(p_init, matched_src, matched_ref, center_radec, for_fitting=False)
    numpy.savetxt("ccmatch.diff2", diff2)

    args = (matched_src, matched_ref, center_radec, True)
    fit = scipy.optimize.leastsq(difference_source_reference_cat, p_init, args=args, full_output=1)

    print "\n\nbefore/after fit"
    print p_init
    print fit[0]
    # Compute uncertainty on the shift and rotation
    uncert = numpy.sqrt(numpy.diag(fit[1]))
    print uncert
    best_shift_rotation_solution = fit[0]


    diff_afterfit = difference_source_reference_cat(fit[0], matched_src, matched_ref, center_radec, for_fitting=False)
    numpy.savetxt("ccmatch.diff_afterfit", diff_afterfit)

    # 
    # Now we have the best fitting rotation and shift position.
    # In a next step, go ahead and refine the distortion in the frame
    # 
    print "\n\n\nOptimizing distortion..."
    print "Using initial guess", best_shift_rotation_solution


    inputframe = sys.argv[3]
    print "optimizing rotation in frame",inputframe
    
    hdulist = pyfits.open(inputframe)
    hdulist.info()

    # Apply the best-fit rotation to all coordinates in source catalog
    full_src_cat = src_raw.copy()
    numpy.savetxt("ccmatch.opt_src", full_src_cat)
    src_rotated = rotate_shift_catalog(full_src_cat[:,0:2], (center_ra, center_dec), 
                                       angle=best_shift_rotation_solution[0], 
                                       shift=best_shift_rotation_solution[1:3],
                                       verbose=True)
    numpy.savetxt("ccmatch.opt_rot", src_rotated)
    full_src_cat[:,0:2] = src_rotated[:,0:2]
    # Now we have the full source catalog, with better matching star coordinates
    # Next, eliminate all stars with flags
    valid_stars = full_src_cat[:,7] == 0
    valid_src_cat = full_src_cat[valid_stars]

    diff = full_src_cat[:,0:2] - src_raw[:,0:2]
    diff2 = src_rotated[:,0:2] - src_raw[:,0:2]
    
    numpy.savetxt("ccmatch.opt_ref", ref_cat)
    numpy.savetxt("ccmatch.opt_valid", valid_src_cat)
    numpy.savetxt("ccmatch.opt_diff", diff)
    numpy.savetxt("ccmatch.opt_diff2", diff2)

    # Next we can do the actual coordinate matching between the source and 
    # reference star catalog
    matched_catalog = kd_match_catalogs(valid_src_cat, ref_cat, matching_radius=(2./3600))

    print matched_catalog[:,8]

    for ext in range(len(hdulist)):
        if (not is_image_extension(hdulist[ext])):
            continue

        ota_extension = hdulist[ext]
        # print ota_extension.header['FPPOS'], ota_extension.header['FPPOS'][2:4]
        ota = int(ota_extension.header['FPPOS'][2:4])
        print "working on OTA %02d ..." %(ota)

        # sources from this OTA
        in_this_ota = (matched_catalog[:,8] == ota)
        print numpy.sum(in_this_ota)
        print in_this_ota
        print matched_catalog[:,8]

        ota_cat = matched_catalog[in_this_ota]

        print "sources in ota %d = %s ..." % (ota, str(ota_cat.shape))

        

    sys.exit(0)

    args = ( matched_src, matched_ref, (center_ra, center_dec) )
    fit = scipy.optimize.leastsq(optimize_rotation_angle, p_init, args=args, full_output=1)
    
    print fit[0]

    d_r = optimize_rotation_angle(fit[0], matched_src, matched_ref, (center_ra, center_dec))
    numpy.savetxt("ccmatch.d_r_afterfit", d_r)
    d_r = optimize_rotation_angle(p_init, matched_src, matched_ref, (center_ra, center_dec))
    numpy.savetxt("ccmatch.d_r_beforefit", d_r)
