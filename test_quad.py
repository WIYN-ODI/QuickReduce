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

def select_brightest(radec, mags, n):
    # print mags

    n = mags.shape[0] if n > mags.shape[0] else n
    si = numpy.argsort(mags[:,0])
    # print si

    output_radec = numpy.zeros(shape=(n, radec.shape[1]))
    output_mags = numpy.zeros(shape=(n, mags.shape[1]))

    for i in range(n):
        # print si[i], mags[si[i],0]
        output_radec[i,:] = radec[si[i],:]
        output_mags[i,:] = mags[si[i],:]

    return output_radec, output_mags

@profile
def compute_quad_ratios(ref_coord):
    areas = numpy.zeros(shape=(4))
    area_ratios = numpy.zeros(shape=(2))

    n_points = ref_coord.shape[0]
    n_quads = n_points * (n_points-1) * (n_points-2) * (n_points-3) / 24

    ratios = numpy.zeros(shape=(n_quads,2))
    indices = numpy.zeros(shape=(n_quads,4))
    point_vector = numpy.zeros(shape=(4))
    all_areas = numpy.zeros(shape=(n_quads,4))

    sorted_areas = numpy.zeros_like(areas)
    sorted_vector = numpy.zeros_like(point_vector)

    #
    # Precompute the area of all triangles we might need later on
    #
    triangle_areas = numpy.zeros(shape=(n_points,n_points,n_points))
    # for (p1,p2,p3) in itertools.combinations(range(n_points),3):
    #     triangle_areas[p1,p2,p3]  = math.fabs((ref_coord[p1,0]-ref_coord[p2,0])*(ref_coord[p2,1]-ref_coord[p3,1]) - (ref_coord[p2,0]-ref_coord[p3,0])*(ref_coord[p1,1]-ref_coord[p2,1]))
        
    for p1 in range(ref_coord.shape[0]-2):
        for p2 in range(p1+1, ref_coord.shape[0]-1):
            for p3 in range(p2+1, ref_coord.shape[0]):
                triangle_areas[p1,p2,p3]  = math.fabs((ref_coord[p1,0]-ref_coord[p2,0])*(ref_coord[p2,1]-ref_coord[p3,1]) - (ref_coord[p2,0]-ref_coord[p3,0])*(ref_coord[p1,1]-ref_coord[p2,1]))
    
    cur_quad = 0
    brightest_ref_radec = ref_coord
    for p1 in range(brightest_ref_radec.shape[0]-3):
        for p2 in range(p1+1, brightest_ref_radec.shape[0]-2):
            for p3 in range(p2+1, brightest_ref_radec.shape[0]-1):
                for p4 in range(p3+1, brightest_ref_radec.shape[0]):

                    #
                    # Compute areas of the four possible triangles
                    #
                    # triangle opposite of p1
                    all_areas[cur_quad,0] = triangle_areas[p2,p3,p4]
                    # triangle opposite of p2
                    all_areas[cur_quad,1] = triangle_areas[p1,p3,p4]
                    # triangle opposite of p3
                    all_areas[cur_quad,2] = triangle_areas[p1,p2,p4]
                    # triangle opposite of p4
                    all_areas[cur_quad,3] = triangle_areas[p1,p2,p3]


                    # Save which points form this triangle
                    indices[cur_quad,0] = p1
                    indices[cur_quad,1] = p2
                    indices[cur_quad,2] = p3
                    indices[cur_quad,3] = p4
                    # The above is faster than the more elegant
                    # point_vector[:] = [p1,p2,p3,p4]

                    cur_quad += 1

    ## Create area ratios by dividing by the largest area
    #area_ratios[0] = sorted_areas[2] / sorted_areas[3]
    #area_ratios[1] = sorted_areas[1] / sorted_areas[3]

    print "all-areas:\n",all_areas[0:5,:]

    sortindex = numpy.argsort(all_areas, axis=1)
    print "sortindex\n",sortindex[0:5,:]

    xxx = numpy.ravel_multi_index((numpy.arange(all_areas.shape[0])[:, None], sortindex), dims=all_areas.shape)
    
    print "xxx=\n",xxx[0:5,:]

    sorted_xxx = (all_areas.ravel()[xxx]).reshape(all_areas.shape)
    all_areas = sorted_xxx

    print "test-sorted=\n", sorted_xxx[0:5,:]

    print "indices before sort\n",indices[0:5,:]
    indices = (indices.ravel()[xxx]).reshape(indices.shape)
    print "indices after sort\n",indices[0:5,:]

    ratios[:,0] = all_areas[:,2] / all_areas[:,3]
    ratios[:,1] = all_areas[:,1] / all_areas[:,3]

    # print ratios
    return ratios, indices, all_areas


@profile
def find_optimal_transformation(reference_coords, catalog_coords, verbose=False, skip_inplausible=True, max_rotation=20):

    ref_coord = reference_coords
    brightest_cat_radec  = catalog_coords

    ref_ratios, ref_indices, ref_areas = compute_quad_ratios(ref_coord)
    if (verbose):
        print "ref-ratios:\n",ref_ratios,"\nend of ref_ratios",ref_ratios.shape

    search_tree = scipy.spatial.cKDTree(ref_ratios)

    # print "\n\n\n\n\ncat before offset"
    cat_ratios, cat_indices, cat_areas = compute_quad_ratios(brightest_cat_radec)
    # print "cat-ratios:\n",cat_ratios,"\nend of cat_ratios",cat_ratios.shape

    # brightest_cat_radec *= [0.8,1.4]
    # print "\n\n\n\n\ncat after offset"
    # cat_ratios, cat_indices, cat_areas = compute_quad_ratios(brightest_cat_radec)
    # print "cat-ratios:\n",cat_ratios,"\nend of cat_ratios",cat_ratios.shape

    coord_tree = scipy.spatial.cKDTree(cat_ratios)

#   matches = search_tree.query_ball_tree(coord_tree, 0.0001, p=2)
    matches = search_tree.query_ball_tree(coord_tree, 0.001, p=2)

    # print search_tree

    # print matches

    if (verbose):
        print "# matches=",len(matches)
        print "# ref quads",ref_ratios.shape
        print "# cat quads",cat_ratios.shape
        print cat_ratios[:10,:]


    # Count how many pais we are going to find
    n_matches = search_tree.count_neighbors(coord_tree, 0.001, p=2)

    all_transformations = numpy.zeros((n_matches,6))

    cur_match = 0
    for i in range(len(matches)):
        if (len(matches[i]) <= 0):
            continue

        if (verbose):
            print 
            print matches[i]

        match_indices = matches[i]
        for j in range(len(match_indices)):

            if (verbose):
                print "input #",i,"/",len(matches),", match #",j+1," --> "
                print "       ",match_indices[j], ref_ratios[i], cat_ratios[match_indices[j]], cat_indices[match_indices[j]], ref_indices[i]
                print "       ",reference_coords[ref_indices[i,0]], brightest_cat_radec[cat_indices[match_indices[j],0]]
                print "       ",reference_coords[ref_indices[i,1]], brightest_cat_radec[cat_indices[match_indices[j],1]]
                print "       ",reference_coords[ref_indices[i,2]], brightest_cat_radec[cat_indices[match_indices[j],2]]
                print "       ",reference_coords[ref_indices[i,3]], brightest_cat_radec[cat_indices[match_indices[j],3]]

            ref_x, ref_y, cat_x, cat_y = numpy.zeros(shape=(4)), numpy.zeros(shape=(4)), numpy.zeros(shape=(4)), numpy.zeros(shape=(4))
            for p in range(4):
                ref_x[p] = reference_coords[ref_indices[i,p],0]
                cat_x[p] = brightest_cat_radec[cat_indices[match_indices[j],p],0]
                ref_y[p] = reference_coords[ref_indices[i,p],1]
                cat_y[p] = brightest_cat_radec[cat_indices[match_indices[j],p],1]

            if (verbose):
                print "       Ref:", ref_x, "    Cat:",cat_x
                print "       Ref:", ref_y, "    Cat:",cat_y

            # Now we have a set of corresponding points, so work out what coordinate transformation we need 
            # Coordinate transformation is as follows:
            # x' = a*x + b*y + c*1
            # y' = d*x * e*y * f*1

            # Create the input matrix for the abc fit
            ABC = numpy.vstack([ref_x, ref_y, numpy.ones(len(ref_x))]).T
            # print ABC
            # solve equation: 
            transformation_x, residuals_x, rank_x, s_x = numpy.linalg.lstsq(ABC, cat_x)
            # print transformation_x

            if (transformation_x[0] > 1.1 or transformation_x[0] < 0.8):
                continue

            DEF = numpy.vstack([ref_x, ref_y, numpy.ones(len(ref_x))]).T
            transformation_y, residuals_y, rank_y, s_y = numpy.linalg.lstsq(DEF, cat_y)
            # print transformation_y

            # Check if this solution is a plausible solution, i.e. 
            # mostly rotation and offset with little shear
            if (skip_inplausible):
                rel_difference = (transformation_x[0] - transformation_y[1]) / (math.fabs(transformation_x[0]) + math.fabs(transformation_y[1]))
                if (math.fabs(rel_difference) > 0.05):
                    continue
                determinant = transformation_x[0] * transformation_y[1] - transformation_x[1] * transformation_y[0]
                if (determinant > 1.1 or determinant < 0.9):
                    continue

            # now combine both x/y transformation into one
            transformation_xy = numpy.append(transformation_x, transformation_y).reshape((1,6))

            all_transformations[cur_match, :] = transformation_xy

            if (verbose):
                print "       ",transformation_xy

            # all_transformations = numpy.append(all_transformations, transformation_xy, axis=0)
            cur_match += 1


        # print cat_ratios[i,:]
        # print [ref_ratios[matches[i][j]] for j in range(len(matches[i]))]
        if (verbose):
            print

    if (verbose):
        print all_transformations

    all_transformations = all_transformations[:cur_match,:]

    #
    # Now we have all matches, go ahead and find the most frequent solution
    # potentially needed: some scaling for the offset 
    #
    all_transformations[:,2] /= 100.
    all_transformations[:,5] /= 100.
    transform_tree = scipy.spatial.cKDTree(all_transformations)

    matching_transforms = transform_tree.query_ball_tree(transform_tree, 0.003, p=2)
    # matching_transforms = transform_tree.count_neighbors(transform_tree, 0.001, p=2)

    # Now loop over all matching transformations and count how many times we find this or a very similar match:
    n_matches = numpy.zeros(shape=(len(matching_transforms)))
    for i_match in range(len(matching_transforms)):
        n_matches[i_match] = len(matching_transforms[i_match])

    # print matching_transforms
    if (verbose):
        print "Number of matches:\n",n_matches

    numpy.savetxt("test_quad.n_matches", n_matches)

    # Find the transformation with the most matches
    idx_best_transformation = numpy.argmax(n_matches)

    numpy.savetxt("test_quad.all_transform", all_transformations)

    best_transformation = all_transformations[idx_best_transformation]
    if (verbose):
        print "Found a best transformation:"
        print "    ",best_transformation

    return best_transformation, n_matches[idx_best_transformation]



@profile
def match_catalog_areas(src, to_match, radius):


    src_radec = src[:,0:2]
    src_tree = scipy.spatial.cKDTree(src_radec)

    match_radec = to_match[:,0:2]
    match_tree = scipy.spatial.cKDTree(match_radec)

    # Now select all neighbors within the search radius
    matches = match_tree.query_ball_tree(src_tree, radius, p=2)

    valid = numpy.ones(shape=to_match.shape[0])

    for i in range(len(matches)):
        if (len(matches[i]) > 0):
            # This means we found a star in the source catalog close enough to this star
            # don't do anything
            pass
        else:
            # We didn't find a nearby neighbor, so this is a source to be removed.
            valid[i] = 0

    matched = to_match[valid == 1]
    return matched

if __name__ == "__main__":

    if (cmdline_arg_isset("-random")):
        # Create a list of random coordinates
        radec = numpy.random.random((100,2))
        mags = numpy.random.random((100,1)) * 5 + 15

        # print mags

        # randomly select 50 stars from the first list
        n_select = 50
        select = numpy.random.random((n_select)) * radec.shape[0]
        # print select


        cat_radec = numpy.zeros(shape=(n_select,2))
        cat_mags = numpy.zeros(shape=(n_select,1))
        for i in range(n_select):
            ii = int(select[i])
            cat_radec[i,:] = radec[ii,:]
            cat_mags[i,:] = mags[ii,:]

        # print cat_mags


        brightest_ref_radec, brightest_ref_mags = select_brightest(radec, mags, 25)
        brightest_cat_radec, brightest_cat_mags = select_brightest(cat_radec, cat_mags, 10)
        # print "cat before offset",brightest_cat_radec


        # brightest_cat_radec *= [0.8,1.4]
        # brightest_cat_radec += [0,0]
        # print "cat after offset",brightest_cat_radec

        brightest_cat_radec += [0.2,0.1]
        scatter = numpy.random.randn(brightest_cat_radec.shape[0], brightest_cat_radec.shape[1]) * (0.3 / 3600.)
        print "scatter=",scatter.shape, brightest_cat_radec.shape
        print scatter * 3600.


        brightest_cat_radec += scatter

        # print brightest_ref_mags
        # print brightest_cat_mags


        n_stars = brightest_ref_radec.shape[0]
        n_triangles = n_stars * (n_stars-1) * (n_stars-2) / 6
        ref_triangle_cat = numpy.zeros(shape=(n_triangles, 10))

        # Create triangle catalog for reference catalog
        i_triangle = 0



        ref_coord = brightest_ref_radec

        best_transformation = find_optimal_transformation(ref_coord, brightest_cat_radec)

        print "Found a best transformation:"
        print "    ",best_transformation

        sys.exit(0)

    if (cmdline_arg_isset("-data")):
        reference_file = get_clean_cmdline()[1]
        catalog_file = get_clean_cmdline()[2]

        ref_cat = numpy.loadtxt(reference_file)
        mag_column = int(cmdline_arg_set_or_default("-refmag", 0))-1
        ref_count = int(cmdline_arg_set_or_default("-refcount", 25))

        ref_radec = ref_cat[:,0:2]
        if (mag_column > 0):
            mags = ref_cat[:,mag_column].reshape((ref_cat.shape[0],1))
            ref_radec, ref_mags = select_brightest(ref_radec, mags, ref_count)
            

        src_cat = numpy.loadtxt(catalog_file)
        mag_column = int(cmdline_arg_set_or_default("-srcmag", 0))-1
        src_count = int(cmdline_arg_set_or_default("-srccount", 10))

        src_radec = src_cat[:,0:2]
        if (mag_column > 0):
            mags = src_cat[:,mag_column].reshape((src_cat.shape[0],1))
            src_radec, src_mags = select_brightest(src_radec, mags, src_count)
            

        # print src_radec
        # print ref_radec
        print "\n\nReference catalog:"
        numpy.savetxt(sys.stdout, ref_radec)
        print "\n\nSource catalog:"
        numpy.savetxt(sys.stdout, src_radec)

        print "\n\nGoing to work"
        best_transformation = find_optimal_transformation(ref_radec, src_radec)

        print "Found a best transformation:"
        print "    ",best_transformation

        sys.exit(0)





    if (cmdline_arg_isset("-data2")):
        catalog_file = get_clean_cmdline()[1]

        src_cat = numpy.loadtxt(catalog_file)

        # Eliminate all stars with flags
        # src_cat = src_cat[src_cat[:,7] == 0]

        mag_column = int(cmdline_arg_set_or_default("-srcmag", 0))-1
        src_count = int(cmdline_arg_set_or_default("-srccount", 10))

        src_radec = src_cat[:,0:2]
        if (mag_column > 0):
            mags = src_cat[:,mag_column].reshape((src_cat.shape[0],1))
            src_radec, src_mags = select_brightest(src_radec, mags, src_count)
            

        # Load the reference catalog
        # ra = numpy.median(src_cat[:,0])
        # dec = numpy.median(src_cat[:,1])
        print src_cat[:,0]
        ra = scipy.stats.nanmedian(src_cat[:,0])
        dec = scipy.stats.nanmedian(src_cat[:,1])
        print "Center of field ~", ra, dec

        # ipp_cat = podi_search_ipprefcat.get_reference_catalog(ra, dec, 0.7, "/Volumes/odifile/Catalogs/2mass_fits/")
        ipp_cat = podi_search_ipprefcat.get_reference_catalog(ra, dec, 0.7, "/datax/2mass_fits/")
        print "Found in 2mass",ipp_cat.shape

        matched_cat = match_catalog_areas(src_radec, ipp_cat, (2./60.))
        numpy.savetxt("ref.cat.raw", ipp_cat) 
        numpy.savetxt("ref.cat.matched", matched_cat) 

        ref_count = int(cmdline_arg_set_or_default("-refcount", 40))

        ref_radec = matched_cat[:,0:2]                                                                                                                 
        ref_mag = matched_cat[:,3].reshape((matched_cat.shape[0],1))

        ref_radec = ref_radec[(ref_mag[:,0] > 13) & (ref_mag[:,0]<16)]
        ref_mag = ref_mag[(ref_mag[:,0] > 13) & (ref_mag[:,0]<16)]
        ref_radec, ref_mag = select_brightest(matched_cat, ref_mag, ref_count)
        print ref_radec.shape, ref_count

        numpy.savetxt("test_quad.ref_radec", ref_radec)
        numpy.savetxt("test_quad.src_radec", src_radec)

        # print src_radec
        # print ref_radec
        print "\n\nReference catalog:"
        numpy.savetxt(sys.stdout, ref_radec)
        print ref_radec.shape

        print "\n\nSource catalog:"
        numpy.savetxt(sys.stdout, src_radec)
        print src_radec.shape

        print "\n\nGoing to work"

        #import cProfile, pstats
        ##cProfile.run("""best_transformation = find_optimal_transformation(ref_radec, src_radec)""", "profiler")
        #cProfile.run("""best_transformation, n_matches = find_optimal_transformation(ref_radec[:,0:2], src_radec)""", "profiler")
        #p = pstats.Stats("profiler")
        #p.strip_dirs().sort_stats('time').print_stats()
        #p.sort_stats('time').print_stats()

        best_transformation, n_matches = find_optimal_transformation(ref_radec[:,0:2], src_radec)

        # best_transformation = find_optimal_transformation(ref_radec[:,0:2], src_radec)

        print "Found a best transformation:"
        print "    ",best_transformation
        print " # matches =",n_matches


        sys.exit(0)
