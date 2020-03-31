#!/usr/local/bin/python

import astropy.io.fits as pyfits
import math
import scipy
import sys
import os

import scipy.ndimage
import scipy.stats
import scipy.signal
import numpy
import bottleneck
from podi_definitions import *
from podi_commandline import *
import itertools

#numpy.seterr(divide='ignore', invalid='ignore')

# def rebin_image(data, binfac):
    
#     if (binfac < 1):
#         stdout_write("Rebinning at the moment only supports binning to larger pixels with binfac>1\n")
#         return None
#     elif (binfac == 1):
#         return data
    
#     out_size_x, out_size_y = int(math.ceil(data.shape[0]*1.0/binfac)), int(math.ceil(data.shape[1]*1.0/binfac))
    
#     if (out_size_x*binfac != data.shape[0] or out_size_y*binfac != data.shape[1]):
#         # The input array size is not a multiple of the new binning
#         # Create a slightly larger array to hold the data to be rebinned
#         container = numpy.zeros(shape=(out_size_x*binfac, out_size_y*binfac))
        
#         # And insert the original data
#         container[0:data.shape[0], 0:data.shape[1]] = data[:,:]
#     else:
#         container = data 
    
#     rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)
#     #rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).nm(axis=-1).nm(axis=1)
    
#     #rs = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
#     #rb1 = bottleneck.nanmean(rs, axis=-1)
#     #rb2 = bottleneck.nanmean(rb1, axis=1)
#     #rebinned = rb2
    
#     return rebinned





def find_center(hdu_data, coord_x, coord_y,
                #lx, ly, 
                prebin=8, 
                #r_minmax=[100,200], x_minmax=[400,600], y_minmax=[400,600], 
                search_x=[0, 300, 10],
                search_y=[0, 300, 10],
                search_r=[800, 1400, 10],
                fixed_radius=1270,
                # dx=2, dy=2, dr=2
                threshold=None,
                debugname=None,
                verbose=False
                ):


    #
    # Bin the data frame and the corresponding coordinate frames
    #
    ot33b = rebin_image(hdu_data, prebin, operation=numpy.mean)
    all_x = rebin_image(coord_x, prebin, operation=numpy.mean)
    all_y = rebin_image(coord_y, prebin, operation=numpy.mean)

    if (debugname != None):
        pyfits.HDUList([pyfits.PrimaryHDU(data=ot33b.T)]).writeto("debug_"+debugname+"___data.fits", overwrite=True)
        pyfits.HDUList([pyfits.PrimaryHDU(data=all_x.T)]).writeto("debug_"+debugname+"___all_x.fits", overwrite=True)
        pyfits.HDUList([pyfits.PrimaryHDU(data=all_y.T)]).writeto("debug_"+debugname+"___all_y.fits", overwrite=True)
    ot33_orig = ot33b.copy()

    ot33b[numpy.isnan(ot33b)] = 0

    #
    # Compute the gradient frame by differentiating the image 
    # frame using a Sobel filter
    #
    x33 = scipy.ndimage.sobel(ot33b, axis=0, mode='constant')
    y33 = scipy.ndimage.sobel(ot33b, axis=1, mode='constant')
    abs33 = numpy.hypot(x33, y33)
    if (debugname != None):
        pyfits.HDUList([pyfits.PrimaryHDU(data=abs33.T)]).writeto("debug_"+debugname+"___sobel.fits", overwrite=True)

    #
    # Create a mask of all valid images
    # This is needed to get rid of all the artificial edges 
    # at the edges of each cell
    # This mask is set to 1 for all NaN-masked pixels in the input frame
    # 
    mask = numpy.array(numpy.isnan(ot33_orig), dtype=numpy.float32)
    mask[3,3] = 1.
    numpy.savetxt("mask", mask)

    # 
    # Convolve the mask with a simple flat kernel to widen all gaps 
    # and get rid of cell edges
    # 
    kernel = numpy.ones(shape=(5,5))
    kernel_norm = kernel
    if (verbose): print("number valid pixels before growing",numpy.sum((mask == False)))
    mask_grown_float = scipy.ndimage.filters.convolve(mask, kernel)
    mask_grown = mask_grown_float > 0
    numpy.savetxt("mask.g", mask_grown)
    if (verbose): print("number valid pixels after growing",numpy.sum((mask_grown == False)))

    if (debugname is not None):
        pyfits.HDUList([pyfits.PrimaryHDU(data=mask_grown_float.T)]).writeto("debug_"+debugname+"___mask.fits", overwrite=True)

    #
    # Apply the now widened mask to the gradient map
    #
    abs33[numpy.isnan(ot33_orig)] = numpy.NaN
    abs33_binary = abs33.copy()
    abs33_binary[abs33 < 0.1] = 0
    abs33_binary[abs33 >= 0.1] = 1
    abs33[mask_grown] = numpy.NaN
    if (debugname != None):
        pyfits.HDUList([pyfits.PrimaryHDU(data=abs33.T)]).writeto("debug_"+debugname+"___sobel_filtered.fits", overwrite=True)
    #
    # Now apply a threshold so we only deal with strong gradients and get rid 
    # of a lot of the underlying background noise
    #
    # Figure out what contrast we need
    #
    valid_pixels = numpy.isfinite(abs33)
    #print "number valid pixels",numpy.sum(valid_pixels)
    top10percent = scipy.stats.scoreatpercentile(abs33[valid_pixels].ravel(), 90)
    if (threshold is None):
        strong_values = abs33 > top10percent
        if (verbose): print("Only using pixels >",top10percent)
    else:
        strong_values = abs33 > threshold
        if (verbose): print("Only using pixels >",threshold)


#    all_y, all_x = numpy.indices(abs33.shape)

    #
    # All we need from now on are the coordinates of pixels with strong gradients
    #
    pixel_x = all_x[strong_values]
    pixel_y = all_y[strong_values]
    pixel_value = abs33[strong_values]
    if (verbose): print(numpy.sum(strong_values),"pixels with enough signal left")

    #
    # setup the array for the pattern recognition
    #
    if (verbose):
        print("search-r=",search_r)
        print("search-x=",search_x)
        print("search-y=",search_y)
    n_x = int(math.ceil((search_x[1] - search_x[0])/search_x[2]))
    x_values_to_try = (numpy.arange(n_x) * search_x[2]) + search_x[0]
    n_y = int(math.ceil((search_y[1] - search_y[0])/search_y[2]))
    y_values_to_try = (numpy.arange(n_y) * search_y[2]) + search_y[0]
    # Make sure each ring is at least one (binned) pixel wide
    if (search_r[2] < prebin): search_r[2] = prebin
    n_r = int(math.ceil((search_r[1] - search_r[0])/search_r[2]))
    r_values_to_try = (numpy.arange(n_r+1) * search_r[2]) + search_r[0]

    #print "search box x=",x_values_to_try
    #print "search box y=",y_values_to_try
    #print "search box r=",r_values_to_try

    #
    # Now do the hough transformation:
    # Loop over all possible center positions and count how 
    # many pixels fall into the radial slices
    #
    bincount = numpy.zeros(shape=(n_x, n_y, n_r))
    bincount_w = numpy.zeros(shape=(n_x, n_y, n_r))
    for i_cx, i_cy in itertools.product(range(n_x), range(n_y)):

        cx = x_values_to_try[i_cx]
        cy = y_values_to_try[i_cy]

        pixel_radius = numpy.sqrt( (cx-pixel_x)**2 + (cy-pixel_y)**2 )

        # Count pixels in each of the rings.
        # Also use intensity as weight to emphasize stronger features
        count_w,edges = numpy.histogram(pixel_radius, bins=r_values_to_try, 
                                      weights=pixel_value)
 
        count, edges = numpy.histogram(pixel_radius, bins=r_values_to_try)
 
        # Compute the area of this ring
        outer_r = edges[1:]**2
        inner_r = edges[:-1]**2
        area_full_circle = outer_r - inner_r

        # Insert the ring count into the overall structure
        bincount[i_cx, i_cy, :] = count[:]
        bincount_w[i_cx, i_cy, :] = count_w[:] / area_full_circle

        
    if (fixed_radius is not None):
        # Find what radial bin is covered by the specified radius
        ir = int(math.floor((fixed_radius - search_r[0]) / search_r[2]))
        if (verbose): print("using fixed radius",fixed_radius, ir)

        # For this radius, find the center position with the strongest signal
        center_only = bincount_w[:,:,ir]
        index = numpy.argmax(center_only)
        #print "argmax=",index
        ix, iy = numpy.unravel_index(index, center_only.shape)
        #print "unraveled:",ix,iy
        #print x_values_to_try[ix], ix*search_x[2]+search_x[0]

    else:

        index = numpy.argmax(bincount_w)
        ix, iy, ir = numpy.unravel_index(index, bincount.shape)


    if (verbose): 
        print("found some results:")
        print("  center-x=",x_values_to_try[ix],"(",ix,")")
        print("  center-y=",y_values_to_try[iy],"(",iy,")")
        print("    radius=",r_values_to_try[ir],"(",ir,")")

    return x_values_to_try[ix], y_values_to_try[iy], r_values_to_try[ir], bincount, bincount_w, abs33




pupilghost_center_guess = {
    "OTA33.SCI": (4050, 4050),
    "OTA34.SCI": (4050, -100),
    "OTA44.SCI": (-100, -100),
    "OTA43.SCI": (-100, 4050),
    }



fiducial_centers = {
    'OTA43.SCI': (-8, 4136, 1268, -8, 4136, 1268), 
    'OTA34.SCI': (4172, -188, 1268, 4172, -188, 1268), 
    'OTA33.SCI': (4152, 4156, 1268, 4152, 4156, 1268), 
    'OTA44.SCI': (32, -190, 1268, 32, -190, 1268)
    }





def find_pupilghost_center(hdu, verbose=False, fixed_radius=1270, debugname=None):
# returns  fixed_r_x, fixed_r_y, fixed_r_r, var_r_x, var_r_y, var_r_r

    
    cx, cy = pupilghost_center_guess[hdu.header["EXTNAME"]]

    rawdata = hdu.data.T.copy()
    # trim the edges generously
    rawdata[3980:,:] = numpy.NaN
    rawdata[:,3980:] = numpy.NaN
    rawdata[:50,:] = numpy.NaN
    rawdata[:,:50] = numpy.NaN
    rawdata[500:1000,0:500] = numpy.NaN
    px, py = numpy.indices(rawdata.shape)

    search_width=250
    search_r = [1240,1320, 2]
    search_x = [cx-search_width, cx+search_width,10]
    search_y = [cy-search_width, cy+search_width,10]
    if (verbose):
        print("search-r=",search_r)
        print("search-x=",search_x)
        print("search-y=",search_y)
    prebin=8
    x, y, r, bincount, bincount_w, edge_frame = \
        find_center(rawdata, px, py,
                    search_r = search_r, search_x = search_x, search_y = search_y,
                    prebin=prebin,
                    debugname=None if (debugname==None) else (debugname + hdu.header["EXTNAME"] + "__lowres__"),
                    fixed_radius=fixed_radius,
                    verbose=verbose,
                    )
    if (verbose): print("Initial rough center (binning",prebin,") ---> ", x, y, r)

    # numpy.savetxt("bincount"+hdu[i].header["EXTNAME"], numpy.sum(numpy.sum(bincount, axis=0), axis=0))
    # log = open("bincount"+hdu[i].header["EXTNAME"]+".full", "w")
    # for x,y,z  in itertools.product(range(bincount.shape[0]), 
    #                                 range(bincount.shape[1]), 
    #                                 range(bincount.shape[2])):
    #     print >>log, x, y, z, bincount[x,y,z]
    # log.close()

    # hdu[i].data = edge_frame
    # print

    # Now we have a rough idea where the center is.
    # Extract a sub-region of the frame close to the center 
    # that contains the signal plus some extra, re-run the 
    # center-finding routine and refine the solution.
    center_x = x#*prebin
    center_y = y#*prebin
    radius = r#*prebin
    margin = 50 # pixels

    min_x, max_x = center_x - radius - margin, center_x + radius + margin
    min_y, max_y = center_y - radius - margin, center_y + radius + margin

    min_valid_x = numpy.max([0, min_x])
    max_valid_x = numpy.min([rawdata.shape[0], max_x])
    min_valid_y = numpy.max([0, min_y])
    max_valid_y = numpy.min([rawdata.shape[1], max_y])

    midres_data = rawdata[min_valid_x:max_valid_x, min_valid_y:max_valid_y]
    midres_filtered = scipy.signal.medfilt2d(midres_data, 7)
    midres_px   = px[min_valid_x:max_valid_x, min_valid_y:max_valid_y]
    midres_py   = py[min_valid_x:max_valid_x, min_valid_y:max_valid_y]

    search_r = [radius-20,radius+20,2]
    search_x = [center_x-30,center_x+30,4]
    search_y = [center_y-30,center_y+30,4]
    prebin=4
    fixed_r_x, fixed_r_y, fixed_r_r, bincount2, bincount2_w, edge_frame2 = \
        find_center(midres_filtered,
                    midres_px, midres_py,
                    search_r = search_r, search_x = search_x, search_y = search_y,
                    prebin=prebin,
                    debugname=None if (debugname==None) else (debugname + hdu.header["EXTNAME"]+"__midres_f__"),
                    fixed_radius=fixed_radius,
                    threshold=0,
                    verbose=verbose,
                    )

    if (verbose): print("Better resolution (binning",prebin,", fixed r=",fixed_radius,"px): ---> ", fixed_r_x, fixed_r_y, fixed_r_r)

    var_r_x, var_r_y, var_r_r, bincount_var, bincount_var_w, edge_frame2 = \
        find_center(midres_filtered,
                    midres_px, midres_py,
                    search_r = search_r, search_x = search_x, search_y = search_y,
                    prebin=prebin,
                    debugname=None if (debugname==None) else (debugname + hdu.header["EXTNAME"]+"__midres_v__"),
                    threshold=0,
                    verbose=verbose,
                    )

    if (verbose): print("Better resolution (binning",prebin,", variable radius):  ---> ", var_r_x, var_r_y, var_r_r)

    return fixed_r_x, fixed_r_y, fixed_r_r, var_r_x, var_r_y, var_r_r


import podi_wcs

def get_reliable_pupilghost_center(hdulist):

    prebin=8
    pg_centers = {}

    for i in range(len(hdulist)):

        if (not is_image_extension(hdulist[i])):
            continue

        print(hdulist[i].header["EXTNAME"])
        extname = hdu[i].header["EXTNAME"]

        if (not extname in fiducial_centers):
            continue

        fx, fy, fr, vx, vy, vr = find_pupilghost_center(hdu[i], verbose=True)
        print("   fixed radius:", fx, fy, fr)
        print("variable radius:", vx, vy, vr)

        pg_centers[extname] = fx, fy, fr, vx, vy, vr

        wcspoly = podi_wcs.header_to_polynomial(hdulist[i].header)
        radec  = podi_wcs.wcs_pix2wcs(numpy.array([[vx*1.0, vy*1.0]]), wcspoly, False)
        print(radec)


    print(pg_centers)

#     # angle = -140
#     pg_centers = {'OTA43.SCI': (-4, 4142, 1268, -4, 4142, 1268), 'OTA34.SCI': (4210, -170, 1268, 4210, -170, 1268), 'OTA33.SCI': (4146, 4144, 1268, 4146, 4144, 1268), 'OTA44.SCI': (16, -170, 1268, 16, -170, 1268)}

#     # angle = -5
# #    pg_centers = {'OTA43.SCI': (8, 4132, 1268, 8, 4132, 1268), 'OTA34.SCI': (4166, -178, 1268, 4166, -178, 1268), 'OTA33.SCI': (4150, 4162, 1268, 4150, 4162, 1268), 'OTA44.SCI': (22, -188, 1268, 22, -188, 1268)}

#     # fiducial angle = +5

# #    pg_centers = fiducial_centers


    # Now compute the shift between the found center 
    # and the fiducial reference solution
    all_shifts = []
    for ota, centers in pg_centers.iteritems():
        
        if (not ota in fiducial_centers):
            continue
        if (ota == 'OTA33.SCI'):
            continue

        # Get the values for this frame
        fx, fy, fr, vx, vy, vr = centers
        # Get the fiducial coordinates
        _fx, _fy, _fr, _vx, _vy, _vr = pg_centers['OTA33.SCI']

        _, _, _, fidx, fidy, _ = fiducial_centers[ota]
        _, _, _, _fidx, _fidy, _ = fiducial_centers['OTA33.SCI']

        # fiducial_centers[ota]
        
        dx = (vx - _vx) - (fidx - _fidx)
        dy = (vy - _vy) - (fidy - _fidy)
        all_shifts.append([dx, dy])

        print(ota, "this=", (vx - _vx), (vy - _vy), "    fiducial=", (fidx - _fidx), (fidy - _fidy))


    # Compute the best fitting shift
    all_shifts = numpy.array(all_shifts)
    numpy.savetxt(sys.stdout, all_shifts, "%d")

    median_shift = numpy.median(all_shifts, axis=0)
    print(median_shift)

    for ota, centers in fiducial_centers.iteritems():

        fx, fy, fr, vx, vy, vr = centers

        print(ota," --> ", vx-median_shift[0], vy-median_shift[1])



if __name__ == "__main__":


    if (cmdline_arg_isset("-template")):
        hdu = pyfits.open(get_clean_cmdline()[1])

        coord_x, coord_y = numpy.indices(hdu[1].data.shape)

        var_r_x, var_r_y, var_r_r, bincount_var, bincount_var_w, edge_frame2 = \
        find_center(hdu[1].data, coord_x, coord_y,
                    #lx, ly, 
                    prebin=8, 
                    #r_minmax=[100,200], x_minmax=[400,600], y_minmax=[400,600], 
                    search_x=[4300, 4700, 10],
                    search_y=[4300, 4700, 10],
                    search_r=[800, 1400, 10],
                    fixed_radius=None,
                    # dx=2, dy=2, dr=2
                    threshold=None,
                    debugname="___debug2____",
                    verbose=True
                    )

        print(var_r_x, var_r_y, var_r_r)




        var_r_x, var_r_y, var_r_r, bincount_var, bincount_var_w, edge_frame2 = \
        find_center(hdu[1].data, coord_x, coord_y,
                    #lx, ly, 
                    prebin=2, 
                    #r_minmax=[100,200], x_minmax=[400,600], y_minmax=[400,600], 
                    search_x=[4480, 4520, 2],
                    search_y=[4480, 4520, 2],
                    search_r=[1265, 1300, 1],
                    fixed_radius=None,
                    # dx=2, dy=2, dr=2
                    threshold=None,
                    debugname="___debug3____",
                    verbose=True
                    )

        print(var_r_x, var_r_y, var_r_r)

    elif (cmdline_arg_isset("-flat")):

        flatfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(flatfile)

        for i in range(5): #len(hdulist)):

            if (not is_image_extension(hdulist[i])):
                continue

#            print hdulist[i].header["EXTNAME"]
            extname = hdulist[i].header["EXTNAME"]
            print("\n\n\n",extname)

            fixed_r_x, fixed_r_y, fixed_r_r, var_r_x, var_r_y, var_r_r = find_pupilghost_center(
                hdulist[i], 
                verbose=True, 
                debugname="__test__",
                fixed_radius=1286)

            print("center position (fixed/var)")
            print(fixed_r_x, fixed_r_y, fixed_r_r)
            print(var_r_x, var_r_y, var_r_r)

# #     get_reliable_pupilghost_center(hdu)

#     prebin=8

#     pg_centers = {}

#     for i in range(1,5):

#         print hdu[i].header["EXTNAME"]
#         extname = hdu[i].header["EXTNAME"]
#         fx, fy, fr, vx, vy, vr = find_pupilghost_center(hdu[i], verbose=True)
#         print "   fixed radius:", fx, fy, fr
#         print "variable radius:", vx, vy, vr

#         pg_centers[extname] = fx, fy, fr, vx, vy, vr


#     print pg_centers

#     hduout = hdu[0:5]
#     hduout.writeto("/scratch/edges.fits", overwrite=True)

