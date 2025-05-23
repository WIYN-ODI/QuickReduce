#! /usr/bin/env python3
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


"""

Experimental routine to fit and subtract a 2-d surface of the sky-background of
a given input file.

"""

import sys
import os
import astropy.io.fits as pyfits
import numpy
import scipy
import scipy.stats
import scipy.interpolate
import math

from podi_definitions import *
from podi_commandline import *
from podi_plotting import *

from astLib import astWCS

import bottleneck

boxwidth = 20


def find_maximum_extent(minmax_radec, wcs, max_xy):

    ra1, dec1 = wcs.pix2wcs(1, 1)
    ra2, dec2 = wcs.pix2wcs(1, max_xy[1])
    ra3, dec3 = wcs.pix2wcs(max_xy[0], 1)
    ra4, dec4 = wcs.pix2wcs(max_xy[0], max_xy[1])

    minmax_radec[0] = numpy.min( [minmax_radec[0], ra1, ra2, ra3, ra4] )
    minmax_radec[1] = numpy.max( [minmax_radec[1], ra1, ra2, ra3, ra4] )

    minmax_radec[2] = numpy.min( [minmax_radec[2], dec1, dec2, dec3, dec4] )
    minmax_radec[3] = numpy.max( [minmax_radec[3], dec1, dec2, dec3, dec4] )

    return minmax_radec

import itertools

def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = numpy.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = numpy.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(numpy.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = numpy.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


min_found = 200
max_tried = 1.5*min_found


def sample_background_using_ds9_regions(hdu, sky_regions):

    wcs = astWCS.WCS(hdu.header, mode='pyfits')
    pixelscale = wcs.getPixelSizeDeg() * 3600.

    data = hdu.data

    center_xy = wcs.wcs2pix(sky_regions[:,0], sky_regions[:,1])
    # print center_xy
    center_xy = numpy.array(center_xy)
    
    cx = center_xy[:,0]
    cy = center_xy[:,1]
    width = sky_regions[:,2]/2. / pixelscale
    height = sky_regions[:,3]/2. / pixelscale

    in_ota = ((cx + width) > 0) & ((cx - width) < data.shape[1]) & \
             ((cy + height) > 0) & ((cy - height) < data.shape[0])
    
    cx = cx[in_ota]
    cy = cy[in_ota]
    w = width[in_ota]
    h = height[in_ota]

    if (cx.size <= 0):
        # no boxes in this OTA
        return None
        
    left = numpy.floor(cx - w).astype(int)
    right = numpy.ceil(cx + w).astype(int)
    top = numpy.ceil(cy + h).astype(int)
    bottom = numpy.floor(cy - h).astype(int)

    left[left < 0] = 0
    bottom[bottom < 0] = 0
        
    results = []
    for box in range(cx.shape[0]):

        cutout = data[bottom[box]:top[box], left[box]:right[box]]
        median = bottleneck.nanmedian(cutout.astype(numpy.float32))
        if (numpy.isfinite(median)):
            results.append([cx[box], cy[box], median])

    #print results
    if (len(results) <= 0):
        return None

    return numpy.array(results)
        
        


def sample_background(data, wcs, starcat, min_found=200, boxwidth=30, 
                      fit_regions=None, box_center=None,
                      min_box_spacing=5,
                      combine_method=bottleneck.nanmedian):

    # Now pick a number of random data points, and keep 
    # searching until we either found 50 per OTA or have tried 100 times
    found = 0
    tried = 0
    max_tried = int(1.5*min_found)

    if (fit_regions is None):
        fit_regions = []

    skip_nan_boxes = True
    if (box_center is None):
        box_center = numpy.zeros(shape=(max_tried,2))
        box_center[:,0] = numpy.random.randint(boxwidth, data.shape[1]-boxwidth, max_tried)
        box_center[:,1] = numpy.random.randint(boxwidth, data.shape[0]-boxwidth, max_tried)
    else:
        min_found = max_tried = box_center.shape[0]
        skip_nan_boxes = False

    # Unpack the x/y coordinates of all known stars/sources in this frame
    if (starcat is not None):
        ota_x, ota_y = starcat

    #
    # Now check the randomly selected regions
    #
    while (found < min_found and tried < max_tried):
        #print box_center[tried,:]

        x1, x2 = int(box_center[tried,0]-boxwidth), int(box_center[tried,0]+boxwidth)
        y1, y2 = int(box_center[tried,1]-boxwidth), int(box_center[tried,1]+boxwidth)

        cutout = numpy.array(data[y1:y2,x1:x2], dtype=numpy.float32)
        #cutout = data[y1:y2,x1:x2]
        if ((numpy.sum(numpy.isfinite(cutout)) != cutout.shape[0]*cutout.shape[1]) and skip_nan_boxes):
            # Contains an illegal value
            tried += 1
            continue

        min_distance = 0
        if (starcat != None and skip_nan_boxes):
            # Check if there's a star in or close to this box
            star_contaminated = False
            dx = box_center[tried,1] - ota_x
            dy = box_center[tried,0] - ota_y
            dr = numpy.sqrt( dx**2 + dy**2 )
            dr_sorted = numpy.sort(dr)
            if (dr_sorted[0] < min_box_spacing*boxwidth):
                # This means there's a star nearby
                tried += 1
                continue
                pass
            min_distance = dr_sorted[0]

        #sky_level = numpy.median(cutout)
        sky_level = combine_method(cutout)

        ra, dec = 0., 0.
        if (wcs != None):
            ra, dec = wcs.pix2wcs(box_center[tried,0], box_center[tried,1])

        sky_point = [ra, dec, box_center[tried,0], box_center[tried,1], sky_level, tried, min_distance]
        fit_regions.append(sky_point)

        tried += 1
        found += 1

    return fit_regions



def fit_background(hdulist, plotname=None, exclude_videocells=True, fit_order=3, makeplots="none"):

    # First of all, get the list of sources in all the frames
    #print hdulist.info()

    try:
        odicat = hdulist['CAT.ODI']
        #print odicat.data

        odi_x = odicat.data.field("X")
        odi_y = odicat.data.field("Y")
        odi_ota = odicat.data.field("OTA")
        odi_ra = odicat.data.field("RA")
        odi_dec = odicat.data.field("DEC")
    except:
        odi_ota = numpy.zeros(shape=(0))
        odi_ra = numpy.zeros(shape=(0))
        odi_dec = numpy.zeros(shape=(0))
        odi_x = numpy.zeros(shape=(0))
        odi_y = numpy.zeros(shape=(0))

    #odi_radec = odicat.data[:,0:2]
    #print odi_radec

    # Now work out what the filter is and what OTAs' we should use for the background fit
    filter = hdulist['PRIMARY'].header['FILTER']
    #print filter

    otas_to_fit = central_3x3
    #otas_to_fit = otas_to_normalize_ff[filter]
    print("Fitting sky in OTAs",otas_to_fit)

    # Now go through each of the OTAs, and pick a number of datapoints to determine 
    # the background level. Measuring points should not be close to any stars and also
    # should not contain any NaN pixels.
    # For each region, save the median level and the sky-coordinates that we, in the last
    # step convert back into pixel coordinates

    fit_regions = []
    #min_ra, max_ra, min_dec, max_dec = numpy.nan, numpy.nan, numpy.nan, numpy.nan
    minmax_radec =  [1e9, -1e9, 1e9, -1e9] 

    all_ra, all_dec = None, None
    for ota in otas_to_fit:

        ext_name = "OTA%02d.SCI" % (ota)
        stdout_write("\rRandom-sampling background of OTA %s ..." % (ext_name))

        cellmode = hdulist[ext_name].header['CELLMODE']
        if (cellmode.find("V") != -1):
            # this OTA contains at least one video cell. This screws up the 
            # background, so let's ignore this OTA alltogether
            continue

        wcs = astWCS.WCS(hdulist[ext_name].header, mode="pyfits")
        
        # Select the X/Y coordinates of all sources in this ota
        ota_select = odi_ota == ota
        ota_x = odi_x[ota_select]
        ota_y = odi_y[ota_select]
        ota_ra = odi_ra[ota_select]
        ota_dec = odi_dec[ota_select]

        if (all_ra is not None):
            all_ra = ota_ra
            all_dec = ota_dec
        else:
            all_ra = numpy.append(all_ra, ota_ra)
            all_dec = numpy.append(all_dec, ota_dec)

        max_xy = hdulist[ext_name].data.shape

        # Determine the maximum and minimum coordinates
        minmax_radec = find_maximum_extent(minmax_radec, wcs, max_xy)

        starcat = (ota_x, ota_y) if numpy.sum(ota_select) > 0 else None
        fit_regions = sample_background(hdulist[ext_name].data, wcs, 
                                        starcat, 
                                        min_found=200, 
                                        boxwidth=20, 
                                        fit_regions=fit_regions)

        continue

        if (False):
            # Now pick a number of random data points, and keep 
            # searching until we either found 50 per OTA or have tried 100 times
            found = 0
            tried = 0
            box_center = numpy.zeros(shape=(max_tried,2))
            box_center[:,0] = numpy.random.randint(boxwidth, max_xy[0]-boxwidth, max_tried)
            box_center[:,1] = numpy.random.randint(boxwidth, max_xy[1]-boxwidth, max_tried)

            while (found < min_found and tried < max_tried):
                #print box_center[tried,:]

                x1, x2 = box_center[tried,0]-boxwidth, box_center[tried,0]+boxwidth
                y1, y2 = box_center[tried,1]-boxwidth, box_center[tried,1]+boxwidth

                cutout = hdulist[ext_name].data[y1:y2,x1:x2]
                if (not numpy.isfinite(numpy.min(cutout))):
                    # Contains an illegal value
                    tried += 1
                    continue

                # Check if there's a star in or close to this box
                star_contaminated = False
                dx = box_center[tried,1] - ota_x
                dy = box_center[tried,0] - ota_y
                dr = numpy.sqrt( dx**2 + dy**2 )
                dr_sorted = numpy.sort(dr)
                if (dr_sorted[0] < 5*boxwidth):
                    # This means there's a star nearby
                    tried += 1
                    continue
                    pass

                sky_level = numpy.median(cutout)

                ra, dec = wcs.pix2wcs(box_center[tried,0], box_center[tried,1])

                sky_point = [ra, dec, box_center[tried,0], box_center[tried,1], sky_level, tried, dr_sorted[0]]
                fit_regions.append(sky_point)

                tried += 1
                found += 1

    stdout_write(" done!\n")
    #dump = open("skyfit.dump", "w")
    #numpy.savetxt(dump, numpy.array(fit_regions))
    #dump.close()

    skypoints = numpy.array(fit_regions)

    #
    # Now we have all points, let's interpolate the grid in Ra/Dec
    #
    #print minmax_radec

    #print "#points=",skypoints.shape
    
    #print "#\n#\n#\nfitorder=",fit_order,"\n\#\n#\n"

    if (fit_order < 1):

        stdout_write("Doing some simple bg-subtraction...\n")

        skylevels = skypoints[:,4]
        valid = numpy.isfinite(skylevels)
        
        # Do some iterative sigma-clipping to get rid of outliers
        median = 0
        for repeat in range(3):
            median = numpy.median(skylevels[valid])
            low_sigma = scipy.stats.scoreatpercentile(skylevels[valid], 16)
            hi_sigma = scipy.stats.scoreatpercentile(skylevels[valid], 84)
            sigma = 0.5 * (hi_sigma - low_sigma)
            min_value = median - 3 * sigma
            max_value = median + 3 * sigma
            valid = (skylevels > min_value) & (skylevels < max_value)
            if (numpy.sum(valid) < 0.1 * skypoints.shape[0]):
                break

        hdulist_out = [hdulist[0]]
        for ota in otas_to_fit:
            ext_name = "OTA%02d.SCI" % (ota)
            stdout_write("\rOTA %s: " % (ext_name))
            cellmode = hdulist[ext_name].header['CELLMODE']
            if (cellmode.find("V") >= 0):
                continue
            hdulist[ext_name].data -= median
            hdulist[ext_name].header["BGLVLCST"] = (median, "constant background level")
            hdulist_out.append(hdulist[ext_name])
        stdout_write(" done!\n")

    else:
        #
        # Fit a polynomial to the sky-background
        # coordinates are still Ra/Dec for now
        #
        stdout_write("Creating global background map ...")
        x, y, z = skypoints[:,0], skypoints[:,1], skypoints[:,4]

        # Fit a 3rd order, 2d polynomial
        m = polyfit2d(x,y,z, order=fit_order)

        # Evaluate it on a grid...
        nx, ny = 50, 50
        xx, yy = numpy.meshgrid(numpy.linspace(x.min(), x.max(), nx), 
                             numpy.linspace(y.min(), y.max(), ny))
        zz = polyval2d(xx, yy, m)

        # Plot
        #print x.max(), x.min()
        #print y.min(), y.max()
        #print xx.min(), xx.max()
        #print zz
        #print x.min(), y.max(), x.max(), y.min()
        #matplotlib.pyplot.imshow(zz, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower')
        if (plotname != None and (makeplots=="global" or makeplots=="all")):
            vmin = zz.min()
            vmax = zz.max()
            matplotlib.pyplot.imshow(zz, 
                                     extent=(minmax_radec[0], minmax_radec[1], minmax_radec[2], minmax_radec[3]), 
                                     origin='lower',
                                     vmin=vmin, vmax=vmax)
            matplotlib.pyplot.colorbar()
            matplotlib.pyplot.scatter(x, y, c=z, linewidth=0, vmin=vmin, vmax=vmax)
            #matplotlib.pyplot.scatter(all_ra, all_dec, s=2, marker=',')
            matplotlib.pyplot.title("Global sky-background fit")
            matplotlib.pyplot.xlabel("RA [degrees]")
            matplotlib.pyplot.ylabel("DEC [degrees]")
            matplotlib.pyplot.show()
            matplotlib.pyplot.savefig(plotname+".globalskyfit.png")
            matplotlib.pyplot.close()

        stdout_write(" done!\n")

        #dump = open("skyfit.dump.fit2", "w")
        #for x in range(nx):
        #    for y in range(ny):
        #        print >>dump, xx[x,y], yy[x,y], zz[x,y]
        #dump.close()


        # 
        # Now that we have the sky as fct of Ra/Dec, convert it to x/y OTA by OTA
        #
        hdulist_out = [hdulist[0]]
        for ota in otas_to_fit:

            ext_name = "OTA%02d.SCI" % (ota)
            stdout_write("\rOTA %s: " % (ext_name))
            wcs = astWCS.WCS(hdulist[ext_name].header, mode="pyfits")
            max_xy = hdulist[ext_name].data.shape


            # Sample the x/y grid with n steps
            n = 30
            # Use 20 points in each axis, sampling the OTA in a total of 20x20=400 points
            #ota_xx, ota_yy = numpy.meshgrid(numpy.linspace(1, max_xy[0], n), 
            #                                numpy.linspace(1, max_xy[1], n))

            overlap = int(0.05 * max_xy[0])
            ota_xx, ota_yy = numpy.meshgrid(numpy.linspace(-overlap, max_xy[0]+overlap, n), 
                                    numpy.linspace(-overlap, max_xy[1]+overlap, n))

            #print "otaxx/yy.shape=",ota_xx.shape, ota_yy.shape
            #print ota_xx

            # Convert x/y into ra/dec ...
            ra, dec = numpy.zeros_like(ota_xx), numpy.zeros_like(ota_yy)
            ij = itertools.product(range(n), range(n))
            for (i,j) in ij:
                ra[i,j], dec[i,j] = wcs.pix2wcs(ota_xx[i,j], ota_yy[i,j])

            # ... and look up the points from the sky-grid
            sky_radec = polyval2d(ra, dec, m)
            #print ra[0,0], dec[0,0]

            #vmin, vmax = sky_radec.min(), sky_radec.max()

            #matplotlib.pyplot.close()
            #matplotlib.pyplot.imshow(xx, yy, zz, extent=(xx.min(), xx.max(), yy.min(), yy.max()), origin='lower')
            #matplotlib.show()

            #matplotlib.pyplot.close()

            x = open("dummy", "w")
            numpy.savetxt(x, sky_radec)
            x.close()

            #print "ra/dec/xy shapes:",ra.shape, dec.shape, sky_radec.shape, xx.shape, yy.shape

            # Now use the new 20x20 grid to interpolate the sky-background as a function of x/y
            stdout_write("Creating interpolation ...")
            f = scipy.interpolate.interp2d(ota_xx, ota_yy, sky_radec, kind='linear') #cubic')

            # Create the full resolution grid
            full_xx = numpy.linspace(1, max_xy[0], max_xy[0])
            full_yy = numpy.linspace(1, max_xy[1], max_xy[1])

            stdout_write(" performing interpolation ...")
            fullres_z = f(full_xx, full_yy)

            if (plotname != None and (makeplots=="ota" or makeplots=="all")):
                stdout_write(" plotting ...")
                #print "ra,dec, sky_radec", ra.shape, dec.shape, sky_radec.shape
                matplotlib.pyplot.scatter(ra, dec, c=sky_radec, vmin=vmin, vmax=vmax)
                #                          extent=(minmax_radec[0], minmax_radec[1], minmax_radec[2], minmax_radec[3]))
                matplotlib.pyplot.colorbar()
                matplotlib.pyplot.title("Local OTA sky-background fit")
                matplotlib.pyplot.xlabel("RA [degrees]")
                matplotlib.pyplot.ylabel("DEC [degrees]")
                matplotlib.pyplot.show(block=True)
                matplotlib.pyplot.savefig(plotname+".skyfit_"+ext_name+".png")
                matplotlib.pyplot.close()


                matplotlib.pyplot.imshow(fullres_z, vmin=vmin, vmax=vmax, origin='lower')
                matplotlib.pyplot.scatter(ota_xx, ota_yy, c=sky_radec, vmin=vmin, vmax=vmax)
                #                          extent=(minmax_radec[0], minmax_radec[1], minmax_radec[2], minmax_radec[3]))
                matplotlib.pyplot.colorbar()
                matplotlib.pyplot.title("Local OTA sky-background fit")
                matplotlib.pyplot.xlabel("X [pixels]")
                matplotlib.pyplot.ylabel("Y [pixels]")
                matplotlib.pyplot.show(block=True)
                matplotlib.pyplot.savefig(plotname+".skyfit_"+ext_name+"_1.png")
                matplotlib.pyplot.close()


            hdulist[ext_name].data -= fullres_z
            stdout_write(" done!\n")
            hdulist_out.append(hdulist[ext_name])


    hdu = pyfits.HDUList(hdulist_out)
    print ("done!")
    return hdu
            


if __name__ == "__main__":

    fit_order = int(cmdline_arg_set_or_default("-fitorder",3))
    plotting = cmdline_arg_set_or_default("-plot", "all")
    noclobber = cmdline_arg_isset("-noclobber")

    if (cmdline_arg_isset("-multi")):

        filelist = get_clean_cmdline()[1:]
        #print filelist
        for filename in filelist:
            #filename = sys.argv[1]

            if (not os.path.isfile(filename)):
                continue

            hdulist = pyfits.open(filename)
            plotname = filename[:-5]
            outfile = filename[:-5]+".skysub.fits"
            if (noclobber and os.path.isfile(outfile)):
                stdout_write("%s already exists, skipping ...\n" % (outfile))
            else:
                stdout_write("#########################\n")
                stdout_write("#\n# Sky-sub for frame %s\n#\n" % filename)
                stdout_write("#########################\n")
                hdu_out = fit_background(hdulist, plotname, fit_order=fit_order, makeplots=plotting)
                stdout_write("Writing output file %s ..." % (outfile))
                clobberfile(outfile)
                hdu_out.writeto(outfile, overwrite=True)
                stdout_write(" done!\n")

    elif (cmdline_arg_isset("-sample")):
        fitsfile = get_clean_cmdline()[1]
        output_base = get_clean_cmdline()[2]

        n_samples = int(cmdline_arg_set_or_default('-nsamples', 750))
        boxsize = int(cmdline_arg_set_or_default('-boxsize', 10))
        hdulist = pyfits.open(fitsfile)
        try:
            src_hdu = hdulist['CAT.ODI']
            src_catalog = src_hdu.data
        except:
            src_catalog = None

        # print src_catalog
        for i in range(len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue

            ota = hdulist[i].header['OTA']
            print(ota)
            if (nsrc_catalog is not None):
                src_this_ota = src_catalog.field('OTA') == ota
                #src_cat = numpy.zeros(shape=(numpy.sum(src_this_ota),2))
                #src_cat[:,0] = src_catalog.field('X')[src_this_ota]
                #src_cat[:,1] = src_catalog.field('Y')[src_this_ota]
                src_cat = (src_catalog.field('X')[src_this_ota], 
                           src_catalog.field('Y')[src_this_ota])
            else:
                src_cat = None


            bgsample = sample_background(data=hdulist[i].data, wcs=None, 
                                         starcat=src_cat, min_found=n_samples, 
                                         boxwidth=boxsize, 
                                         fit_regions=[], box_center=None,
                                         combine_method=numpy.average)
            numpy.savetxt(output_base+".OTA%02d" % ota, bgsample)
            
    else:
        filename = sys.argv[1]
        outfile = sys.argv[2]

        hdulist = pyfits.open(filename)
        plotname = filename[:-5]
        hdu_out = fit_background(hdulist, plotname, fit_order=fit_order, makeplots=plotting)
        stdout_write("Writing output file %s ..." % (outfile))
        hdu_out.writeto(outfile, overwrite=True)
        stdout_write(" done!\n")
