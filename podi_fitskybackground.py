#! /usr/bin/env python

#
# (c) Ralf Kotulla for WIYN/pODI
#

import sys
import os
import pyfits
import numpy
import scipy
import scipy.interpolate
import math

from podi_definitions import *
from podi_plotting import *

from astLib import astWCS

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

def fit_background(hdulist, plotname=None, exclude_videocells=True, fit_order=3, makeplots="none"):

    # First of all, get the list of sources in all the frames
    #print hdulist.info()

    odicat = hdulist['CAT.ODI']
    #print odicat.data

    odi_x = odicat.data.field("X")
    odi_y = odicat.data.field("Y")
    odi_ota = odicat.data.field("OTA")
    odi_ra = odicat.data.field("RA")
    odi_dec = odicat.data.field("DEC")

    #odi_radec = odicat.data[:,0:2]
    #print odi_radec

    # Now work out what the filter is and what OTAs' we should use for the background fit
    filter = hdulist['PRIMARY'].header['FILTER']
    #print filter

    otas_to_fit = central_3x3
    #otas_to_fit = otas_to_normalize_ff[filter]
    print "Fitting sky in OTAs",otas_to_fit

    # Now go through each of the OTAs, and pick a number of datapoints to determine 
    # the background level. Measuring points should not be close to any stars and also
    # should not contain any NaN pixels.
    # For each region, save the median level and the sky-coordinates that we, in the last
    # step convert back into pixel coordinates

    fit_regions = []
    #min_ra, max_ra, min_dec, max_dec = numpy.NaN, numpy.NaN, numpy.NaN, numpy.NaN
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

        if (all_ra == None):
            all_ra = ota_ra
            all_dec = ota_dec
        else:
            all_ra = numpy.append(all_ra, ota_ra)
            all_dec = numpy.append(all_dec, ota_dec)

        max_xy = hdulist[ext_name].data.shape

        # Determine the maximum and minimum coordinates
        minmax_radec = find_maximum_extent(minmax_radec, wcs, max_xy)

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
    

    #
    # Fit a polynomial to the sky-background
    # coordinates are still Ra/Dec for now
    #
    stdout_write("Creating global background map ...")
    x, y, z = skypoints[:,0], skypoints[:,1], skypoints[:,4]
    
    # Fit a 3rd order, 2d polynomial
    m = polyfit2d(x,y,z, order=3)

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
    print  "done!"
    return hdu
            


if __name__ == "__main__":

    fit_order = cmdline_arg_set_or_default("-fitorder",3)
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
                hdu_out.writeto(outfile, clobber=True)
                stdout_write(" done!\n")
    else:
        filename = sys.argv[1]
        outfile = sys.argv[2]

        hdulist = pyfits.open(filename)
        plotname = filename[:-5]
        hdu_out = fit_background(hdulist, plotname, fit_order=fit_order, makeplots=plotting)
        stdout_write("Writing output file %s ..." % (outfile))
        hdu_out.writeto(outfile, clobber=True)
        stdout_write(" done!\n")
