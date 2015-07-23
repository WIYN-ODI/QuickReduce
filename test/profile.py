#!/usr/bin/env python


import pyfits
import numpy
import math
import scipy
import os, sys

_dir, _ = os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append(_dir+"/../")

import matplotlib.pyplot

from podi_commandline import *
from podi_definitions import clobberfile


def get_profile(data, center_x, center_y, mx, my, width, mode='radial', 
                normalize=False,
                dz2_limit=None):

    x, y = int(center_x), int(center_y)

    minx = numpy.max([0, x-width])
    maxx = numpy.min([data.shape[0], x+width])

    miny = numpy.max([0, y-width])
    maxy = numpy.min([data.shape[1], y+width])


    #print minx, maxx, miny, maxy

    cx = center_x - minx
    cy = center_y - miny

    cutout = data[minx:maxx, miny:maxy].copy()

    idxx, idxy = numpy.indices(cutout.shape)

    # shift coordinates by half a pixel to put coordinates into the center of each pixel
    idxx = idxx.astype(numpy.float32) + 0.5
    idxy = idxy.astype(numpy.float32) + 0.5

    #print idxx, idxy

    

    if (mode == 'radial'):
        _r = numpy.hypot(idxx-cx, idxy-cy)
        #print _r.shape

        si = numpy.argsort(_r.ravel())

        r = _r.ravel()[si]
        data = cutout.ravel()[si]

        # out = numpy.zeros((r.shape[0]*r.shape[1], 2))
        # out[:,0] = (r.ravel())[si]
        # out[:,1] = (cutout.ravel())[si]

        
        norm_factor = data[0]
        if (normalize):
            data /= norm_factor

        return r, data, norm_factor

    elif (mode == 'crosscut'):

        a = my / mx
        b = 1.
        c = -center_y - a*center_x
        #print a,b,c

        c = c + miny + minx*a
        #print "in cutout", a,b,c

        a2b2 = math.sqrt(a**2+b**2)
        distance_z = numpy.fabs(a*idxx + idxy + c)/a2b2
        #print distance_z

        _a = -1./a
        _b = 1.
        _c = -center_y - _a*center_x
        _c = _c + miny + minx*_a
        #print _a, _b, _c

        _a2b2 = math.sqrt(_a**2+_b**2)
        distance_z2 = numpy.fabs(_a*idxx + idxy + _c)/_a2b2

        norm_factor = numpy.max(cutout)
        if (normalize):
            cutout /= norm_factor

        if (not dz2_limit == None):
            within_limit = distance_z2 < dz2_limit
            return distance_z[within_limit], distance_z2[within_limit], cutout[within_limit], norm_factor

        return distance_z, distance_z2, cutout, norm_factor

if __name__ == "__main__":


    radial_mode = cmdline_arg_isset("-radial")
    profile_mode = cmdline_arg_isset("-profile")
    if (not radial_mode and not profile_mode):
        print "Need to specify either -radial or -profile"
        sys.exit(0)

    savefile = cmdline_arg_set_or_default("-save", None)
    if (not savefile == None):
        clobberfile(savefile)

    normalize = cmdline_arg_isset("-normalize")
    calibrate = cmdline_arg_isset("-calibrate")
    if (normalize and calibrate):
        print "Combining -normalize and -calibrate is not allowed!"
        sys.exit(0)

    infilename = get_clean_cmdline()[1]

    mot_x, mot_y = 0, 0
    mxy = cmdline_arg_set_or_default("-nonsidereal", None)
    if (not mxy == None):
        items = mxy.split(",")
        mot_x = float(items[0])
        mot_y = float(items[1])
        angle = numpy.arctan2(mot_x, mot_y)
        print "angle: ",angle, numpy.degrees(angle)

    hdu = pyfits.open(infilename)
    arcsec_per_pixel = hdu[0].header['CD2_2'] * 3600
    width = 15
    data = hdu[0].data.T
    magzero = hdu[0].header['MAGZERO'] if ('MAGZERO' in hdu[0].header and calibrate) else None

    # Convert zeopoint from mag/pixel to mag/arcsec
    pixelarea = arcsec_per_pixel**2
    if (not magzero == None): magzero += 2.5*math.log10(pixelarea)
    

    dz2_limit = float(cmdline_arg_set_or_default("-dz2", 0.5))
    if (not dz2_limit == None):
        dz2_limit /= arcsec_per_pixel
        print "Using dz2 limit of",dz2_limit,"pixels"

    for fxy in get_clean_cmdline()[2:]:
        items = fxy.split(",")
        print fxy, items
        fx = float(items[0]) - 1.
        fy = float(items[1]) - 1.

        # fx = float(get_clean_cmdline()[2]) - 1.
        # fy = float(get_clean_cmdline()[3]) - 1.
        # subtract 1 because fits stars counting at 1, while numpy starts at 0
        x, y = int(fx), int(fy)



        # if (len(get_clean_cmdline()) > 5):
        #     mot_x = float(get_clean_cmdline()[4])
        #     mot_y = float(get_clean_cmdline()[5])

        #     angle = numpy.arctan2(mot_x, mot_y)

        #     print numpy.degrees(angle)


        # minx = numpy.max([0, x-width])
        # maxx = numpy.min([data.shape[0], x+width])

        # miny = numpy.max([0, y-width])
        # maxy = numpy.min([data.shape[1], y+width])

        # print minx, maxx, miny, maxy

        # cx = fx - minx
        # cy = fy - miny

        # cutout = data[minx:maxx, miny:maxy]

        # idxx, idxy = numpy.indices(cutout.shape)
        # print idxx, idxy

        radial, crosscut = True, False
        # radial, crosscut = False, True

        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)

        if (radial_mode):
            # r = numpy.hypot(idxx-cx, idxy-cy)
            # print r.shape

            # si = numpy.argsort(r.ravel())

            r, cutout = get_profile(data, center_x=fx, center_y=fy, 
                                    mx=mot_x,  my=mot_y, width=15, 
                                    mode='radial', 
                                    normalize=normalize
                                )

            valid = numpy.isfinite(cutout)
            if (not magzero == None):
                _c = cutout.copy()
                _c[_c <= 0] = 100
                mag = -2.5*numpy.log10(_c) + magzero
                # mag = -2.5*numpy.log10(cutout) + magzero
                mag[cutout <= 0] = 99
                valid = cutout > 0
                cutout = mag

            if (not savefile == None):
                out = numpy.zeros((r.shape[0], 2))
                out[:,0] = r * arcsec_per_pixel
                out[:,1] = cutout
                numpy.savetxt(savefile, out)

            ax.set_xlabel("Radius [arcsec]")
            ax.set_ylabel("Intensity [flux-normalized counts]")
            ax.scatter((r*arcsec_per_pixel)[valid], cutout[valid])
            ax.set_xlim((0, 1.4*width*arcsec_per_pixel))
            ax.set_title("center position: %.2f %.2f" % (fx, fy))

        if (profile_mode):


            distance_z, distance_z2, cutout = get_profile(
                data, center_x=fx, center_y=fy, 
                mx=mot_x,  my=mot_y, width=15, 
                mode='crosscut', 
                normalize=normalize,
                dz2_limit=dz2_limit,
            )

            valid = numpy.isfinite(cutout)
            if (not magzero == None):
                _c = cutout.copy()
                _c[_c <= 0] = 100
                mag = -2.5*numpy.log10(_c) + magzero
                mag[cutout <= 0] = 99
                valid = cutout > 0
                cutout = mag

            if (not savefile == None):
                f = open(savefile, "a")
                out = numpy.zeros((distance_z.shape[0], 3))
                out[:,0] = distance_z * arcsec_per_pixel
                out[:,1] = cutout
                out[:,2] = distance_z2 * arcsec_per_pixel
                numpy.savetxt(f, out)
                print >>f, "\n\n\n\n\n"
                f.close()


            to_plot = valid
            ax.scatter((distance_z*arcsec_per_pixel)[to_plot], 
                       cutout[to_plot])
            ax.set_xlim((0, 1.4*width*arcsec_per_pixel))
            
            ax.scatter((distance_z2*arcsec_per_pixel)[to_plot], 
                       cutout[to_plot], 
                       c='red', marker='x')
        
            pass


        fig.show()
        matplotlib.pyplot.show()
