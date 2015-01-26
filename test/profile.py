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


def get_profile(data, center_x, center_y, mx, my, width, mode='radial', normalize=False):

    x, y = int(center_x), int(center_y)
    width = 15

    minx = numpy.max([0, x-width])
    maxx = numpy.min([data.shape[0], x+width])

    miny = numpy.max([0, y-width])
    maxy = numpy.min([data.shape[1], y+width])


    print minx, maxx, miny, maxy

    cx = center_x - minx
    cy = center_y - miny

    cutout = data[minx:maxx, miny:maxy].copy()

    idxx, idxy = numpy.indices(cutout.shape)
    print idxx, idxy



    if (mode == 'radial'):
        _r = numpy.hypot(idxx-cx, idxy-cy)
        print _r.shape

        si = numpy.argsort(_r.ravel())

        r = _r.ravel()[si]
        data = cutout.ravel()[si]

        # out = numpy.zeros((r.shape[0]*r.shape[1], 2))
        # out[:,0] = (r.ravel())[si]
        # out[:,1] = (cutout.ravel())[si]

        
        if (normalize):
            data /= data[0]

        return r, data

    elif (mode == 'crosscut'):

        a = my / mx
        b = 1.
        c = -fy - a*fx
        print a,b,c

        c = c + miny + minx*a
        print "in cutout", a,b,c

        a2b2 = math.sqrt(a**2+b**2)
        distance_z = numpy.fabs(a*idxx + idxy + c)/a2b2
        print distance_z

        _a = -1./a
        _b = 1.
        _c = -fy - _a*fx
        _c = _c + miny + minx*_a
        print _a, _b, _c

        _a2b2 = math.sqrt(_a**2+_b**2)
        distance_z2 = numpy.fabs(_a*idxx + idxy + _c)/_a2b2

        if (normalize):
            cutout /= numpy.max(cutout)

        return distance_z, distance_z2, cutout

if __name__ == "__main__":


    radial_mode = cmdline_arg_isset("-radial")
    profile_mode = cmdline_arg_isset("-profile")
    if (not radial_mode and not profile_mode):
        print "Need to specify either -radial or -profile"
        sys.exit(0)
    savefile = cmdline_arg_set_or_default("-save", None)

    normalize = cmdline_arg_isset("-normalize")

    infilename = get_clean_cmdline()[1]

    mot_x, mot_y = 0, 0
    mxy = cmdline_arg_set_or_default("-nonsidereal", None)
    if (not mxy == None):
        items = mxy.split(",")
        mot_x = float(items[0])
        mot_y = float(items[1])

    hdu = pyfits.open(infilename)

    fxy = get_clean_cmdline()[2]
    items = fxy.split(",")
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

    arcsec_per_pixel = hdu[0].header['CD2_2'] * 3600

    width = 15

    # transpose so pixel order is x,y
    data = hdu[0].data.T

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
                                mode='radial', normalize=normalize)

        if (not savefile == None):
            out = numpy.zeros((r.shape[0], 2))
            out[:,0] = (r)
            out[:,1] = (cutout)
            numpy.savetxt(savefile, out)

        ax.set_xlabel("Radius [arcsec]")
        ax.set_ylabel("Intensity [flux-normalized counts]")
        ax.scatter(r*arcsec_per_pixel, cutout)
        ax.set_xlim((0, 1.4*width*arcsec_per_pixel))


    if (profile_mode):


        distance_z, distance_z2, cutout = get_profile(
            data, center_x=fx, center_y=fy, 
            mx=mot_x,  my=mot_y, width=15, 
            mode='crosscut', normalize=normalize,
        )

        # # find equation for line along the trail
        # # must fulfill condition ax + by + c = 0
        # a = mot_y / mot_x
        # b = 1.
        # c = -fy - a*fx
        # print a,b,c

        # c = c + miny + minx*a
        # print "in cutout", a,b,c

        # a2b2 = math.sqrt(a**2+b**2)
        # distance_z = numpy.fabs(a*idxx + idxy + c)/a2b2
        # print distance_z

        # _a = -1./a
        # _b = 1.
        # _c = -fy - _a*fx
        # _c = _c + miny + minx*_a
        # print _a, _b, _c

        # _a2b2 = math.sqrt(_a**2+_b**2)
        # distance_z2 = numpy.fabs(_a*idxx + idxy + _c)/_a2b2

        if (not savefile == None):
            out = numpy.zeros((distance_z.shape[0]*distance_z.shape[1], 3))
            out[:,0] = (distance_z.ravel())
            out[:,1] = (distance_z2.ravel())
            out[:,2] = (cutout.ravel())
            numpy.savetxt(savefile, out)


        ax.scatter((distance_z*arcsec_per_pixel)[distance_z2<2], cutout[distance_z2<2])
        ax.set_xlim((0, 1.4*width*arcsec_per_pixel))

        ax.scatter((distance_z2*arcsec_per_pixel)[distance_z<2], cutout[distance_z<2], c='red', marker='x')
        pass


    fig.show()
    matplotlib.pyplot.show()
