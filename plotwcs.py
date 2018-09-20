#!/usr/bin/env python


import sys
import numpy
import os
from podi_definitions import *
import pyfits
import datetime

from astLib import astWCS
from podi_wcs import *
import bottleneck
import podi_logging

import matplotlib
import matplotlib.pyplot


wcs_offset = numpy.array([2,2])


corners = numpy.array([[4096, 4096],
                       [4096,    0],
                       [   0,    0],
                       [   0, 4096],
                       [4096, 4096]])
corners += [1,1]

if __name__ == "__main__":

    filename = sys.argv[1]
    hdulist = pyfits.open(filename)

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    try:
        extname = sys.argv[2]
    except:
        extname = None

    for idx, ext in enumerate(hdulist):
        if (not type(ext)== pyfits.hdu.image.ImageHDU):
            continue

        # Read input WCS
        ext.header['NAXIS'] = 2
        ext.header['NAXIS1'] = 4096
        ext.header['NAXIS2'] = 4096
        if (cmdline_arg_isset("-zero")):
            ext.header['CRVAL1'] = 0.
            ext.header['CRVAL2'] = 0.
        ext.header['CRVAL1'] += wcs_offset[0]
        ext.header['CRVAL2'] += wcs_offset[1]
        

        wcs = astWCS.WCS(ext.header, mode='pyfits')

        # compute coords of the 4 corners
        softbin = ext.header['SOFTBIN'] if 'SOFTBIN' in ext.header else 1
        # correct OTA-corner positions for software binning
        _corners = corners / softbin
        radec = numpy.array(wcs.pix2wcs(_corners[:,0], _corners[:,1]))
        radec -= wcs_offset
        if (extname is not None and ext.name == extname):
            print("\n",radec[:4,:])
        #print radec, radec.shape


        # Compute all individual cell coordinates
        x,y = numpy.indices((8,8))

        x1 = (x * 508).reshape((-1,1))
        x2 = x1+480
        y1 = (505*(7-y)).reshape((-1,1))
        y2 = y1 + 494

        # correct cell-corner positions for software binning
        _x1, _x2, _y1, _y2 = x1/softbin, x2/softbin, y1/softbin, y2/softbin

        radec_c1 = numpy.array(wcs.pix2wcs(_x1,_y1))
        radec_c2 = numpy.array(wcs.pix2wcs(_x1,_y2))
        radec_c3 = numpy.array(wcs.pix2wcs(_x2,_y2))
        radec_c4 = numpy.array(wcs.pix2wcs(_x2,_y1))

        cell_corners = numpy.empty((64,5,2))
        cell_corners[:,0,:] = radec_c1
        cell_corners[:,1,:] = radec_c2
        cell_corners[:,2,:] = radec_c3
        cell_corners[:,3,:] = radec_c4
        cell_corners[:,4,:] = radec_c1
        cell_corners -= wcs_offset

        coll = matplotlib.collections.PolyCollection(
            cell_corners,
            facecolor='none',
            edgecolor='#808080',
            linestyle='-'
        )
        ax.add_collection(coll)

        ax.plot(radec[:,0], radec[:,1], "b-")
        
        center = numpy.average(radec[:-1,:], axis=0)
        #print center
        ax.text(center[0], center[1], ext.name, ha='center', va='center', color='black')

#        break

    ax.set_title(filename)

    fig.show()
    matplotlib.pyplot.show()

