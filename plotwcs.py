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

    
    for idx, ext in enumerate(hdulist):
        if (not type(ext)== pyfits.hdu.image.ImageHDU):
            continue

        # Read input WCS
        ext.header['NAXIS'] = 2
        ext.header['NAXIS1'] = 4096
        ext.header['NAXIS2'] = 4096
        ext.header['CRVAL1'] += wcs_offset[0]
        ext.header['CRVAL2'] += wcs_offset[1]


        wcs = astWCS.WCS(ext.header, mode='pyfits')

        # compute coords of the 4 corners
        
        radec = numpy.array(wcs.pix2wcs(corners[:,0], corners[:,1]))
        radec -= wcs_offset
        #print radec, radec.shape

        # coll = matplotlib.collections.PolyCollection(
        #     radec,
        #     facecolor='none',edgecolor='#808080', 
        #     linestyle='-')

        # ax.add_collection(coll)

        ax.plot(radec[:,0], radec[:,1], "b-")
        
        center = numpy.average(radec[:-1,:], axis=0)
        #print center
        ax.text(center[0], center[1], ext.name, ha='center', va='center', color='black')

        #break

    fig.show()
    matplotlib.pyplot.show()

