#!/usr/bin/env python


import sys
import numpy
import os
from podi_definitions import *
import pyfits
import datetime
import scipy
import scipy.stats
import scipy.ndimage

from astLib import astWCS


import matplotlib.pyplot as plot

filename = sys.argv[1]
hdulist = pyfits.open(filename)

img = hdulist[13].data


# figure out the maximum extensions of the input frame
wcs = astWCS.WCS(hdulist[13].header, mode="pyfits")

minRA, maxRA, minDec, maxDec = wcs.getImageMinMaxWCSCoords()
print minRA, maxRA
print minDec, maxDec


xd, yd = 4000, 4000

pixelscale = 0.13 / 3600.

#outputRA, outputDec = numpy.indices((4000,4000))
outputRA, outputDec = numpy.indices((xd,yd))

print "computing output Ra/dec for each pixel"
outputDec = numpy.array(outputDec, dtype=numpy.float32) * pixelscale + minDec
outputRA = numpy.array(outputRA, dtype=numpy.float32) * pixelscale / numpy.cos(numpy.radians(outputDec)) + minRA

print outputRA[:5,:5]
print outputDec[:5,:5]

#pixelcoords = numpy.zeros((4000,4000,2))

print "converting ra/dec to x/y for each output pixel"

print outputRA.ravel().shape
#ret = wcs.pix2wcs(outputRA.ravel(), outputDec.ravel())
ret = wcs.wcs2pix(outputRA.ravel(), outputDec.ravel())
#print ret

retarray = numpy.array(ret)
print retarray.shape
print retarray[0:10,:]
x = retarray[:,0]
y = retarray[:,1]

#xy_2d = numpy.reshape(retarray, (500,500,2))
#print xy_2d[0,0:10,:]

#print x[:10]
#print y[:10]

print "deprojecting"

#deproject1 = scipy.ndimage.interpolation.map_coordinates(input=img, coordinates=[[5,5],[7,7]], 
#                                            output=None, order=1, mode='constant', 
#                                            cval=0.0, prefilter=True)
#print deproject1

xy_coords = numpy.swapaxes(retarray,0,1)
print xy_coords

deproject = scipy.ndimage.interpolation.map_coordinates(input=img, coordinates=xy_coords, 
                                            output=None, order=1, mode='constant', 
                                            cval=0.0, prefilter=True)


print deproject.shape

deproject_2d = numpy.reshape(deproject, (xd,yd))

pyfits.PrimaryHDU(data=img[:xd,:yd]).writeto("input.fits", clobber=True)
primhdu = pyfits.PrimaryHDU(data=deproject_2d)
primhdu.header.update("CRPIX1", 1.0)
primhdu.header.update("CRPIX2", 1.0)
primhdu.header.update("CRVAL1", minRA)
primhdu.header.update("CRVAL2", minDec)
primhdu.header.update("CD1_1", pixelscale)
primhdu.header.update("CD2_2", pixelscale)
primhdu.header.update("CUNIT1", 'deg')
primhdu.header.update("CUNIT2", 'deg')
primhdu.header.update("CTYPE1", 'RA---TAN')
primhdu.header.update("CTYPE2", 'DEC--TAN')

primhdu.writeto("output.fits", clobber=True)
