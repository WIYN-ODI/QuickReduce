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

inputhdu = hdulist['OTA61.SCI']

img = inputhdu.data.T


# figure out the maximum extensions of the input frame
wcs = astWCS.WCS(inputhdu.header, mode="pyfits")

def fill_entries(hdr, ordering, keyformat):
    retarray = numpy.zeros(shape=ordering.shape)
    for y in range(retarray.shape[0]):
        for x in range(retarray.shape[1]):
            keyname = keyformat % (ordering[y,x])
            try:
                #retarray[y,x] = ordering[y,x] #hdr[keyname]
                retarray[y,x] = hdr[keyname]
            except:
                pass
    return retarray


def header_to_polynomial(hdr):
    
    xi = numpy.zeros(shape=(6,6))
    xi_r = numpy.zeros(shape=(1,6))

    eta = numpy.zeros(shape=(6,6))
    eta_r = numpy.zeros(shape=(1,6))

    ordering = numpy.array([
        [ 0,  1,  4,  7, 12, 17],
        [ 2,  5,  8, 13, 18, -1],
        [ 6,  9, 14, 19, -1, -1],
        [10, 15, 20, -1, -1, -1],
        [16, 21, -1, -1, -1, -1],
        [22, -1, -1, -1, -1, -1]
        ])
    ordering_r = numpy.array([
            [-1, 3, 11, -1, 23]
            ])

    xi = fill_entries(hdr, ordering, "PV1_%d")
    xi_r = fill_entries(hdr, ordering_r, "PV1_%d")
    eta = fill_entries(hdr, ordering, "PV2_%d")
    eta_r = fill_entries(hdr, ordering_r, "PV2_%d")

    if (not 'PV1_1' in hdr):  xi[0,1] = 1.0
    if (not 'PV2_1' in hdr): eta[0,1] = 1.0

    cd = numpy.zeros(shape=(2,2))
    if ('CD1_1' in hdr): cd[0,0] = hdr['CD1_1']
    if ('CD1_2' in hdr): cd[0,1] = hdr['CD1_2']
    if ('CD2_1' in hdr): cd[1,0] = hdr['CD2_1']
    if ('CD2_2' in hdr): cd[1,1] = hdr['CD2_2']

    crval = numpy.zeros(shape=(2))
    if ('CRVAL1' in hdr): crval[0] = hdr['CRVAL1']
    if ('CRVAL2' in hdr): crval[1] = hdr['CRVAL2']

    crpix = numpy.zeros(shape=(2))
    if ('CRPIX1' in hdr): crpix[0] = hdr['CRPIX1']
    if ('CRPIX2' in hdr): crpix[1] = hdr['CRPIX2']

    return xi, xi_r, eta, eta_r, cd, crval, crpix

import itertools
def polyval2d(x, y, m):
    ij = itertools.product(range(m.shape[0]), range(m.shape[1]))
    z = numpy.zeros_like(x)
    for (i,j) in ij:
        #print i,j
        z += m[i,j] * x**j * y**i
        print i,j,"-->",z
    return z

def my_pix2wcs(xy, wcs_polynomials):
        
    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcs_polynomials

    #use the input x,y and compute xi and eta first

    print "crpix=",crpix
    print "xy=\n",xy
    xy_relative = xy[:,0:2] - crpix
    print "xy relative=\n",xy_relative

    xi  = xy_relative[:,0] * cd[0,0] + xy_relative[:,1] * cd[0,1]
    eta = xy_relative[:,1] * cd[1,0] + xy_relative[:,1] * cd[1,1]

    print "xi=\n",xi
    print "eta=\n",eta

    # compute r = sqrt(xi**2 + eta**2)
    r = numpy.hypot(xi, eta)
    print "r=\n",r

    print "c_xi_r", c_xi_r
    print "c_eta_r", c_eta_r

    xi_prime = polyval2d(xi, eta, c_xi) \
        + numpy.polynomial.polynomial.polyval(r, c_xi_r[0])

    eta_prime = polyval2d(eta, xi, c_eta) \
        + numpy.polynomial.polynomial.polyval(r, c_eta_r[0])

    print "xi_prime=\n",xi_prime
    print "eta_prime=\n",eta_prime

    print "crval:",crval
    output = numpy.zeros_like(xy) #xi.shape[0],2)
    output[:,0] = xi_prime
    output[:,1] = eta_prime
    output[:,0:2] += crval


    return output


print "running soem stupid testing"
ra, dec = wcs.pix2wcs(2000.0, 2000.0)
print ra,dec
x,y = wcs.wcs2pix(ra, dec)
print x,y

wcs.header['PV2_18'] = wcs.header['PV2_18'] #0 #100e6
wcs.updateFromHeader()
x,y = wcs.wcs2pix(ra, dec)
print x,y

minRA, maxRA, minDec, maxDec = wcs.getImageMinMaxWCSCoords()
print minRA, maxRA
print minDec, maxDec


xd, yd = 4000, 4000

pixelscale = 0.11 / 3600. 

#output_data = numpy.zeros(shape=(4000,4000))

cosine_correction = math.cos(math.radians(inputhdu.header['CRVAL2']))
width = (maxRA - minRA) / pixelscale * cosine_correction
height = (maxDec - minDec) / pixelscale
print "\nnew size: %d x %d\n" % (width, height)


primhdu = pyfits.PrimaryHDU(data=numpy.zeros(shape=(width,height)))
# primhdu.header.update("CRPIX1", 0) #inputhdu.header['CRPIX1'])
# primhdu.header.update("CRPIX2", 0) #inputhdu.header['CRPIX2'])
# primhdu.header.update("CRVAL1", inputhdu.header['CRVAL1'])
# primhdu.header.update("CRVAL2", inputhdu.header['CRVAL2'])
# primhdu.header.update("CD1_1",  inputhdu.header['CD1_1'])
# primhdu.header.update("CD1_2",  inputhdu.header['CD1_2'])
# primhdu.header.update("CD2_1",  inputhdu.header['CD2_1'])
# primhdu.header.update("CD2_2",  inputhdu.header['CD2_2'])

primhdu.header.update("CRPIX1", 0.0)
primhdu.header.update("CRPIX2", 0.0)
primhdu.header.update("CRVAL1", 24.) #inputhdu.header['CRVAL1'])
primhdu.header.update("CRVAL2", 31.) #inputhdu.header['CRVAL2'])
primhdu.header.update("CD1_1", pixelscale)
primhdu.header.update("CD1_2", 0.0)
primhdu.header.update("CD2_1", 0.0)
primhdu.header.update("CD2_2", pixelscale)

primhdu.header.update("CUNIT1", 'deg')
primhdu.header.update("CUNIT2", 'deg')
primhdu.header.update("CTYPE1", 'RA---TAN')
primhdu.header.update("CTYPE2", 'DEC--TAN')
primhdu.header.update("EQUINOX", 2000.0)
primhdu.header.update("RADESYS", 'ICRS')


# Now use this preliminary coordiante system to 
# compute the position of the reference pixel
print "Creating output Ra/Dec system"
outputwcs = astWCS.WCS(primhdu.header, mode="pyfits")


#wcs_polynomials = header_to_polynomial(inputhdu.header)
wcs_polynomials = header_to_polynomial(primhdu.header)
xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_polynomials
print xi

print "=============================="
print "=============================="
print "=============================="
radec = my_pix2wcs(numpy.array([[2000.,1000.]]), wcs_polynomials)
print "my=",radec
print "wcs=",outputwcs.pix2wcs(2000., 1000.)
print "=============================="
print "=============================="
print "=============================="

sys.exit(0)



# Now for each pixel in the output frame, determine the RA and DEC sky coordinates
print "Translating x/y to Ra/Dec"
y, x = numpy.indices(primhdu.data.shape)
#print x[:10,:10]

output_ra = x * (pixelscale / cosine_correction) + minRA
output_dec = y * pixelscale + minDec

# Compute the pixel coordiante of the reference pixel for the current wcs solution
outputwcs.header['CRPIX1'] = 0
outputwcs.header['CRPIX2'] = 0
outputwcs.updateFromHeader()
crpix_x, crpix_y = outputwcs.wcs2pix(output_ra[0,0], output_dec[0,0])
print crpix_x, crpix_y

outputwcs.header['CRPIX1'] = crpix_x * -1
outputwcs.header['CRPIX2'] = crpix_y * -1
primhdu.header['CRPIX1'] = crpix_x * -1
primhdu.header['CRPIX2'] = crpix_y * -1
outputwcs.updateFromHeader()

newtry = outputwcs.wcs2pix(output_ra[0,0], output_dec[0,0])
print newtry

# origin_ra, origin_dec = wcs.pix2wcs(0,0)
# # Now compute the pixel of the original reference in the new frame
# crpix_x, crpix_y = outputwcs.wcs2pix(origin_ra, origin_dec)
# print crpix_x, crpix_y
# outputwcs.header['CRPIX1'] = crpix_x * -1
# outputwcs.header['CRPIX2'] = crpix_y * -1
# primhdu.header['CRPIX1'] = crpix_x * -1
# primhdu.header['CRPIX2'] = crpix_y * -1
# outputwcs.updateFromHeader()
# newtry = outputwcs.wcs2pix(origin_ra, origin_dec)
# print newtry

print "============================"



#output_radec = numpy.empty(shape=(output_ra.shape[0], output_ra.shape[1], 2))
#output_radec[:,:,0] = output_ra
#output_radec[:,:,1] = output_dec

# Using the sky-coordinates, convert them into X/Y in the original frame

# xy_coords_in_inputframe = numpy.array(wcs.wcs2pix(output_RaDec[:,0].ravel(), output_RaDec[:,1].ravel()))

print "converting sky-coordinates of the new frame into pixel x/y from the old frame"
xy_coords_in_inputframe = numpy.array(wcs.wcs2pix(output_ra.ravel(), output_dec.ravel()))

print xy_coords_in_inputframe
print xy_coords_in_inputframe.shape


print "swapping axes to make things compatible with numpy"

yx_coords_in_inputframe = numpy.swapaxes(xy_coords_in_inputframe,0,1)
print yx_coords_in_inputframe

print "deprojecting"

deproject = scipy.ndimage.interpolation.map_coordinates(
    input=img, coordinates=numpy.swapaxes(xy_coords_in_inputframe,0,1), 
    output=None, order=1, mode='constant', 
    cval=0.0, prefilter=True)

print deproject.shape


print "converting image back to 2-d"
deprojected_2d = numpy.reshape(deproject, primhdu.data.shape)
print deprojected_2d.shape


print "Saving image back to fits"
primhdu.data = deprojected_2d

primhdu.writeto("deprojected.fits", clobber=True)


sys.exit(0)




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
