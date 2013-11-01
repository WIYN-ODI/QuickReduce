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

    # print hdr
    xi = numpy.zeros(shape=(6,6))
    xi_r = numpy.zeros(shape=(1,6))

    eta = numpy.zeros(shape=(6,6))
    eta_r = numpy.zeros(shape=(1,6))

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



def update_etaxi(etaxi, etaxi_r, p):

    skip = (ordering  <  0) | (ordering >= p.shape[0])
    valid = (ordering >= 0) & (ordering <  p.shape[0])

    order2d = ordering.copy()
    order2d[skip] = 0

    etaxi_change = p[order2d]
    etaxi[valid] = etaxi_change[valid]



    skip_r  = (ordering_r <  0) | (ordering_r >= p.shape[0])
    valid_r = (ordering_r >= 0) | (ordering_r <  p.shape[0])
    order1d = ordering_r.copy()
    order1d[skip_r] = 0

    etaxi_r_change = p[order1d]
    etaxi_r[valid_r] = etaxi_r_change[valid_r]

    return etaxi, etaxi_r


def update_polynomial(wcs_polynomials, px, py):

    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_polynomials

    xi, xi_r = update_etaxi(xi, xi_r, px)
    eta, eta_r = update_etaxi(eta, eta_r, py)

    return xi, xi_r, eta, eta_r, cd, crval, crpix

