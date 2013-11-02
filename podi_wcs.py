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
import itertools

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


def wcs_wcspoly_to_header(wcs_polynomials, hdr):

    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_polynomials

    # write CRPIX and CRVALs
    for i in range(2):
        hdr['CRPIX%d' % (i+1)] = crpix[i]
        hdr['CRVAL%d' % (i+1)] = crval[i]

    # write CDx_y
    for (i,j) in itertools.product(range(2),range(2)):
        hdr['CD%d_%d' % (i,j)] = cd[i,j]

    # Write all non-radial distortion terms
    for (i,j) in itertools.product(range(ordering.shape[0]), range(ordering.shape[1])):
        x = ordering[i,j]
        hdr['PV1_%d' % (x)] = xi[i,j]
        hdr['PV2_%d' % (x)] = eta[i,j]

    # And finally, finish the radial terms
    for i in range(ordering_r.shape[1]):
        x = ordering_r[0,i]
        hdr['PV1_%d' % (x)] = xi_r[0,i]
        hdr['PV2_%d' % (x)] = eta_r[0,i]

    return

def polyval2d(x, y, m):
    ij = itertools.product(range(m.shape[0]), range(m.shape[1]))
    z = numpy.zeros_like(x)
    for (i,j) in ij:
        #print i,j
        z += m[i,j] * x**j * y**i
        # print i,j,"-->",z
    return z

def wcs_pix2wcs(xy, wcs_polynomials, debug=False):
        
    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcs_polynomials

    #use the input x,y and compute xi and eta first
    xy_relative = xy[:,0:2] - crpix
    xi  = xy_relative[:,0] * cd[0,0] + xy_relative[:,1] * cd[0,1]
    eta = xy_relative[:,1] * cd[1,0] + xy_relative[:,1] * cd[1,1]

    if (debug):
        print "crpix=",crpix
        print "xy=\n",xy
        print "xy relative=\n",xy_relative
        print "xi=\n",xi
        print "eta=\n",eta

    # compute r = sqrt(xi**2 + eta**2)
    r = numpy.hypot(xi, eta)

    if (debug):
        print "r=\n",r
        print "c_xi_r", c_xi_r
        print "c_eta_r", c_eta_r

    xi_prime = polyval2d(xi, eta, c_xi) \
        + numpy.polynomial.polynomial.polyval(r, c_xi_r[0])

    eta_prime = polyval2d(eta, xi, c_eta) \
        + numpy.polynomial.polynomial.polyval(r, c_eta_r[0])

    if (debug):
        print "xi_prime=\n",xi_prime
        print "eta_prime=\n",eta_prime
        print "crval:",crval

    output = numpy.zeros_like(xy) #xi.shape[0],2)
    output[:,0] = xi_prime
    output[:,1] = eta_prime
    output[:,0:2] += crval

    return output



def update_etaxi(etaxi, etaxi_r, p, debug=False):

    if (debug):
        print "\n\n\nUpdateing ext/xi:\n"
        print "eta/xi="
        numpy.savetxt(sys.stdout, etaxi, "%9.2e")
        print "eta/xi radial="
        numpy.savetxt(sys.stdout, etaxi_r, "%9.2e")

    skip = (ordering  <  0) | (ordering >= p.shape[0])
    valid = (ordering >= 0) & (ordering <  p.shape[0])

    order2d = ordering.copy()
    order2d[skip] = 0

    etaxi_change = p[order2d]
    etaxi[valid] = etaxi_change[valid]



    skip_r  = (ordering_r <  0) | (ordering_r >= p.shape[0])
    valid_r = (ordering_r >= 0) & (ordering_r <  p.shape[0])
    order1d = ordering_r.copy()
    order1d[skip_r] = 0

    etaxi_r_change = p[order1d]
    etaxi_r[valid_r] = etaxi_r_change[valid_r]

    if (debug):
        print "order_1d=\n",order1d
        print "skip_1d=\n",skip_r
        print "valid_1d=\n",valid_r

        print "eta/xi radial change="
        numpy.savetxt(sys.stdout, etaxi_r_change, "%9.2e")

        print "eta/xi radial after="
        numpy.savetxt(sys.stdout, etaxi_r, "%9.2e")

    return etaxi, etaxi_r


def update_polynomial(wcs_polynomials, px, py):

    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_polynomials

    xi, xi_r = update_etaxi(xi, xi_r, px)
    eta, eta_r = update_etaxi(eta, eta_r, py)

    return xi, xi_r, eta, eta_r, cd, crval, crpix



def wcs_poly_to_arrays(wcs_poly, debug=False):


    n_dim_xy = numpy.max(ordering)
    n_dim_r  = numpy.max(ordering_r)
    n_dim = numpy.max([n_dim_xy, n_dim_r])

    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_poly

    xi_1d = numpy.zeros(shape=(n_dim+1))
    eta_1d = numpy.zeros(shape=(n_dim+1))

    for xy in itertools.product(range(ordering.shape[0]), range(ordering.shape[1])):
        x, y = xy
        # print xy
        
        i = ordering[x,y]
        xi_1d[i] = xi[x,y]
        eta_1d[i] = eta[x,y]

    for x in range(xi_r.shape[0]):
        i = ordering_r[x]
        xi_1d[i] = xi_r[x]
        eta_1d[i] = eta_r[x]



#    for ((x,y) for x in A for y in B)

    if (debug):
        print "wcs_poly_to_arrays -> n_dim =",n_dim
        print "xi-2d:"
        numpy.savetxt(sys.stdout, xi, "%9.2e")
        print "xi-1d:"
        numpy.savetxt(sys.stdout, xi_1d, "%9.2e")

    return xi_1d, eta_1d



def wcs_apply_rotation(wcs_poly, angle, debug=True):

    # Extract the individual elements
    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_poly

    # Rotation only affects the CD matrix
    
    # compute rotation matrix
    angle_rad = numpy.radians(angle)
    cos_angle = math.cos(angle_rad)
    sin_angle = math.sin(angle_rad)

    rot_matrix = numpy.array([
            [cos_angle, -sin_angle],
            [sin_angle, cos_angle]
            ])

    rotated_cd = cd.dot(rot_matrix)

    if (debug):
        print "input cd matrix"
        numpy.savetxt(sys.stdout, cd, "%+13.6e")
        print "rotation matrix"
        numpy.savetxt(sys.stdout, rot_matrix, "%+13.6e")
        print "input cd matrix"
        numpy.savetxt(sys.stdout, rotated_cd, "%+13.6e")

        print 

    return xi, xi_r, eta, eta_r, rotated_cd, crval, crpix



def wcs_apply_shift(wcs_poly, shift, debug=True):

    # Extract the individual elements
    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_poly

    # Rotation only affects the crval values
    crval_shifted = crval + shift

    if (debug):
        print "input crval matrix"
        numpy.savetxt(sys.stdout, crval, "%+13.6e")

        print "shift"
        numpy.savetxt(sys.stdout, shift, "%+13.6e")
        print "output crval"
        numpy.savetxt(sys.stdout, crval_shifted, "%+13.6e")

        print

    return xi, xi_r, eta, eta_r, cd, crval_shifted, crpix
