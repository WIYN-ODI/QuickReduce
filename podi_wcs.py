#!/usr/bin/env python

"""

This module contains a set of routines to read, manipulaten and write the
scamp-generated WCS system in ODI frames.

"""

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
        [-1, 3, -1, 11, -1, 23]
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
        hdr['CD%d_%d' % (i+1,j+1)] = cd[i,j]

    # Write all non-radial distortion terms
    for (i,j) in itertools.product(range(ordering.shape[0]), range(ordering.shape[1])):
        x = ordering[i,j]
        if (x <= 0):
            continue
        hdr['PV1_%d' % (x)] = xi[i,j]
        hdr['PV2_%d' % (x)] = eta[i,j]

    # And finally, finish the radial terms
    for i in range(ordering_r.shape[1]):
        x = ordering_r[0,i]
        if (x <= 0):
            continue
        hdr['PV1_%d' % (x)] = xi_r[0,i]
        hdr['PV2_%d' % (x)] = eta_r[0,i]

    return


def wcs_clear_distortion(wcs_polynomials, min=0):

    print "clearing distortion"
    xi, xi_r, eta, eta_r, cd, crval, crpix = wcs_polynomials

    xi[ordering>=min] = 0
    xi_r[:] = 0

    eta[ordering>=min] = 0
    eta_r[:] = 0

    xi[0,1] = 1.
    eta[0,1] = 1
    print ordering[0,1]

    wcs_polynomials = xi, xi_r, eta, eta_r, cd, crval, crpix

    return wcs_polynomials

    

def polyval2d(x, y, m):
    ij = itertools.product(range(m.shape[0]), range(m.shape[1]))
#    ij = itertools.product(range(2), range(2))
#    numpy.savetxt(sys.stdout, m)
    # print
    # print x[:3]
    # print y[:3]
    z = numpy.zeros_like(x)
    for (i,j) in ij:
        #print i,j
        z += m[i,j] * x**j * y**i
#        print i,j," [",m[i,j],"]   -->",z[:3]
    return z

def wcs_pix2wcs(xy, wcs_polynomials, include_distortion=True, debug=False):
        
    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcs_polynomials

    #use the input x,y and compute xi and eta first
    xy_relative = xy[:,0:2] - crpix
    xi  = xy_relative[:,0] * cd[0,0] + xy_relative[:,1] * cd[0,1]
    eta = xy_relative[:,1] * cd[1,0] + xy_relative[:,1] * cd[1,1]
    r = numpy.hypot(xi, eta)

    if (debug):
        print "crpix=",crpix
        print "xy=\n",xy
        print "xy relative=\n",xy_relative
        print "xi=\n",xi
        print "eta=\n",eta

    # compute r = sqrt(xi**2 + eta**2)
    _xi = numpy.radians(xi)
    _eta = numpy.radians(eta)
    _r = numpy.hypot(_xi, _eta)

    if (debug):
        print "r=\n",_r
        print "c_xi_r", c_xi_r
        print "c_eta_r", c_eta_r
        
    # xi_prime = polyval2d(_xi, _eta, c_xi) \
    #     + numpy.polynomial.polynomial.polyval(_r, c_xi_r[0])

    # eta_prime = polyval2d(_eta, _xi, c_eta) \
    #     + numpy.polynomial.polynomial.polyval(_r, c_eta_r[0])

    if (include_distortion):
        xi_prime = polyval2d(xi, eta, c_xi) \
            + numpy.polynomial.polynomial.polyval(r, c_xi_r[0])
        eta_prime = polyval2d(eta, xi, c_eta) \
            + numpy.polynomial.polynomial.polyval(r, c_eta_r[0])

        # _xi_prime = polyval2d(_xi, _eta, c_xi) \
        #     + numpy.polynomial.polynomial.polyval(_r, c_xi_r[0])
        # _eta_prime = polyval2d(_eta, _xi, c_eta) \
        #     + numpy.polynomial.polynomial.polyval(_r, c_eta_r[0])
    else:
        xi_prime = xi
        eta_prime = eta

    # x = numpy.zeros(shape=(xy.shape[0],4))
    # x[:,0] = xi
    # x[:,1] = xi_prime
    # x[:,2] = eta
    # x[:,3] = eta_prime
    # dummy_x = open("distortion","a")
    # numpy.savetxt(dummy_x, x)
    # print >>dummy_x, "\n\n\n\n\n\n"
    # dummy_x.close()

    # Uncomment the next line to disable distortion.
    # xi_prime, eta_prime = xi, eta

    CRVAL2 = crval[1]

    cos_dec0 = math.cos(math.radians(crval[1]))
    sin_dec0 = math.sin(math.radians(crval[1]))


    _xi_prime = numpy.radians(xi_prime)
    _eta_prime = numpy.radians(eta_prime)



    # _ra = numpy.arctan2(_xi_prime,
    #                    cos_dec0 - _eta_prime*sin_dec0)
    # ra = numpy.degrees(_ra) + crval[0]
    # ra[ra < 0] += 360.
    
    # _dec_xxx = numpy.sqrt((cos_dec0 - _eta_prime*sin_dec0)**2 + _xi_prime**2)
    # _dec = numpy.arctan2(_eta_prime*cos_dec0 + sin_dec0, _dec_xxx)
    
    # dec = numpy.degrees(_dec)

    
    a1 = _xi_prime / math.cos(math.radians(CRVAL2))
    a2 = 1. - _eta_prime * math.tan(math.radians(CRVAL2))
    alpha_prime = numpy.arctan2(a1, a2)
    ra = numpy.degrees(alpha_prime) + crval[0]

    d1 = (_eta_prime + math.tan(math.radians(CRVAL2))) * numpy.cos(alpha_prime)
    d2 = 1. - _eta_prime * math.tan(math.radians(CRVAL2))
    dec = numpy.degrees(numpy.arctan2(d1, d2))


    rp = numpy.sqrt(_xi_prime**2 + _eta_prime**2)
        
    phi = numpy.arctan2(_xi_prime, -_eta_prime)
    phi[rp == 0] = 0
        
    theta = numpy.arctan2(math.radians(CRVAL2), rp)






    # if (debug):
    #     print "xi_prime=\n",xi_prime
    #     print "eta_prime=\n",eta_prime
    #     print "crval:",crval

    output = numpy.zeros_like(xy) #xi.shape[0],2)

    # r_p = (xi_prime**2 + eta_prime**2)
    
    # phi = numpy.arctan2(xi_prime, -eta_prime)

    # theta = numpy.arctan2( (180./math.pi), r_p)

    # output[:,0] = phi #xi_prime
    # output[:,1] = theta #eta_prime
    # output[:,0:2] += crval

    output[:,0] = ra
    output[:,1] = dec

    return output





def wcs_pix2wcs_astlib(xy, wcs_polynomials, include_distortion=True, debug=False):
        
    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcs_polynomials

    #use the input x,y and compute xi and eta first
    xy_relative = xy[:,0:2] - crpix

    # Apply CD matrix
    dx  = xy_relative[:,0] * cd[0,0] + xy_relative[:,1] * cd[0,1]
    dy = xy_relative[:,1] * cd[1,0] + xy_relative[:,1] * cd[1,1]

    r = numpy.sqrt(dx**2 + dy**2)

    dxr = numpy.radians(dx)
    dyr = numpy.radians(dy)
    rr = numpy.radians(r)

    if (include_distortion):
        
        dx_prime = polyval2d(dx, dy, c_xi) \
            + numpy.polynomial.polynomial.polyval(r, c_xi_r[0])
        dy_prime = polyval2d(dx, dy, c_eta) \
            + numpy.polynomial.polynomial.polyval(r, c_eta_r[0])

        # dx_prime = polyval2d(dxr, dyr, c_xi) \
        #     + numpy.polynomial.polynomial.polyval(rr, c_xi_r[0])
        # dy_prime = polyval2d(dxr, dyr, c_eta) \
        #     + numpy.polynomial.polynomial.polyval(rr, c_eta_r[0])


        dx, dy = dx_prime, dy_prime


    # l, m = dx, dy

    l = numpy.radians(dx)
    m = numpy.radians(dy)

    cos0 = math.cos(math.radians(crval[1]))
    sin0 = math.sin(math.radians(crval[1]))

    ra0 = math.radians(crval[0])
    dec0 = math.radians(crval[1])

    dect = cos0 - m * sin0

    rat = ra0 + numpy.arctan2(l, dect)

    dect = numpy.arctan( numpy.cos(rat-ra0) * (m * cos0 + sin0) / dect)

    ra_out = rat
    dec_out = dect

    ra_out[ra_out < 0] += 2 * math.pi


    output = numpy.zeros_like(xy) #xi.shape[0],2)

    output[:,0] = numpy.degrees(ra_out)
    output[:,1] = numpy.degrees(dec_out)

    return output



def wcs_pix2wcs_2(xy, wcs_polynomials, debug=False):
        
    # Now create a fits header and save all WCS keywords
    hdr = pyfits.core.Header()
    wcs_wcspoly_to_header(wcs_polynomials, hdr)
    # print hdr

    hdr['NAXIS'] = 2
    hdr['NAXIS1'] = 4096
    hdr['NAXIS2'] = 4096

    hdr['CTYPE1'] = "RA---TPV"
    hdr['CTYPE2'] = "DEC--TPV"
    hdr['EQUINOX'] = 2000.0

    wcs = astWCS.WCS(hdr, mode="pyfits")

    #print wcs.getCentreWCSCoords()

    #print "x=\n",xy[:,0]
    #print "y=\n",xy[:,1]
    
    #print "runnign the wcs.pix2wcs ..."
    radec = wcs.pix2wcs(xy[:,0], xy[:,1])
    #print radec
    radec = numpy.array(radec)
    #print radec
    #print radec.shape

    return radec


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



def wcs_apply_rotation(wcs_poly, angle, debug=False):

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



def wcs_apply_shift(wcs_poly, shift, debug=False):

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



def rotate_polynomial(poly, rotate):

    out = numpy.zeros_like(poly)

    # print "x'=", rotate[0]
    # print "y'=", rotate[1]

    # poly_coeffs = [
    #     [1],
    #     [1,2,1],
    #     [1,3,3,1],
    #     [1,4,6,4,1],
    #     [1,5,10,10,5,1],
    #     [1,6,15,20,15,6,1]
    #     ]

    rotate_x = numpy.array([[0, math.cos(math.radians(angle))],
                            [math.sin(math.radians(angle)), 0]])

    rotate_y = numpy.array([[0, math.cos(math.radians(angle))],
                            [-math.sin(math.radians(angle)), 0]])


    for (i,j) in itertools.product(range(poly.shape[0]), range(poly.shape[1])):
        print "\n\n",i,j

        new_matrix = ax_plus_by_multipower(rotate_x, i, rotate_y, j)
        
        sx, sy = new_matrix.shape[0], new_matrix.shape[1]

        out[:sx, :sy] += poly[i,j] * new_matrix

    return out

def ax_plus_by_multiply(in1, in2):

    size = in1.shape[0] + in2.shape[0] - 1
    out = numpy.zeros(shape=(size, size))

    # Loop over all items in array 1
    for (x,y) in itertools.product(range(in1.shape[0]), range(in1.shape[1])):

        # now perform multiplication
        to_add = in1[x,y] * in2

        out[x:x+in2.shape[0], y:y+in2.shape[1]] += to_add

    return out

def ax_plus_by_power(in1, power):

    if (power == 0):
        return numpy.array([[1]])
    elif (power == 1):
        return in1
    else:
        out = in1
        for i in range(power-1):
            out = ax_plus_by_multiply(out, in1)

    return out


def ax_plus_by_multipower(in1, power1, in2, power2):

    out1 = ax_plus_by_power(in1, power1)
    out2 = ax_plus_by_power(in2, power2)
    return ax_plus_by_multiply(out1, out2)


if __name__ == '__main__':


    if (cmdline_arg_isset("-rotate")):

        # poly = numpy.array([
        #         [0. ,  1.0],
        #         [1.0,  0.0]
        #         ])

        # rotation = numpy.array([
        #         [0.9, 0.1], 
        #         [0.1, 0.9]
        #         ])

        # print rotation[0]
        # result = rotate_polynomial(poly, rotation)

        # print result


        in1 = numpy.array([[0,1],[1,0]])
        in2 = numpy.array([[0,0.5],[1,0]])
        out = ax_plus_by_multiply(in1, in2)
        print out


        out3 = ax_plus_by_multiply(in1, out)
        print out3

        out4 = ax_plus_by_multiply(in1, out3)
        print out4

        out4b = ax_plus_by_multiply(out3, in1)
        print out4b

        out6 = ax_plus_by_multiply(out3, out3)
        print out6

        print "\n\n\n\n"
        print ax_plus_by_power(in1, 0)
        print ax_plus_by_power(in1, 1)
        print ax_plus_by_power(in1, 2)
        print ax_plus_by_power(in1, 3)
        print ax_plus_by_power(in1, 4)
        print ax_plus_by_power(in1, 5)

        print "\n\n"
        print ax_plus_by_multipower(in1,0,in2,0)
        print ax_plus_by_multipower(in1,1,in2,0)
        print ax_plus_by_multipower(in1,0,in2,1)
        print ax_plus_by_multipower(in1,1,in2,1)

        print "\n\n\ntrying rotation ...\n"
        rot_angle = 0
        new_x = numpy.array([[0, math.cos(math.radians(rot_angle))],
                             [math.sin(math.radians(rot_angle)), 0],
                             ])
        new_y = numpy.array([[0, math.cos(math.radians(rot_angle))],
                             [-math.sin(math.radians(rot_angle)), 0],
                             ])

        pure_x = numpy.array([[0, 1,],
                              [0, 0]])

        rotated_x = ax_plus_by_multiply(pure_x, new_x)
        rotated_y = ax_plus_by_multiply(pure_x, new_y)
        print rotated_x
        print rotated_y

        sys.exit(0)

    print "# hello!"

    fitsfile = sys.argv[1]

    catfile = sys.argv[2]

    cat = numpy.loadtxt(catfile)

    xy = cat[:,2:4]
    
    
    # Load the fits file and get the header
    hdulist = pyfits.open(fitsfile)
    wcspoly = header_to_polynomial(hdulist[0].header)
    
    radec_cat = cat[:,0:2]

    radec_computed = wcs_pix2wcs(xy, wcspoly, include_distortion=True)

    radec_computed2 = wcs_pix2wcs_astlib(xy, wcspoly, include_distortion=True)

    output = numpy.append(radec_cat, radec_computed, axis=1)
    output = numpy.append(output, radec_computed2, axis=1)

    # print output.shape
    numpy.savetxt(sys.stdout, output)
