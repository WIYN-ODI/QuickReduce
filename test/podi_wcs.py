#!/usr/local/bin/python

import numpy
import pyfits
import sys

def pix2sky(x, y, hdr):
    # Figure out what type of coordinate system we have
    
    wcs_type = hdr['CTYPE1'][-3:]
    if (wcs_type in ('TAN', 'TPV')):
        # Ok, we got this one
        return pix2sky_tan(x, y, hdr)
    else:
        print "This WCS system is yet unsupported"
        return None
        

def read_pvs(hdr, axis):

    pvs = numpy.zeros(shape=(40))

    for i in range(0, 40):
        name = "PV%d_%d" % (axis, i)
        if (name in hdr): pvs[i] = hdr[name]

    name0 = "PV%d_1" % (axis)
    pvs[1] = hdr[name0] if name0 in hdr else 1.

    print pvs
    
    return pvs

def compute_pvs(pv, xi, eta):

    r = numpy.sqrt(xi*xi + eta*eta)

    xi_out = pv[0] \
        + pv[1]*xi + pv[2]*eta + pv[3]*r \
        \
        + pv[4]*xi*xi + pv[5]*xi*eta + pv[6]*eta*eta \
        \
        + pv[7]*xi*xi*xi + pv[8]*xi*xi*eta + pv[9]*xi*eta*eta + pv[10]*eta*eta*eta + pv[11]*r*r*r \
        \
        + pv[12]*xi*xi*xi*xi + pv[13]*xi*xi*xi*eta + pv[14]*xi*xi*eta*eta + pv[15]*xi*eta*eta*eta + pv[16]*eta*eta*eta*eta \
        \
        + pv[17]*xi*xi*xi*xi*xi + pv[18]*xi*xi*xi*xi*eta + pv[19]*xi*xi*xi*eta*eta + pv[20]*xi*xi*eta*eta*eta \
        + pv[21]*xi*eta*eta*eta*eta + pv[22]*eta*eta*eta*eta*eta + pv[23]*r*r*r*r*r \
        \
        + pv[24]*xi*xi*xi*xi*xi*xi + pv[25]*xi*xi*xi*xi*xi*eta + pv[26]*xi*xi*xi*xi*eta*eta + pv[27]*xi*xi*xi*eta*eta*eta \
        + pv[28]*xi*xi*eta*eta*eta*eta + pv[29]*xi*eta*eta*eta*eta*eta + pv[30]*eta*eta*eta*eta*eta*eta \
        \
        + pv[31]*xi*xi*xi*xi*xi*xi*xi + pv[32]*xi*xi*xi*xi*xi*xi*eta + pv[33]*xi*xi*xi*xi*xi*eta*eta \
        + pv[34]*xi*xi*xi*xi*eta*eta*eta + pv[35]*xi*xi*xi*eta*eta*eta*eta + pv[36]*xi*xi*eta*eta*eta*eta*eta \
        + pv[37]*xi*eta*eta*eta*eta*eta*eta + pv[38]*eta*eta*eta*eta*eta*eta*eta + pv[39]*r*r*r*r*r*r*r
    
    return xi_out
    

def pix2sky_tan(x, y, hdr):

    cd_11 = hdr['CD1_1'] if 'CD1_1' in hdr else 0.
    cd_12 = hdr['CD1_2'] if 'CD1_2' in hdr else 0.
    cd_21 = hdr['CD2_1'] if 'CD2_1' in hdr else 0.
    cd_22 = hdr['CD2_2'] if 'CD2_2' in hdr else 0.

    crpix1 = hdr['CRPIX1'] if 'CRPIX1' in hdr else hdr['NAXIS1']/2
    print "CRPIX-check=",hdr['CRPIX1'], crpix1
    crpix2 = hdr['CRPIX2'] if 'CRPIX2' in hdr else hdr['NAXIS2']/2

    xi  = cd_11 * (x - crpix1) + cd_12 * (y - crpix2)
    eta = cd_21 * (x - crpix1) + cd_22 * (y - crpix2)

    pvs = [None]*3
    pvs[1] = read_pvs(hdr, 1)
    pvs[2] = read_pvs(hdr, 2)
    
    #print pvs

    xi_  = compute_pvs(pvs[1], xi, eta)
    eta_ = compute_pvs(pvs[2], eta, xi)
    
    #print xi, eta

    ra  = xi_  + hdr['CRVAL1']
    dec = eta_ + hdr['CRVAL2']
    
    return ra, dec
    
if __name__ == "__main__":
    hdulist = pyfits.open(sys.argv[1])

    x = numpy.random.randint(4000, size=2)
    y = numpy.random.randint(4000, size=2)

    pix2sky(x, y, hdulist[1].header)
    
