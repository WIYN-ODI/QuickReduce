#!/usr/bin/env python


import sys
import numpy
import os
from podi_definitions import *
import pyfits
import datetime

from astLib import astWCS
from podi_wcs import *


#
# This tool follows the initial part of the instructions found in 
# David Shupe, SPIE 8551, 84511M
#

def delete_distortion_headers(hdr):
    # now delete all PV headers if requested
    print "deleting headers:", len(hdr.cards)

    deleted_one = True
    while (deleted_one):
        deleted_one = False
        for card in hdr.cards:
            # print card
            key = card.keyword
            if (key[:2] == "PV"):
                # print key, key[:2]
                del hdr[key]
                deleted_one = True
                break

                # print hdr[key]
    print "done deleting headers:", len(hdr.cards)



def update_wcspoly(p, wcspoly, min_fit_order):

    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcspoly

    # Update the crpix values
    crpix[:] = p[0:2]

    # Now loop over the xi/eta matrices and update the values
    p_idx = 2
    for (i,j) in itertools.product( range(ordering.shape[0]), range(ordering.shape[1]) ):
        if ((i+j) >= min_fit_order and ordering[i,j] > 0):
            c_xi[i,j] = p[p_idx]
            p_idx += 1
    for (i,j) in itertools.product( range(ordering.shape[0]), range(ordering.shape[1]) ):
        if ((i+j) >= min_fit_order and ordering[i,j] > 0):
            c_eta[i,j] = p[p_idx]
            p_idx += 1

    wcs_poly_for_fitting = c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix

    return wcs_poly_for_fitting


def wcspoly_to_p_array(wcspoly, min_fit_order):

    c_xi, c_xi_r, c_eta, c_eta_r, cd, crval, crpix = wcspoly
    
    p = [crpix[0], crpix[1]]

    # Now loop over the xi/eta matrices and update the values
    for (i,j) in itertools.product( range(ordering.shape[0]), range(ordering.shape[1]) ):
        if ((i+j) >= min_fit_order and ordering[i,j] > 0):
            p.append(c_xi[i,j])
    for (i,j) in itertools.product( range(ordering.shape[0]), range(ordering.shape[1]) ):
        if ((i+j) >= min_fit_order and ordering[i,j] > 0):
            p.append(c_eta[i,j])

    return numpy.array(p)



def fit_wcs_model(p, xy, radec, wcspoly, min_fit_order = 2):

    # Given are the 
    # - updated CD matrix <-- fixed
    # - CRPIX <-- can change
    # - CRVAL <-- fixed
    # - polynomial coefficients <-- can change above linear
    # - radial terms are not allowed

    wcs_poly_for_fitting = update_wcspoly(p, wcspoly, min_fit_order)

    computed_radec = wcs_pix2wcs(xy, wcspoly, include_distortion=True, debug=False)

    diff = radec - computed_radec

    return diff.ravel()


if __name__ == "__main__":


    input_wcs = sys.argv[1]

    output_wcs = sys.argv[2]

    hdulist = pyfits.open(input_wcs)
    d_crval1, d_crval2 = None, None

    dummyfile = open("simplify.debug", "w")

    for ext in range(len(hdulist)):
        if (not type(hdulist[ext])== pyfits.hdu.image.ImageHDU):
            continue

        print hdulist[ext].header['EXTNAME']
        hdr = hdulist[ext].header

        wcs_in = astWCS.WCS(hdr, mode='pyfits')
        print "CHIP CENTER:",wcs_in.getCentreWCSCoords()

        wcspoly_before = header_to_polynomial(hdr)

        # Now read all values we might need
        crval_1 = hdr['CRVAL1']
        crval_2 = hdr['CRVAL2']

        crpix_1 = hdr['CRPIX1']
        crpix_2 = hdr['CRPIX2']

        cd_11 = hdr['CD1_1']
        cd_12 = hdr['CD1_2']
        cd_21 = hdr['CD2_1']
        cd_22 = hdr['CD2_2']

        # Also read the first three distortion factors if available
        pv_10, pv_20 = 0., 0.
        pv_11, pv_21 = 1., 1.,
        pv_12, pv_22 = 0., 0.,
        if ('PV1_0' in hdr): pv_10 = hdr['PV1_0']
        if ('PV1_1' in hdr): pv_11 = hdr['PV1_1']
        if ('PV1_2' in hdr): pv_12 = hdr['PV1_2']
        if ('PV2_0' in hdr): pv_20 = hdr['PV2_0']
        if ('PV2_1' in hdr): pv_21 = hdr['PV2_1']
        if ('PV2_2' in hdr): pv_22 = hdr['PV2_2']

        cd = numpy.matrix([[cd_11, cd_12], [cd_21, cd_22]])

        pv0 = numpy.array([[pv_10],[pv_20]])
        print "\nPV_0=\n", pv0

        #print cd.dot(pv_10_20)

        pv = numpy.matrix([[pv_11, pv_12], [pv_21, pv_22]])
        print "\nPV=",pv

        cd_prime = cd.dot(pv)

        print "\nCD =\n",cd
        print "\nCD'=\n",cd_prime

        # Compute the new CD matrix, folding in the lowest order PV factors
        # (11  12)'  (11  12)   (11  12)   (cd11.pv11+cd12.pv22   cd11.pv12+cd12.pv21)
        # (  CD  ) = (  CD  ) x (  PV  ) = (                                         )
        # (21  22)   (21  22)   (22  21)   (cd21.pv11+cd22.pv22   cd21.pv12+cd22.pv21)

        # product (cd) x (pv)
        cd_11_prime = cd_11 * pv_11 + cd_12 * pv_22
        cd_12_prime = cd_11 * pv_12 + cd_12 * pv_21
        cd_21_prime = cd_21 * pv_11 + cd_22 * pv_22
        cd_22_prime = cd_21 * pv_12 + cd_22 * pv_21

        # # product (pv) x (cd)
        # cd_11_prime = cd_11 * pv_11 + cd_21 * pv_12
        # cd_12_prime = cd_12 * pv_11 + cd_22 * pv_12
        # cd_21_prime = cd_11 * pv_22 + cd_21 * pv_21
        # cd_22_prime = cd_12 * pv_22 + cd_22 * pv_21

        cd_prime = numpy.matrix([[cd_11_prime, cd_12_prime], [cd_21_prime, cd_22_prime]])
        print "\nCD' (#2) =\n",cd_prime


        # invert CD matrix
        cd_inverse = cd.I

        UV_diff = cd_inverse.dot(pv0)
        print "\nUV-diff=\n",UV_diff

        # # Also compute the offset to the CRVAL values
        # crval_1_prime = crval_1 + cd_11 * pv_10 + cd_12 * pv_20
        # crval_2_prime = crval_2 + cd_21 * pv_10 + cd_22 * pv_20

        # Compute XY_diff
        XY_diff = cd_prime.dot(UV_diff)
        print "\nXY_diff=\n",XY_diff

        crval_1_shift = XY_diff[0,0]/math.cos(math.radians(crval_2))
        crval_2_shift = XY_diff[1,0]

        # Compute a delta_CRVAL from the first encountered OTA header
        d_crval1 = crval_1_shift if (d_crval1 == None) else d_crval1
        d_crval2 = crval_2_shift if (d_crval2 == None) else d_crval2

        # and apply it to the output header
        crval_1_prime = crval_1 + d_crval1 
        crval_2_prime = crval_2 + d_crval2 

        #  _
        # /!\  So far that part is largely taken directly from
        #  |   the SPIE paper, the next part is more ODI specific
        #  |
        #  |   We need to ensure all CRVAL parameters are the same
        #  |   small differences in the relative position should go 
        #  |   into differences in the CRPIX keywords.
        #

        # Now write all headers back to file
        hdr['CRVAL1'] = crval_1_prime
        hdr['CRVAL2'] = crval_2_prime

        hdr['CD1_1' ] = cd_prime[0,0]
        hdr['CD1_2' ] = cd_prime[0,1]
        hdr['CD2_1' ] = cd_prime[1,0]
        hdr['CD2_2' ] = cd_prime[1,1]

        # And reset all distortion parameters
        hdr['PV1_0' ] = 0.
        hdr['PV1_1' ] = 1.
        hdr['PV1_2' ] = 0.
        hdr['PV2_0' ] = 0.
        hdr['PV2_1' ] = 1.
        hdr['PV2_2' ] = 0.

        # Create a wcs-poly structure for the simplified WCS system
        wcspoly_simplified = header_to_polynomial(hdr)        

        # Now create a random pattern of points so we can re-fit 
        # the polynomial distortion factors

        xy = numpy.random.random((1000,2)) * 4000.
        # Compute the Ra/Dec values using the input WCS system
        radec_in = wcs_pix2wcs(xy, wcspoly_before, False)

        # The introduced rotation is likely small, so let's obtain
        # a starting value for the fit of the new coefficients from 
        # the current distortion coefficients
        min_fit_order = 2
        p_init = wcspoly_to_p_array(wcspoly_before, min_fit_order)



        # And finally, start the optimization
        # fit_wcs_model(p, xy, radec, wcspoly, min_fit_order = 2)
        print "Starting fitting the new distortion..."
        args = (xy, radec_in, wcspoly_simplified, 2)
        fit = scipy.optimize.leastsq(fit_wcs_model,
                                     p_init, 
                                     args=args, 
                                     full_output=1)
        print "Yippie, done fitting"

        p_final = fit[0]
        print "before-fit:\n",p_init
        print "after fit:\n",fit[0]

        # Apply the best-fit parameter to the simplified wcs-poly structure
        wcspoly_after = update_wcspoly(p_final, wcspoly_simplified, min_fit_order)

        # Update the header
        wcs_wcspoly_to_header(wcspoly_after, hdr)

        # And put header back into hdulist
        hdulist[ext].header = hdr

        wcs_out = astWCS.WCS(hdr, mode='pyfits')
        print "CHIP CENTER:",wcs_out.getCentreWCSCoords()
        print wcs_out.pix2wcs(4000,4000)
        # print hdr

#        wcspoly_after = header_to_polynomial(hdr)


        xy = numpy.random.random((1000,2)) * 4000
        radec_before = wcs_pix2wcs(xy, wcspoly_before)
        radec_after  = wcs_pix2wcs(xy, wcspoly_after)
        radecs = numpy.append(radec_before, radec_after, axis=1)

        numpy.savetxt(dummyfile, radecs)
        print >>dummyfile, "\n\n\n\n\n"


    hdulist.writeto(output_wcs, clobber=True)















def dummy():
    # print "CRVAL before", crval_1, crval_2
    # print "CRVAL after ", crval_1_prime, crval_2_prime

    # The above step might have over- or under-corrected CRVAL
    # Compute how much offset we need to compensate via a change to CRPIX
    d_crval1_crpix = d_crval1 - crval_1_shift
    d_crval2_crpix = d_crval2 - crval_2_shift
    print "Left to compensate via crpix:", d_crval1_crpix, d_crval2_crpix

    # Invert the new CD' matrix
    cd_prime_inverse = cd_prime.I

    # print cd_prime

    print "\nCD'^-1=\n",cd_prime_inverse

    # Just as a sanity check, the product of CD and cd_inverse should be one
    unity = cd_prime.dot(cd_prime_inverse)
    print unity

    d_crval_for_crpix = numpy.array([[d_crval1_crpix], [d_crval2_crpix]])
    print "\nD-CRVAL=\n",d_crval_for_crpix

    d_crpix = cd_prime_inverse.dot(d_crval_for_crpix)
    print "delta crpix=\n",d_crpix

    print "CRPIX (unchanged):", crpix_1, crpix_2

