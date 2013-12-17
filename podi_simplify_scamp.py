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

import dev_ccmatch

if __name__ == "__main__":


    input_wcs = sys.argv[1]

    output_wcs = sys.argv[2]

    hdulist = pyfits.open(input_wcs)
    d_crval1, d_crval2 = None, None

    dummyfile = open("simplify.debug", "w")
    n_stars = 5000

    for ext in range(len(hdulist)):
        if (not type(hdulist[ext])== pyfits.hdu.image.ImageHDU):
            continue

        imghdu = pyfits.ImageHDU(header=hdulist[ext].header, data=numpy.array((4000,4000)))
        imghdu.header['NAXIS'] = 2
        imghdu.header['NAXIS1'] = 4000
        imghdu.header['NAXIS2'] = 4000
        imghdu.header['CTYPE1'] = "RA---TPV"
        imghdu.header['CTYPE2'] = "DEC--TPV"
        imghdu.header['CRVAL1'] += 70
        imghdu.header['CRVAL2'] += 30

        newhdu = pyfits.ImageHDU(header=hdulist[ext].header, data=numpy.array((4000,4000)))
        newhdu.header['NAXIS'] = 2
        newhdu.header['NAXIS1'] = 4000
        newhdu.header['NAXIS2'] = 4000
        newhdu.header['CTYPE1'] = "RA---TPV"
        newhdu.header['CTYPE2'] = "DEC--TPV"
        newhdu.header['CRVAL1'] += 70
        newhdu.header['CRVAL2'] += 30

        #print imghdu.header
        #print "="*100

        hdr = hdulist[ext].header

        extname = hdr['EXTNAME']
        print extname

        #hdr['NAXIS'] = 2
        #hdr['NAXIS1'] = 4000
        #hdr['NAXIS2'] = 4000
        
        wcs_in = astWCS.WCS(imghdu.header, mode='pyfits')

        # Create a bunch of random x/y coordinates
        xy = numpy.random.random((n_stars,2)) * 4000.
        print xy.shape
        print xy[:5,:]


        # From the x/y, compute Ra/decs
        radec = numpy.array(wcs_in.pix2wcs(xy[:,0]-1, xy[:,1]-1))

        # Now create a new WCS, with the most simple PV distortion 
        # parameters disabled
        newhdu.header['PV1_0'] = 0.
        newhdu.header['PV1_1'] = 1.
        newhdu.header['PV1_2'] = 0.
        newhdu.header['PV2_0'] = 0.
        newhdu.header['PV2_1'] = 1.
        newhdu.header['PV2_2'] = 0.

        # if ('PV1_0' in newhdu.header): del newhdu.header['PV1_0']
        # if ('PV1_1' in newhdu.header): del newhdu.header['PV1_1']
        # if ('PV1_2' in newhdu.header): del newhdu.header['PV1_2']

        # if ('PV2_0' in newhdu.header): del newhdu.header['PV2_0']
        # if ('PV2_1' in newhdu.header): del newhdu.header['PV2_1']
        # if ('PV2_2' in newhdu.header): del newhdu.header['PV2_2']

        # print "before=",newhdu.header['CRPIX1']
        # # newhdu.header['CRPIX1'] += 100
        # # newhdu.header['CRVAL1'] += 0.2
        # print "after =",newhdu.header['CRPIX1']

        # print "\n"*5,"="*140
        # print "newhdu="
        # print newhdu.header
        # print "\n"*5,"="*140
        # print "oldhdu="
        # print imghdu.header

        wcs_out = astWCS.WCS(newhdu.header, mode='pyfits')

        # wcs_out.header['PV1_0'] = 0.
        # wcs_out.header['PV1_1'] = 1.
        # wcs_out.header['PV1_2'] = 0.
        # wcs_out.header['PV2_0'] = 0.
        # wcs_out.header['PV2_1'] = 1.
        # wcs_out.header['PV2_2'] = 0.

        # wcs_out.updateFromHeader()
        radec_new = numpy.array(wcs_out.pix2wcs(xy[:,0]-1, xy[:,1]-1))

        # Now fit the new WCS system to absorb the changed values into the 
        # CRVAL and CD matrix
        catalog = numpy.zeros((n_stars, 6))

        # Columns of catalog to be compatible with optimize_wcs_solution
        # 3/4: x/y
        # last 2 columns: reference ra/dec

        catalog[:,0:2] = radec[:,:]
        catalog[:,2:4] = xy[:,:]
        catalog[:,4:6] = radec[:,:]

        numpy.savetxt("optimize.%s.x" % (extname), catalog)

        print "imghdu -> \n",radec[:5,:]
        print "newhdu -> \n",radec_new[:5,:]

#        header_keywords  = ('CRPIX1', 'CRPIX2',
#                            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',)
#'CRVAL1', 'CRVAL2',
                            
        header_keywords  = ('CRPIX1', 'CRPIX2',
                            'CRVAL1', 'CRVAL2',
                            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                            'PV1_4', 'PV1_5', 'PV1_6', 'PV1_7', 'PV1_8', 'PV1_9',
                            'PV1_10', 'PV1_12', 'PV1_13', 'PV1_14', 'PV1_15', 'PV1_16', 
                            'PV1_17', 'PV1_18', 'PV1_19', 'PV1_20', 'PV1_21', 'PV1_22', 
                            'PV2_4', 'PV2_5', 'PV2_6', 'PV2_7', 'PV2_8', 'PV2_9',
                            'PV2_10', 'PV2_12', 'PV2_13', 'PV2_14', 'PV2_15', 'PV2_16', 
                            'PV2_17', 'PV2_18', 'PV2_19', 'PV2_20', 'PV2_21', 'PV2_22', 
        )
        header_keywords  = ('CRPIX1', 'CRPIX2',
                            'CRVAL1', 'CRVAL2',
                            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',)

        # Find the best solution.
        # The new header values are still in the wcs_out class
        #print "before xxx", wcs_out.header['CRPIX1']
        #wcs_out.header['CRPIX1'] += 100
        print "before fit", wcs_out.header['CRPIX1']
        dev_ccmatch.optimize_wcs_solution(catalog, wcs_out.header, header_keywords)
        print "after  fit", wcs_out.header['CRPIX1']

        # Activate the changes to the WCS headers
        # print "\n"*20
        # print "updateFromHeader\n"
        # print wcs_out.header
        wcs_out.updateFromHeader()
        # print "="*140
        # print wcs_out.header
        
        # compute new catalog
        radec_new = numpy.array(wcs_out.pix2wcs(xy[:,0]-1, xy[:,1]-1))
        #print radec[:5,:]

        catalog[:,0:2] = radec_new[:,:]

        catalog_out = numpy.append(radec, radec_new, axis=1)
        numpy.savetxt("optimize.%s" % (extname), catalog_out)

        dradec = radec - radec_new #catalog[:,0:2]  - catalog[:,4:6]
        print dradec[:5,:] * 3600. * 1000.

        # Update all headers with the new fitted values
        for keyname in hdr:
            # print keyname, hdr[keyname]
            if (keyname in wcs_out.header):
                hdulist[ext].header[keyname] = wcs_out.header[keyname]
                if (wcs_out.header[keyname] != wcs_in.header[keyname]):
                    print "found change in key",keyname,", was:",wcs_in.header[keyname],"  now: ",wcs_out.header[keyname]


        # hdulist[ext].data = numpy.ones((4000,4000))
        
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

