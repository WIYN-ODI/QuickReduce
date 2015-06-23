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





def minimize_wcs_error(p, xy, ref_radec, astwcs, optimize_header_keywords):

    # Transfer all fitting parameters to astWCS
    for i in range(len(optimize_header_keywords)):
        astwcs.header[optimize_header_keywords[i]] = p[i]
    # and update astWCS so the changes take effect
    astwcs.updateFromHeader()

    # Now compute all Ra/Dec values based on the new WCS solution
    src_radec = numpy.array(astwcs.pix2wcs(xy[:,0], xy[:,1]))

    # This gives us the Ra/Dec values as 2-d array
    # compute difference from the Ra/Dec of the reference system
    src_ref = src_radec - ref_radec

    # return the 1-d version for optimization
    return src_ref.ravel()




if __name__ == "__main__":


    input_wcs = sys.argv[1]

    output_wcs = sys.argv[2]

    hdulist = pyfits.open(input_wcs)
    out_hdulist = pyfits.open(input_wcs)
    
    d_crval1, d_crval2 = None, None

    dummyfile = open("simplify.debug", "w")
    n_stars = 5000

    wcs_offset = numpy.array([10.0, 0.0])


    #
    # Now compute the center position of the 3/3-4/4 OTAs to get the new 
    # reference position 
    #
    radec_ref = numpy.empty((4,2))
    radec_ref[:,:] = numpy.NaN

    corners = {
        'OTA33.SCI': (4096, 4096),
        'OTA34.SCI': (4096,    0),
        'OTA43.SCI': (   0, 4096),
        'OTA44.SCI': (   0,    0),
    }

    for idx, ota in enumerate(corners):
        try:
            ext = hdulist[ota]
            wcs = astWCS.WCS(ext.header, mode='pyfits')
            wcs.header['NAXIS'] = 2
            wcs.header['NAXIS1'] = 4096
            wcs.header['NAXIS2'] = 4096
            wcs.updateFromHeader()
            x,y = corners[ota]
            print numpy.array(wcs.pix2wcs(x,y))
            radec_ref[idx,:] = numpy.array(wcs.pix2wcs(x,y))
        except:
            # ignore OTAs that are missing
            # this implies we need at least ONE of the 4 central OTAs to 
            # make this correction
            podi_logging.log_exception()
            pass

    ref_point = bottleneck.nanmean(radec_ref, axis=0)
    print "REF:",ref_point

    #ref_point -= wcs_offset

    #
    # Now go through each OTA and correct the CRPIX values accordingly 
    # to match the new CRVAL values we just computed
    #
    for idx, ext in enumerate(out_hdulist):
        if (not type(ext)== pyfits.hdu.image.ImageHDU):
            continue

        # Create WCS system for this OTA
        hdr = pyfits.Header(ext.header)
        hdr['NAXIS'] = 2
        hdr['NAXIS1'] = 4096
        hdr['NAXIS2'] = 4096
        wcs = astWCS.WCS(hdr, mode='pyfits')
        
        # compute the pixel position of the reference Ra/Dec point
        #
        #print ref_point[0], ref_point[1]
        print ext.name, wcs.wcs2pix(ref_point[0], ref_point[1])

#    sys.exit(0)



    for ext in range(1, len(hdulist)):
        if (not type(hdulist[ext])== pyfits.hdu.image.ImageHDU):
            continue

        extname = hdulist[ext].name

        # Read input WCS
        hdulist[ext].header['NAXIS'] = 2
        hdulist[ext].header['NAXIS1'] = 4096
        hdulist[ext].header['NAXIS2'] = 4096
        hdulist[ext].header['CRVAL1'] += wcs_offset[0]
        hdulist[ext].header['CRVAL1'] += wcs_offset[1]

        in_wcs = astWCS.WCS(hdulist[ext].header, mode='pyfits')

        # generate random coordinates
        xy = numpy.random.random((n_stars,2))*4096

        print xy[:10]

        # convert to Ra/Dec
        radec = numpy.array(in_wcs.pix2wcs(xy[:,0], xy[:,1]))
        print radec[:10]
        
#        continue

        #break

        print "------------------\n"*5
        # continue

        # now change the initial PV values to 1/0, and prepare re-fitting
        out_wcs = astWCS.WCS(hdulist[ext].header, mode='pyfits')
        out_wcs.header['PV1_0'] = 0.0
        out_wcs.header['PV1_1'] = 1.0
        out_wcs.header['PV2_0'] = 0.0
        out_wcs.header['PV2_1'] = 1.0
        out_wcs.updateFromHeader()

        max_pv = 0
        for i in range(0,50):
            if "PV1_%d" % (i) in out_wcs.header:
                max_pv = i

        # Load initial values
        n_parameters = 2*max_pv + 4
        p_start = numpy.zeros(n_parameters)

        header_names = ['CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
#        header_names = ['CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
        for i in range(2,max_pv):
            header_names.append('PV1_%d' % i)
            header_names.append('PV2_%d' % i)

        for idx, kw in enumerate(header_names):
            p_start[idx] = out_wcs.header[kw] if kw in out_wcs.header else 0.0
        # p_start[0] = out_wcs.header['CD1_1']
        # p_start[1] = out_wcs.header['CD1_2']
        # p_start[2] = out_wcs.header['CD2_1']
        # p_start[3] = out_wcs.header['CD2_2']

        # for i in range(2,max_pv):
        #     kw = 'PV1_%d' % i
        #     p_start[2*i]   = out_wcs.header[kw] if kw in out_wcs.header else 0.0
        #     kw = 'PV2_%d' % i
        #     p_start[2*i+1] = out_wcs.header[kw] if kw in out_wcs.header else 0.0
        #     if kw in out_wcs.header:
        #         print kw, out_wcs.header[kw]

        print max_pv, p_start

        # Now prepare for re-fitting the WCS with the updated PV values
        fit_args = (xy, radec, out_wcs, header_names)
        try:
            fit = scipy.optimize.leastsq(minimize_wcs_error,
                                         p_start, 
                                         args=fit_args,
                                         #maxfev=1000,
                                         full_output=1)
            #print fit
            p_final = fit[0]
        except:
            print "fit failed:",extname
            p_final = p_start

        #print xy.shape
        #print xy[:10]


        for idx, kw in enumerate(header_names):
            out_wcs.header[kw] = p_final[idx]
        out_wcs.updateFromHeader()
        out_radec = numpy.array(in_wcs.pix2wcs(xy[:,0], xy[:,1]))

        numpy.savetxt("XY_%s" % extname, xy)
        numpy.savetxt("RADEC_IN_%s" % extname, radec)
        numpy.savetxt("RADEC_OUT_%s" % extname, out_radec)

        out_hdulist[ext].header['PV1_0'] = 0.0
        out_hdulist[ext].header['PV1_1'] = 1.0
        out_hdulist[ext].header['PV2_0'] = 0.0
        out_hdulist[ext].header['PV2_1'] = 1.0
        for idx, kw in enumerate(header_names):
            out_hdulist[ext].header[kw] = p_final[idx]

    out_hdulist.writeto(output_wcs, clobber=True)


    sys.exit(0)

