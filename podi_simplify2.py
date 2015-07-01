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

from podi_commandline import *




def minimize_wcs_error(p, xy, ref_radec, astwcs, optimize_header_keywords):

    # Transfer all fitting parameters to astWCS
    for i in range(len(optimize_header_keywords)):
        astwcs.header[optimize_header_keywords[i]] = p[i]
    # and update astWCS so the changes take effect
    astwcs.updateFromHeader()
#    print astwcs.header

    # Now compute all Ra/Dec values based on the new WCS solution
    src_radec = numpy.array(astwcs.pix2wcs(xy[:,0], xy[:,1]))

    # This gives us the Ra/Dec values as 2-d array
    # compute difference from the Ra/Dec of the reference system
    src_ref = src_radec - ref_radec

    # return the 1-d version for optimization
    return src_ref.ravel()




if __name__ == "__main__":

    zero_crval = True

    if (cmdline_arg_isset("-merge")):
        recompute_distortions = False
        zero_crval = cmdline_arg_isset("-zero")

        hdus = [pyfits.PrimaryHDU()]
        extname_list = []

        for fn in get_clean_cmdline()[1:-1]:
            print "Adding", fn
            hdulist = pyfits.open(fn)
            for ext in hdulist:
                if (not type(ext)== pyfits.hdu.image.ImageHDU):
                    continue
                extname = ext.name
                if (extname in extname_list):
                    # this is a duplicate OTA
                    pass
                else:
                    if (not cmdline_arg_isset("-keepdata")):
                        ext.data = None
                    hdus.append(ext)
                    extname_list.append(extname)
        hdulist = pyfits.HDUList(hdus)
        out_hdulist = pyfits.HDUList(hdus)
 
        output_wcs = get_clean_cmdline()[-1]
    else:
        input_wcs = get_clean_cmdline()[1]
        output_wcs = get_clean_cmdline()[2]
        hdulist = pyfits.open(input_wcs)
        out_hdulist = pyfits.open(input_wcs)
        recompute_distortions = True #False

    
    d_crval1, d_crval2 = None, None

    dummyfile = open("simplify.debug", "w")
    n_stars = 500

    wcs_offset = numpy.array([10.0, 0.0])

    #
    # Make sure to have reasonable CRVALs
    # 
    for ext in hdulist:
        if (not type(ext)== pyfits.hdu.image.ImageHDU):
            continue
        if ('CRVAL1' in ext.header):
            crval = ext.header['CRVAL1']
            ext.header['CRVAL1'] = math.fmod(crval-math.floor(crval/360.0)*360,360.0)

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
            hdr = hdulist[ota].header.copy()
            hdr['NAXIS'] = 2
            hdr['NAXIS1'] = 4096
            hdr['NAXIS2'] = 4096
            hdr['CRVAL1'] = math.fmod(hdr['CRVAL1'] + wcs_offset[0], 360.0)
            hdr['CRVAL2'] += wcs_offset[1]
            wcs = astWCS.WCS(hdr, mode='pyfits')
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
    # ref_point += wcs_offset

    print "REF:",ref_point

    crval = numpy.array([hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']])

    d_crval = ref_point - crval
    print "current CRVAL:", hdulist[1].header['CRVAL1'], hdulist[1].header['CRVAL2']

    #
    #
    # Now change the WCS solution so as to put the CRPIX of all OTAs to match
    # up with this center of the focal plane
    #
    #
    


    # #
    # # Now go through each OTA and correct the CRPIX values accordingly 
    # # to match the new CRVAL values we just computed
    # #
    # for idx, ext in enumerate(out_hdulist):
    #     if (not type(ext)== pyfits.hdu.image.ImageHDU):
    #         continue

    #     # Create WCS system for this OTA
    #     hdr = pyfits.Header(ext.header)
    #     hdr['NAXIS'] = 2
    #     hdr['NAXIS1'] = 4096
    #     hdr['NAXIS2'] = 4096
    #     wcs = astWCS.WCS(hdr, mode='pyfits')

    #     #
    #     # compute the pixel position of the reference Ra/Dec point
    #     #
    #     #print ref_point[0], ref_point[1]
    #     print ext.name, wcs.wcs2pix(ref_point[0], ref_point[1])


    #recompute_distortions = True #False
    if (recompute_distortions):

        for ext in range(len(hdulist)):
            if (not type(hdulist[ext])== pyfits.hdu.image.ImageHDU):
                continue
            # if (not hdulist[ext].name in ["OTA44.SCI", "OTA33.SCI",
            #                               "OTA34.SCI", "OTA43.SCI"]):
            #    continue

            extname = hdulist[ext].name
            print extname

            # Read input WCS
            hdr = hdulist[ext].header.copy()
            hdr['NAXIS'] = 2
            hdr['NAXIS1'] = 4096
            hdr['NAXIS2'] = 4096
            hdr['CRVAL1'] = math.fmod(hdr['CRVAL1'] + wcs_offset[0], 360.0)
            hdr['CRVAL2'] += wcs_offset[1]

            print "input:", hdr['CRVAL1'], hdr['CRVAL2']
            # hdulist[ext].header['NAXIS'] = 2
            # hdulist[ext].header['NAXIS1'] = 4096
            # hdulist[ext].header['NAXIS2'] = 4096
            # hdulist[ext].header['CRVAL1'] += wcs_offset[0]
            # hdulist[ext].header['CRVAL2'] += wcs_offset[1]


            #
            # generate random coordinates
            #
            in_wcs = astWCS.WCS(hdr, mode='pyfits')
            xy = numpy.random.random((n_stars,2))*4096

            #
            # convert to Ra/Dec
            #
            radec = numpy.array(in_wcs.pix2wcs(xy[:,0], xy[:,1]))
            #print radec[:10]


            #
            # compute the pixel position of the reference Ra/Dec point
            #
            ref_crpix = in_wcs.wcs2pix(ref_point[0], ref_point[1])
            print hdulist[ext].name, ref_point, ref_crpix
            # Set these coordinates to align with the reference point 
            # in_wcs.header['CRPIX1'] = ref_crpix[0]
            # in_wcs.header['CRPIX2'] = ref_crpix[1]

            # now change the initial PV values to 1/0, and prepare re-fitting
            out_hdr = hdr.copy()
            out_hdr['PV1_0'] = 0.0
            out_hdr['PV1_1'] = 1.0
            out_hdr['PV2_0'] = 0.0
            out_hdr['PV2_1'] = 1.0
            out_hdr['CRPIX1'] = ref_crpix[0]
            out_hdr['CRPIX2'] = ref_crpix[1]
            print "output pre-fit:", out_hdr['CRVAL1'], out_hdr['CRVAL2']
            out_wcs = astWCS.WCS(out_hdr, mode='pyfits')

    #        out_hdr['CRVAL1'] -= d_crval[0]
    #        out_hdr['CRVAL2'] -= d_crval[1]
            
            # out_wcs = astWCS.WCS(hdulist[ext].header, mode='pyfits')
    #         out_wcs.header['PV1_0'] = 0.0
    #         out_wcs.header['PV1_1'] = 1.0
    #         out_wcs.header['PV2_0'] = 0.0
    #         out_wcs.header['PV2_1'] = 1.0

    #         out_wcs.header['CRPIX1'] = ref_crpix[0]
    #         out_wcs.header['CRPIX2'] = ref_crpix[1]

    # #        out_wcs.header['CRVAL1'] -= d_crval[0]
    # #        out_wcs.header['CRVAL2'] -= d_crval[1]

            out_wcs.updateFromHeader()

            max_pv = 0
            for i in range(0,50):
                if (("PV1_%d" % (i) in out_wcs.header) or
                    ("PV2_%d" % (i) in out_wcs.header)):
                    max_pv = i

            #
            # Load initial values
            #
            # header_names = ['CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
            header_names = ['CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
            for i in range(2,max_pv):
                header_names.append('PV1_%d' % i)
                header_names.append('PV2_%d' % i)

            n_parameters = len(header_names)
            p_start = numpy.zeros(n_parameters)

            for idx, kw in enumerate(header_names):
                p_start[idx] = out_wcs.header[kw] if kw in out_wcs.header else 0.0

            # Now prepare for re-fitting the WCS with the updated PV values
            fit_args = (xy, radec, out_wcs, header_names)
            try:
                fit = scipy.optimize.leastsq(minimize_wcs_error,
                                             p_start, 
                                             args=fit_args,
                                             maxfev=500,
                                             full_output=1)
                #print fit
                p_final = fit[0]
            except:
                print "fit failed:",extname
                p_final = p_start

            #print xy.shape
            #print xy[:10]


            #
            # Check the resulting WCS solution - this is mostly for debugging
            #
            for idx, kw in enumerate(header_names):
                out_wcs.header[kw] = p_final[idx]
            out_wcs.updateFromHeader()
            out_radec = numpy.array(out_wcs.pix2wcs(xy[:,0], xy[:,1]))
            numpy.savetxt("XY_%s" % extname, xy)
            numpy.savetxt("RADEC_IN_%s" % extname, radec)
            numpy.savetxt("RADEC_OUT_%s" % extname, out_radec)

            for idx, kw in enumerate(header_names):
                out_hdulist[ext].header[kw] = p_final[idx]

            out_hdulist[ext].header['CRPIX1'] = ref_crpix[0]
            out_hdulist[ext].header['CRPIX2'] = ref_crpix[1]

            out_hdulist[ext].header['PV1_0'] = 0.0
            out_hdulist[ext].header['PV1_1'] = 1.0
            out_hdulist[ext].header['PV2_0'] = 0.0
            out_hdulist[ext].header['PV2_1'] = 1.0

            # subtract the wcs offset, since CRVAL1/2 was a free parameter
            out_hdulist[ext].header['CRVAL1'] -= wcs_offset[0]
            out_hdulist[ext].header['CRVAL2'] -= wcs_offset[1]

    #
    # Now that all OTAs have re-computed WCS distortions, re-set the CRVAL 
    # values to have the center of the WCS system be close to the center of 
    # the focal plane
    # 
    radec_ref[:,:] = numpy.NaN
    for idx, ota in enumerate(corners):
        try:
            hdr = out_hdulist[ota].header.copy()
            hdr['NAXIS'] = 2
            hdr['NAXIS1'] = 4096
            hdr['NAXIS2'] = 4096
            hdr['CRVAL1'] = math.fmod(hdr['CRVAL1'] + wcs_offset[0], 360.0)
            hdr['CRVAL2'] += wcs_offset[1]
            wcs = astWCS.WCS(hdr, mode='pyfits')
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

    print radec_ref
    ref_point = bottleneck.nanmean(radec_ref, axis=0) - wcs_offset

    print ref_point 

    if (zero_crval):
        # Now go through each OTA and subtract the reference position
        for idx, ext in enumerate(out_hdulist):
            if (not type(ext)== pyfits.hdu.image.ImageHDU):
                continue
            crval = out_hdulist[idx].header['CRVAL1'] - ref_point[0]
            out_hdulist[idx].header['CRVAL1'] = math.fmod((crval - math.floor(crval/360.)/360), 360.0)
            out_hdulist[idx].header['CRVAL2'] -= ref_point[1]

    # Finally, write the new WCS file to disk
    out_hdulist.writeto(output_wcs, clobber=True)

    sys.exit(0)

