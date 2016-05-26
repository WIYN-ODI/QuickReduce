#!/usr/bin/env python

import numpy, scipy, pyfits
import os, sys
import scipy.ndimage

save_products = True #False

def mask_outliers_in_stack(median_frame,
                           singles):

    ref_hdu = pyfits.open(median_frame)
    
    ref_weight = median_frame[:-5]+".weight.fits"
    weight_hdu = pyfits.open(ref_weight)

    ref_smoothed = scipy.ndimage.filters.gaussian_filter(
        input=ref_hdu[0].data, 
        sigma=1, 
        order=0, 
        output=None, 
        mode='constant', 
        cval=0.0, 
        )

    ref_median = scipy.ndimage.filters.median_filter(
        input=ref_hdu[0].data, 
        size=3, 
        footprint=None, output=None, 
        mode='constant', cval=0.0)
    #
    # compute noise in median image
    #
    # ref_noise = #numpy.sqrt(
    #     (ref_hdu[0].data + ref_hdu[0].header['SKYLEVEL']) * weight_hdu[0].data)
    ref_electrons = (ref_hdu[0].data + ref_hdu[0].header['SKYLEVEL']) \
                / weight_hdu[0].data 
    ref_noise = numpy.sqrt(ref_electrons) * weight_hdu[0].data
    #* ref_hdu[0].header['GAIN']

        #ref_hdu[0].header['GAIN'])
    
    single_frame_bias = 1.3

    for fn in singles:
        hdulist = pyfits.open(fn)
        
        data = numpy.array(hdulist[0].data)
        #diff = data - ref_hdu[0].data
        diff = data - ref_smoothed

        #
        # compute noise level
        #
        fn_weight = fn[:-5]+".weight.fits"
        weight_hdu = pyfits.open(fn_weight, mode='update')
        # img_electrons = (hdulist[0].data + ref_hdu[0].header['SKYLEVEL']) \
        #             * ref_hdu[0].header['GAIN'] \
        #             / ref_hdu[0].header['NCOMBINE'] 
        img_electrons = (hdulist[0].data + ref_hdu[0].header['SKYLEVEL']) \
                    / weight_hdu[0].data * single_frame_bias
        img_noise = numpy.sqrt(img_electrons) * weight_hdu[0].data / single_frame_bias
        #                    / weight_hdu[0].data
        #

        diff_noise = numpy.hypot(img_noise, ref_noise)

        if (save_products):
            hdulist[0].data = diff
            hdulist.writeto(fn[:-5]+".diff.fits", clobber=True)
            hdulist[0].data = img_electrons
            hdulist.writeto(fn[:-5]+".electrons.fits", clobber=True)
            hdulist[0].data = img_noise
            hdulist.writeto(fn[:-5]+".imgnoise.fits", clobber=True)
            hdulist[0].data = img_noise
            hdulist.writeto(fn[:-5]+".diffnoise.fits", clobber=True)
        
        abs_diff = numpy.fabs(diff)
        #bad = abs_diff > 3*img_noise
        bad = (abs_diff > 3*diff_noise) & (abs_diff > 0.6*ref_median)
        #bad = (abs_diff > 3*diff_noise) & (abs_diff > 0.3*ref_smoothed)
        data[bad] = numpy.NaN

        if (save_products):
            hdulist[0].data = data
            hdulist.writeto(fn[:-5]+".fixed.fits", clobber=True)
        
        weight_hdu[0].data[bad] = 0.
        weight_hdu.flush()
        
        weight_hdu.close()
        hdulist.close()

    if (save_products or True):
        ref_hdu[0].data = ref_noise
        ref_hdu.writeto(median_frame[:-5]+".noise.fits", clobber=True)
        ref_hdu[0].data = ref_smoothed
        ref_hdu.writeto(median_frame[:-5]+".smoothed.fits", clobber=True)
        ref_hdu[0].data = ref_smoothed
        ref_hdu.writeto(median_frame[:-5]+".median.fits", clobber=True)
        ref_hdu[0].data = ref_electrons
        ref_hdu.writeto(median_frame[:-5]+".electrons.fits", clobber=True)


if __name__ == "__main__":

    ref_file = sys.argv[1]

    mask_outliers_in_stack(median_frame=ref_file,
                           singles=sys.argv[2:])
