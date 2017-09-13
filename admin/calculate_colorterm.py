#!/usr/bin/env python

import os
import sys
import pyfits
import numpy

import scipy.odr
import argparse

if __name__ == "__main__":

    #
    # Handle all command line stuff
    #
    parser = argparse.ArgumentParser(
        description='Create light-curve database for a given SExtractor configuration.')
    parser.add_argument(
        'filter_1', type=str, #nargs=1,
        metavar='filter_1',
        help='shsorter-wavelength filter')
    parser.add_argument(
        'filter_2', metavar='filter_2', type=str, #nargs=1,
        help='redder wavelength filter')
    parser.add_argument(
        'filenames', metavar='sex.param', type=str, nargs="+",
        help='input files')
    parser.add_argument('--order', dest='poly_order',
                        default=1, type=int, help='polygon order')
    parser.add_argument('-o', '--output', dest='output',
                        default="dummy.out", type=str, help='output filename')
    parser.add_argument('--all', dest='all_stars',
                        default=False, action='store_true', help='use all stars')
    parser.add_argument('--maxerr', dest='max_error',
                        default=0.05, type=float, help='use all stars')
    # parser.add_argument('--all', dest='all_stars',
    #                     default=False, action='store_true', help='use all stars')
    args = parser.parse_args()

    filter1 = args.filter_1 #sys.argv[1]
    filter2 = args.filter_2 #sys.argv[2]


    full_catalog = None


    for i_fn, fn in enumerate(args.filenames): #sys.argv[3:]):
        # print fn

        hdulist = pyfits.open(fn)
        try:
            photcalib_tbhdu = hdulist['CAT.PHOTREF']
        except:
            try:
                photcalib_tbhdu = hdulist['CAT.PHOTCALIB']
            except:
                print "No PHOTCAL data found, aborting"
                continue


        mag0size = hdulist[0].header['MAG0SIZE']

        _odi_col = 'ODI_MAG_D%02d' % (mag0size*10)
        _odi_err = 'ODI_ERR_D%02d' % (mag0size*10)
        photfilt = hdulist[0].header['PHOTFILT'].upper()

        used = photcalib_tbhdu.data.field('USED4CALIB')
        print "Read %5d stars from %s" % (used.shape[0], fn)
        # print used.shape

        odi_mag = photcalib_tbhdu.data.field(_odi_col) #- 2.5*numpy.log10(hdulist[0].header['EXPMEAS'])
        odi_err = photcalib_tbhdu.data.field(_odi_err)

        ref_mag = photcalib_tbhdu.data.field("REF_%s" % (photfilt))
        ref_err = photcalib_tbhdu.data.field("REF_ERR_%s" % (photfilt))

        ref_color = photcalib_tbhdu.data.field("REF_%s" % (filter1)) - \
                    photcalib_tbhdu.data.field("REF_%s" % (filter2))
        ref_color_error = numpy.hypot(
            photcalib_tbhdu.data.field("REF_ERR_%s" % (filter1)),
            photcalib_tbhdu.data.field("REF_ERR_%s" % (filter2)))

        rel_zp = ref_mag - odi_mag - hdulist[0].header['PHOTZP_X']
        rel_zp_err = numpy.hypot(odi_err, ref_err)


        catalog = numpy.array([
            ref_color, ref_color_error, rel_zp, rel_zp_err
        ]).T
        # print catalog.shape

        if (args.all_stars):
            select = (ref_color_error < args.max_error) & (rel_zp_err < args.max_error)
            catalog = catalog[select]
            print "  down-selecting to %5d stars with errors < %.3f mag" % (
                catalog.shape[0], args.max_error
            )
        else:
            catalog = catalog[used]
            print "  down-selecting to %5d stars that were used for photometric calibration" % (
                catalog.shape[0]
            )


        # print catalog.shape

        full_catalog = catalog if full_catalog is None \
            else numpy.append(full_catalog, catalog, axis=0)


    #
    # Now we have the full catalog
    #
    print "\nSaving combined catalog to %s" % (args.output)
    numpy.savetxt(args.output, full_catalog)

    print "\nStarting fit using polygon of order %d for total of %d stars" % (
        args.poly_order, full_catalog.shape[0]
    )
    polynom = scipy.odr.polynomial(args.poly_order)

    data = scipy.odr.RealData(
        x=full_catalog[:, 0],
        y=full_catalog[:, 2],
        sx = full_catalog[:,1],
        sy = full_catalog[:,3]
    )

    myodr = scipy.odr.ODR(data, polynom)
    fit = myodr.run()
    fit.pprint()
    # print fit
    print fit.beta

    #
    # compare quality of photometric calibration with and without color
    # correction
    #
    before_correction = numpy.percentile(full_catalog[:,2], [16,84,50])
    before_median = before_correction[2]
    before_sigma = (before_correction[1]-before_correction[0])/2.
    print "BEFORE: median=%.4f sigma=%f" % (
        before_median, before_sigma)

    corrected_catalog = full_catalog.copy()
    correction = numpy.polyval(fit.beta[::-1], corrected_catalog[:,0])
    corrected_catalog[:,2] -= correction

    print "\nSaving corrected catalog to %s" % (args.output+".corr")
    numpy.savetxt(args.output+".corr", corrected_catalog)

    after_correction = numpy.percentile(corrected_catalog[:,2], [16,84,50])
    after_median = after_correction[2]
    after_sigma = (after_correction[1]-after_correction[0])/2.
    print "after: median=%.4f sigma=%f" % (
        after_median, after_sigma)
