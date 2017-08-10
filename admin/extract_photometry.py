#!/usr/bin/env python

import os
import sys
import pyfits
import numpy


if __name__ == "__main__":

    input_fn = sys.argv[1]
    output_fn = sys.argv[2]

    hdulist = pyfits.open(input_fn)
    photcalib_tbhdu = hdulist['CAT.ODI']

    magzero = hdulist[0].header['PHOTZP_X']
    photzp_sd = hdulist[0].header['PHOTZPSD']

    #
    # Now extract all sources for each OTA
    #
    tfields = photcalib_tbhdu.header['TFIELDS']

    columns = []
    hdr = photcalib_tbhdu.header
    for field in range(tfields):

        _key = 'TTYPE%d' % (field + 1)
        ttype = hdr[_key] if _key in hdr else None
        _key = 'TFORM%d' % (field + 1)
        tform = hdr[_key] if _key in hdr else None
        _key = 'TDISP%d' % (field + 1)
        tdisp = hdr[_key] if _key in hdr else None
        _key = 'TUNIT%d' % (field + 1)
        tunit = hdr[_key] if _key in hdr else None

        data = photcalib_tbhdu.data.field(field)

        if (ttype is not None and
            ttype.startswith('MAG_')):

            # this is a magnitude column, so add the mag. zeropoint to calibrate
            # photometry

            data += magzero

        columns.append(
            pyfits.Column(name=ttype,
                          format=tform,
                          unit=tunit,
                          disp=tdisp,
                          array=data)
        )

        if (ttype.startswith('MAGERR_')):
            data = numpy.hypot(data, photzp_sd)
            columns.append(
                pyfits.Column(name=ttype+"_COMB",
                              format=tform,
                              unit=tunit,
                              disp=tdisp,
                              array=data)
            )

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = photcalib_tbhdu.name

    output_hdu = pyfits.HDUList([
        pyfits.PrimaryHDU(header=hdulist[0].header),
        tbhdu,
    ])
    output_hdu.writeto(output_fn, clobber=True)
    print "done calibrating %s --> %s" % (input_fn, output_fn)