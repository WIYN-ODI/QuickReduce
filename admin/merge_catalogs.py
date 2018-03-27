#!/usr/bin/env python

import os
import sys
import pyfits
import numpy

if __name__ == "__main__":

    output_fn = sys.argv[1]
    cat_name = sys.argv[2]

    full_catalog = None

    ttypes = []
    tunits = []
    tforms = []
    tdisps = []

    first_file = True

    for i_fn, fn in enumerate(sys.argv[3:]):
        # print fn

        hdulist = pyfits.open(fn)
        photcalib_tbhdu = hdulist[cat_name]

        #
        # Now extract all sources for each OTA
        #
        tfields = photcalib_tbhdu.header['TFIELDS']
        nstars = photcalib_tbhdu.header['NAXIS2']
        catalog = numpy.empty((nstars, tfields))

        for field in range(tfields):
            catalog[:, field] = photcalib_tbhdu.data.field(field)

            hdr = photcalib_tbhdu.header
            if (first_file):
                _key = 'TTYPE%d' % (field + 1)
                ttypes.append(hdr[_key] if _key in hdr else None)
                _key = 'TFORM%d' % (field + 1)
                tforms.append(hdr[_key] if _key in hdr else None)
                _key = 'TDISP%d' % (field + 1)
                tdisps.append(hdr[_key] if _key in hdr else None)
                _key = 'TUNIT%d' % (field + 1)
                tunits.append(hdr[_key] if _key in hdr else None)


        # add some info to each catalog
        addtl = numpy.ones((catalog.shape[0], 3)) * \
            [i_fn+1, hdulist[0].header['PHOTZP_X'], hdulist[0].header['EXPMEAS']]
        catalog = numpy.append(catalog, addtl, axis=1)
        ttypes.extend(['FRAME_ID', 'PHOTZP', 'EXPMEAS'])
        tdisps.extend(['I3.3', 'F8.4', 'F8.4'])
        tforms.extend(['I', 'E', 'E'])
        tunits.extend(['', 'mag', 'second'])

        full_catalog = catalog if full_catalog is None else \
            numpy.append(full_catalog, catalog, axis=0)
        print(catalog.shape, full_catalog.shape)

        print("read %d entries from %s" % (nstars, fn))
        first_file = False

    #
    # Now create the combined table
    #
    print(len(ttypes), len(tforms), len(tdisps), len(tunits))
    print(full_catalog.shape[0])

    columns = []
    for col in range(full_catalog.shape[1]):

        columns.append(
            pyfits.Column(name=ttypes[col],
                      format=tforms[col],
                      unit=tunits[col],
                      disp=tdisps[col],
                      array=full_catalog[:, col])
        )

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = cat_name

    output_hdu = pyfits.HDUList([
        pyfits.PrimaryHDU(),
        tbhdu,
    ])
    output_hdu.writeto(output_fn, clobber=True)
