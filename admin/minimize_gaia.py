#!/usr/bin/env python

import os
import sys
import pyfits

if __name__ == "__main__":

    for fn in sys.argv[1:]:

        out_fn = fn[:-5]+".mini.fits"
        if (os.path.isfile(out_fn)):
            continue

        print "Minimizing %s to %s" % (fn, out_fn)

        hdulist = pyfits.open(fn)

        col_names = ['ra', 'dec', 'ra_error', 'dec_error', 'phot_g_mean_mag']
        columns = []
        for col in col_names:
            for tbcol in hdulist[1].columns:
                if (tbcol.name == col):
                    columns.append(pyfits.Column(name=col,
                                                 format=tbcol.format,
                                                 disp=tbcol.disp,
                                                 unit=tbcol.unit,
                                                 array=hdulist[1].data.field(col)))

        flux = hdulist[1].data.field('phot_g_mean_flux')
        flux_error = hdulist[1].data.field('phot_g_mean_flux_error')
        mag_error = flux_error / flux
        mag_error[(flux<=0) | (flux_error<=0)] = -99.9
        columns.append(pyfits.Column(
            name='phot_g_mean_mag_error', array=mag_error, unit='Magnitude[mag]', format='D'))

        coldefs = pyfits.ColDefs(columns)
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'GAIA_CAT'

        # assemble the output FITS file
        # copy primary extension
        out_list = [hdulist[0], tbhdu]

        out_hdulist = pyfits.HDUList(out_list)
        out_hdulist.writeto(out_fn, clobber=True)