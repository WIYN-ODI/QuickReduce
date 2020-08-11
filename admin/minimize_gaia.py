#!/usr/bin/env python

import os
import sys
import astropy.io.fits as pyfits

if __name__ == "__main__":

    for fn in sys.argv[1:]:

        out_fn = fn[:-5]+".mini.fits"
        if (os.path.isfile(out_fn)):
            continue

        print("Minimizing %s to %s" % (fn, out_fn))

        hdulist = pyfits.open(fn)

        col_names = ['ra', 'dec', 'ra_error', 'dec_error', 'phot_g_mean_mag']

        flux = hdulist[1].data.field('phot_g_mean_flux')
        flux_error = hdulist[1].data.field('phot_g_mean_flux_error')
        mag_error = flux_error / flux
        mag_error[(flux<=0) | (flux_error<=0)] = -99.9

        mas_to_deg = (1/1000.) / 3600.
        columns = [
            pyfits.Column(name='ra', format='D', unit='Angle[deg]',
                          array=hdulist[1].data.field('ra')),
            pyfits.Column(name='dec', format='D', unit='Angle[deg]',
                          array=hdulist[1].data.field('dec')),

            pyfits.Column(name='ra_error', format='D', unit='Angle[deg]',
                          array=hdulist[1].data.field('ra_error')*mas_to_deg),
            pyfits.Column(name='dec_error', format='D', unqit='Angle[deg]',
                          array=hdulist[1].data.field('dec_error')*mas_to_deg),

            pyfits.Column(name='phot_g_mean_mag', format='D', unit='Magnitude[mag]',
                          array=hdulist[1].data.field('phot_g_mean_mag')),
            pyfits.Column(name='phot_g_mean_mag_error', format='D', unit='Magnitude[mag]',
                          array=mag_error),
        ]
        # for col in col_names:
        #     for tbcol in hdulist[1].columns:
        #         if (tbcol.name == col):
        #             columns.append(pyfits.Column(name=col,
        #                                          format=tbcol.format,
        #                                          disp=tbcol.disp,
        #                                          unit=tbcol.unit,
        #                                          array=hdulist[1].data.field(col)))
        #
        # columns.append(pyfits.Column(
        #     name='phot_g_mean_mag_error', array=mag_error, unit='Magnitude[mag]', format='D'))

        coldefs = pyfits.ColDefs(columns)
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
        tbhdu.name = 'GAIA_CAT'

        # assemble the output FITS file
        # copy primary extension
        out_list = [hdulist[0], tbhdu]

        out_hdulist = pyfits.HDUList(out_list)
        out_hdulist.writeto(out_fn, overwrite=True)