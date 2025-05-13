#!/usr/bin/env python3


import astropy.io.fits as pyfits
import numpy
import sys

import wiyn_filters



if __name__ == "__main__":

    # read stuff from user
    cont_fn = sys.argv[1]
    nb_fn = sys.argv[2]

    out_fn = sys.argv[3]

    # open images
    cont_hdu = pyfits.open(cont_fn)
    nb_hdu = pyfits.open(nb_fn)

    cont = cont_hdu[0].data
    nb = nb_hdu[0].data

    # get filter names
    cont_filter = cont_hdu[0].header['FILTER']
    nb_filter = nb_hdu[0].header['FILTER']

    # look up filter definitions
    (_, _, _, cont_fwhm, _, _, _, _, cont_area, _, _) = wiyn_filters.filter_bandpass[cont_filter]
    (_, _, _, nb_fwhm, _, _, _, _, nb_area, _, _) = wiyn_filters.filter_bandpass[nb_filter]

    print("Cont: ", cont_fwhm, cont_area)
    print("NB:   ", nb_fwhm, nb_area)

    # calculate continuum-band image without contribution from narrowband image
    cont_only = (cont * cont_area - nb * nb_area) / (cont_area - nb_area)

    # subtract continuum level from narrowband image
    nb_only = nb - cont_only

    # save results
    nb_hdu[0].data = nb_only
    nb_hdu.writeto(out_fn, overwrite=True)



