#!/usr/bin/env python

import os
import sys
import numpy
from podi_definitions import is_guide_ota, is_image_extension
import pyfits
from podi_commandline import *


if __name__ == "__main__":

    debug = cmdline_arg_isset("-debug")
    print(get_clean_cmdline())

    for fn in get_clean_cmdline()[1:]:

        guide_otas = []

        print(fn)
        hdulist  = pyfits.open(fn)
        for ext in hdulist:
            if (not is_image_extension(ext)):
                continue
            # if (ext.name != 'OTA41.SCI'):
            #     continue

            # is_guide, excesses, _mean, _median, skynoise, corners, centers =  is_guide_ota(hdulist[0], ext, debug=True)
            # pyfits.HDUList(centers).writeto("guideota_centers.fits", clobber=True)
            # pyfits.HDUList(corners).writeto("guideota_corners.fits", clobber=True)
            #
            # print(excesses)
            # print(skynoise)
            # print(_mean, _median)

            is_guide = is_guide_ota(hdulist[0], ext)
            # print(ext.name, is_guide)
            if (is_guide):
                guide_otas.append(ext.name)

        print("%s ==> %s" % (fn, ",".join(guide_otas)))
