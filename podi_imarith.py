#! /usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#


"""

Provide basic support for image manipulations. 

The syntax is largely based on the IRAF task imarith:

    ``podi_imarith.py frame1.fits (op) frame2.fits output.fits

Note that instead of frame2.fits you can also specify a numeric value.

Operations currently supported are:
* **+** addition
* **-** subtraction
* **x** multiplication
* **/** division
* **^** raise to the n-th power. This only works with numerical values as 
  second parameter.

"""

import sys
import os
import pyfits
import numpy
import scipy
import scipy.stats

from podi_definitions import *
from podi_commandline import *


def hdu_imarith(hdu1, op, hdu2, simple=False, numeric_2=None, rebin_fac=1.):

    logger = logging.getLogger("ImArith")

    out_hdulist = [pyfits.PrimaryHDU(header=hdu1[0].header)]

    # Now go though each extension and perform the operation
    for idx_img1 in range(0, len(hdu1)):
        if (not is_image_extension(hdu1[idx_img1])):
            continue

        img1 = pyfits.ImageHDU(data=hdu1[idx_img1].data,
                               header=hdu1[idx_img1].header,
                               name=hdu1[idx_img1].name)
        img1.data = img1.data.astype(numpy.float)

        fppos1 = img1.header['EXTNAME'] if ('EXTNAME' in img1.header and not simple) else idx_img1
        # stdout_write("\rComputing extension %s (%2d of %2d) ..." % (str(fppos1), idx_img1, len(hdu1) - 1))
        if (numeric_2 is not None):
            if (op == "+"):
                img1.data += numeric_2
            elif (op == "-"):
                img1.data -= numeric_2
            elif (op == "/"):
                img1.data /= numeric_2
            elif (op == "x"):
                img1.data *= numeric_2
            elif (op == "^"):
                img1.data = numpy.pow(img1.data, numeric_2)
            else:
                logger.warning("Unknown operation %s" % (op))

        else:
            for img2_idx in range(len(hdu2)):  # [0:]:
                fppos2 = hdu2[img2_idx].header['EXTNAME'] if (
                'EXTNAME' in hdu2[img2_idx].header and not simple) else img2_idx
                img2 = hdu2[fppos2]
                if (fppos2 == fppos1):
                    # This is the one
                    img2.data = img2.data.astype(numpy.float)

                    if (op == "+"):
                        img1.data += img2.data

                    elif (op == "-"):
                        img1.data -= img2.data

                    elif (op == "/"):
                        img1.data /= img2.data

                    elif (op == "x"):
                        img1.data *= img2.data

                    else:
                        logger.warning("Unknown operation %s" % (op))

        if (rebin_fac > 1):
            img1.data = rebin_image(img1.data, rebin_fac)

        out_hdulist.append(img1)

    return pyfits.HDUList(out_hdulist)

    # Now all the work is done, all final data is stored in img2, write results to new file
    #if (output == None):
    #    return hdu1



def imarith(input_1, op, input_2, output, simple):

    stdout_write("\nOpening input files ...")
    # Open both input fits files
    hdu_1 = pyfits.open(input_1)

    numeric_2 = None
    hdu_2 = None
    if (not os.path.isfile(input_2)):
        numeric_2 = float(input_2)
        print(numeric_2)
    else:
        hdu_2 = pyfits.open(input_2)

    stdout_write(" done!\n")

    rebin_fac = int(cmdline_arg_set_or_default("-bin", 1))
    
    output_hdu = hdu_imarith(hdu_1, op, hdu_2,
                         simple=simple,
                         numeric_2=numeric_2,
                         rebin_fac=rebin_fac)

    stdout_write(" writing output ...")
    clobberfile(output)
    output_hdu.writeto(output, output_verify='fix+ignore')
    stdout_write(" done!\n\n")

    return

if __name__ == "__main__":

    # Read in the input parameters
    input_1 = sys.argv[1]
    op = sys.argv[2]
    input_2 = sys.argv[3]
    output = sys.argv[4]

    try:
        simple = sys.argv[5] == "simple"
    except:
        simple=False

    imarith(input_1, op, input_2, output, simple)
