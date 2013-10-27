#!/usr/bin/env python


import sys
import numpy
import os
from podi_definitions import *
import pyfits
import datetime
import scipy
import scipy.stats
import scipy.ndimage

import matplotlib.pyplot as plot

filename = sys.argv[1]
hdulist = pyfits.open(filename)

img = hdulist[1].data
data = rebin_image(hdulist[1].data, 4)

img = plot.imshow(data, interpolation='nearest', origin='lower')
plot.show(block=False)
plot.close()

mask = numpy.zeros(shape=data.shape)
mask[numpy.isnan(data)] = 1

data[numpy.isnan(data)] = 0


rotated = scipy.ndimage.interpolation.rotate(input=data, angle=30, axes=(1, 0), reshape=False, 
                                             mode='constant', cval=0, )

rotated_mask = scipy.ndimage.interpolation.rotate(input=mask, angle=30, axes=(1, 0), reshape=False, 
                                             mode='constant', cval=1, )

img = plot.imshow(rotated, interpolation='nearest', origin='lower', vmin=0., vmax=0.15)
plot.show(block=False)

pyfits.PrimaryHDU(data=rotated).writeto("rotated.fits", clobber=True)
pyfits.PrimaryHDU(data=rotated_mask).writeto("rotated_mask.fits", clobber=True)

rotated[rotated_mask > 0.1] = numpy.NaN
pyfits.PrimaryHDU(data=rotated).writeto("rotated_masked.fits", clobber=True)
