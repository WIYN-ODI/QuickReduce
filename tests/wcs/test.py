#!/usr/bin/env python

import os, sys, numpy, pyfits, time
import astLib.astWCS as astWCS


cat = numpy.loadtxt(sys.argv[1])

hdr = pyfits.Header().fromtextfile(sys.argv[2])

astwcs = astWCS.WCS(hdr, mode='pyfits')
xy  = cat[:, 2:4] - [1., 1.]
radec = numpy.array(astwcs.pix2wcs(xy[:, 0], xy[:, 1]))
print radec[0]



