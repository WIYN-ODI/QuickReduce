#!/usr/bin/env python

import pyfits
import numpy
import os, sys
import itertools
import bottleneck

sys.path.append("/work/podi_devel/")
from podi_definitions import *



if __name__ == "__main__":
    
    filename = sys.argv[1]

    hdulist = pyfits.open(filename)

    binning = hdulist[0].header['BINNING']
    skylevel = hdulist[0].header['SKYLEVEL']
    gain = hdulist[0].header['GAIN']

    print skylevel, gain
    skynoise = hdulist[0].header['SKYNOISE'] #math.sqrt(skylevel*gain)
    print skynoise

    w=20/binning
    for ext in hdulist:
        if (not is_image_extension(ext)):
            continue

        excesses = numpy.empty((8,8))
        excesses[:,:] = numpy.NaN
        # skyvalue = numpy.empty((8,8))
        # skyvalue[:,:] = numpy.NaN

        for cx, cy in itertools.product(range(8), repeat=2):

            #
            # Get pixel coord for this cell
            #
            #print binning

            x1,x2,y1,y2 = cell2ota__get_target_region(cx, cy, binning=binning, trimcell=0)
            x21 = (x2-x1)/2

            # extract the mean value in the bottom corner
            corner = bottleneck.nanmean(ext.data[y1:y1+w, x1:x1+w].astype(numpy.float32))

            # also get the value in the bottom center
            center = bottleneck.nanmean(ext.data[y1:y1+w, x1+x21-w/2:x1+x21+w//2].astype(numpy.float32))

            excess = corner - center
            #print ext.name, cx, cy, corner, center, excess
            
            excesses[cx,cy] = excess

        _mean = bottleneck.nanmean(excesses)
        _median = bottleneck.nanmedian(excesses)

        guideota = (_median > 10*skynoise)
        print ext.name, _mean, _median, skylevel, guideota


        #break
