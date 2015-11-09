#!/usr/bin/env python

import sys, pyfits
import podi_makefake56

hdu = pyfits.open(sys.argv[1])


hdu.info()

hduout = [hdu[0]]

for (old,new) in podi_makefake56.old_new:

    old_extname = "OTA%02d.SCI" % (old)
    new_extname = "OTA%02d.SCI" % (new)

    ext = pyfits.ImageHDU(header=hdu[old_extname].header)
    
    ext.name = new_extname
    ext.header['EXTNAME'] = new_extname

    hduout.append(ext)

print "writing output"

hdux = pyfits.HDUList(hduout)
hdux.info()
hdux.writeto(sys.argv[2], clobber=True)
    
