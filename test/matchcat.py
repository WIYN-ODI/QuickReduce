#!/usr/local/bin/python
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

import sys
import os
import pyfits
import numpy
import scipy
import scipy.optimize


import matplotlib.pyplot as pl
import podi_fixwcs

# Create reference catalog
reference = numpy.random.sample((50,2))
print reference.shape


starcat = numpy.zeros((20,2))

random_shift = numpy.random.sample((1,2)) * 0.15
print "Random shift=", random_shift

for i in range(starcat.shape[0]):
    star = numpy.random.randint(0, reference.shape[0]-1)

    starcat[i,:] = reference[star,:] - random_shift[0,:]


scatter_max = 0.5 / 3600.
scattered_starcat = starcat + numpy.random.randn(starcat.shape[0], starcat.shape[1])*scatter_max

#pl.plot(reference[:,0], reference[:,1], "ro",
#        starcat[:,0], starcat[:,1], 'b^')
#pl.title("Ref/OTA catalogs raw  ")
#pl.show(block=True)
#pl.close()

x, y, z = podi_fixwcs.shift_align_wcs(starcat[:,0], starcat[:,1], reference[:,0], reference[:,1], 
                                      max_offset=0.2,
                                      verbose=False)
print x, y, z
print "MISALIGNMENT (no scatter): ",(random_shift-[x,y])*3600

x, y, z = podi_fixwcs.shift_align_wcs(scattered_starcat[:,0], scattered_starcat[:,1], reference[:,0], reference[:,1], 
                                      max_offset=0.2,
                                      verbose=False)
print x, y, z
print "MISALIGNMENT (w/ scatter): ",(random_shift-[x,y])*3600


# Now correct the catalog with the initial shift

scattered_corrected = scattered_starcat + numpy.array([x, y])





#pl.plot(reference[:,0], reference[:,1], "ro",
#        scattered_corrected[:,0], scattered_corrected[:,1], 'b.')
#pl.title("Ref/OTA catalogs with best-guess applied ")
#pl.show(block=True)


median, matched, matches = podi_fixwcs.refine_wcs_shift(ref_x=reference[:,0], ref_y=reference[:,1],
                                                        ota_x=scattered_starcat[:,0], ota_y=scattered_starcat[:,1], 
                                                        best_guess=numpy.array([x,y]),
                                                        verbose=True
                                                        )
                  
print median, median*3600

#print matched
pl.close()
pl.plot(reference[:,0], reference[:,1], "ro",
        matched[:,0], matched[:,1], 'b.')
pl.title("Ref/OTA after best-guess shift applied")
pl.show(block=True)

pl.close()
pl.plot(matches[:,0]*3600., matches[:,1]*3600., "ro")
pl.grid(True)
pl.show(block=True)


total_shift = [x,y] + median
print "MISALIGNMENT (after corr): ",(random_shift-total_shift)*3600.


           
