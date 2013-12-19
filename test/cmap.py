#!/usr/bin/env python

import pylab as m
import time

cdict = {
'red'  :  ((0., 0., 0.), (0.5, 0.25, 0.25), (1., 1., 1.)),
'green':  ((0., 1., 1.), (0.7, 0.0, 0.5), (1., 1., 1.)),
'blue' :  ((0., 1., 1.), (0.5, 0.0, 0.0), (1., 1., 1.))
}
#generate the colormap with 1024 interpolated values
my_cmap = m.matplotlib.colors.LinearSegmentedColormap('my_colormap', cdict, 1024)

#create a gaussian
x = m.arange(0, 50 , 1 , m.Float)
y = x[:,m.NewAxis]
x0 = y0 = x.size // 2
fwhm= x0/1.2
z = m.exp(-4*m.log(2)*((x-x0)**2+(y-y0)**2)/fwhm**2)

pcolormesh(z, cmap = my_cmap)
colorbar()

time.sleep(5)
