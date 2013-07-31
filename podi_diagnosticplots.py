#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy


import matplotlib
import matplotlib.pyplot
import matplotlib.colors
from podi_plotting import *

gain_correct_frames = False

import matplotlib.lines


def wcsdiag_scatter(matched_cat, filename):

    print "Creating some diagnostic plots"
    # Create some plots for WCS diagnosis
    fig = matplotlib.pyplot.figure()
    # matches_zeroed = matches - wcs_shift_refinement
    matches_zeroed = matched_cat

    count, xedges, yedges = numpy.histogram2d(matches_zeroed[:,0]*3600., matches_zeroed[:,1]*3600.,
                                              bins=[60,60], range=[[-3,3], [-3,3]])
    img = matplotlib.pyplot.imshow(count.T, 
                                   extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), 
                                   origin='lower', 
                                   cmap=cmap_bluewhite)
    # interpolation='nearest', 
    fig.colorbar(img)

    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    d_ra  = good_matches[:,0] - good_matches[:,2]
    d_dec = good_matches[:,1] - good_matches[:,3]
    matplotlib.pyplot.plot(d_ra*3600., d_dec*3600., "b,", linewidth=0)
    matplotlib.pyplot.title("WCS Scatter")
    matplotlib.pyplot.xlabel("error RA ['']")
    matplotlib.pyplot.ylabel("error DEC ['']")
    matplotlib.pyplot.xlim((-3,3))
    matplotlib.pyplot.ylim((-3,3))
    matplotlib.pyplot.grid(True)
    matplotlib.pyplot.axes().set_aspect('equal')

    if (filename == None):
        fig.show()
    else:
        png_wcsscatter = filename
        fig.savefig(png_wcsscatter)

    matplotlib.pyplot.close()


def wcsdiag_shift(matched_cat, filename):
    matches_zeroed = matched_cat
    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    d_ra  = good_matches[:,0] - good_matches[:,2]
    d_dec = good_matches[:,1] - good_matches[:,3]

    fig = matplotlib.pyplot.figure()
    #matched_zeroed = matched_cat
    #matched_zeroed[:,0:2] -= wcs_shift_refinement
    #matplotlib.pyplot.plot(matched_zeroed[:,0], matched_zeroed[:,1], ",", color=(1,1,1))
    #matching_radius_arcsec = 3. / 3600.
    #valid = (numpy.fabs(matched_zeroed[:,0]-matched_zeroed[:,2]) < matching_radius_arcsec) & \
    #    (numpy.fabs(matched_zeroed[:,1]-matched_zeroed[:,3]) < matching_radius_arcsec)
    #matched_zeroed = matched_zeroed[valid]
    
    ramin, ramax = numpy.min(matches_zeroed[:,0]), numpy.max(matches_zeroed[:,0])
    decmin, decmax = numpy.min(matches_zeroed[:,1]), numpy.max(matches_zeroed[:,1])

    dimension = numpy.min([ramax-ramin, decmax-decmin])
    vector_scaling = 2 * dimension/100 * 3600. # size of 1 arcsec in percent of screen size
    #vector_scaling = 0.02 * 3600.
    matplotlib.pyplot.quiver(good_matches[:,0], good_matches[:,1], 
                             d_ra*vector_scaling, d_dec*vector_scaling,
                             linewidth=0, angles='xy', scale_units='xy', scale=1, pivot='middle')
    # Determine min and max values
    #ramin, ramax = numpy.min(good_matches[:,0]), numpy.max(good_matches[:,0])
    #decmin, decmax = numpy.min(good_matches[:,1]), numpy.max(good_matches[:,1])
    matplotlib.pyplot.title("WCS misalignment")
    matplotlib.pyplot.xlim((ramin-0.02, ramax+0.02))
    matplotlib.pyplot.ylim((decmin-0.02, decmax+0.02))
    matplotlib.pyplot.xlabel("RA [degrees]")
    matplotlib.pyplot.xlabel("DEC [degrees]")
    # draw some arrow to mark how long the other arrows are
    arrow_x = ramin  + 0.05 * (ramax - ramin)
    arrow_y = decmax - 0.05 * (decmax - decmin)
    arrow_length = 1 / 3600. # arcsec
    matplotlib.pyplot.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
                             angles='xy', scale_units='xy', scale=1, pivot='middle',
                             headwidth=0)


    # add a line
    x,y = numpy.array([[arrow_x-arrow_length*vector_scaling, arrow_x+arrow_length*vector_scaling], [arrow_y, arrow_y]])
    matplotlib.pyplot.plot(x,y, linewidth=3, color='black')

    # add label saying "2''"
    matplotlib.pyplot.text(arrow_x, arrow_y-2*vector_scaling/3600., "%d''" % (2*arrow_length*3600), 
                           horizontalalignment='center')
#    matplotlib.pyplot.text(arrow_x, arrow_y, "%d''" % (2*arrow_length*3600), 
#                           horizontalalignment='center', verticalalignment='top')

    if (filename == None):
        fig.show()
    else:
        png_wcsdirection = filename
        fig.savefig(png_wcsdirection)

    matplotlib.pyplot.close()



    return 0



if __name__ == "__main__":

    matched_cat = numpy.loadtxt("odi+2mass.matched")

    filename = "output"

#    make_plots(matched_cat, filename)
    wcsdiag_scatter(matched_cat, filename+".wcs1.png")
    wcsdiag_shift(matched_cat, filename+".wcs2.png")
#    wcsdiag_shift(matched_cat, None)
