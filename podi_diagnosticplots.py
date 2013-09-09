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

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
import math

import matplotlib
import matplotlib.pyplot
import matplotlib.colors
from podi_plotting import *

gain_correct_frames = False

import matplotlib.lines

from podi_definitions import *

def wcsdiag_scatter(matched_cat, filename):

    stdout_write("Creating the WCS scatter plot ...")
    # Create some plots for WCS diagnosis
    fig = matplotlib.pyplot.figure()
    # matches_zeroed = matches - wcs_shift_refinement
    matches_zeroed = matched_cat

    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    d_ra  = good_matches[:,0] - good_matches[:,2]
    d_dec = good_matches[:,1] - good_matches[:,3]

    count, xedges, yedges = numpy.histogram2d(d_ra*3600., d_dec*3600.,
                                              bins=[60,60], range=[[-3,3], [-3,3]])
    img = matplotlib.pyplot.imshow(count.T, 
                                   extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), 
                                   origin='lower', 
                                   cmap=cmap_bluewhite)
    # interpolation='nearest', 
    fig.colorbar(img)

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
    stdout_write(" done!\n")


def wcsdiag_shift(matched_cat, filename):

    stdout_write("Creating the WCS offset/shift plot ...")

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
    stdout_write(" done!\n")

    return 0


def photocalib_zeropoint(odi_mag, odi_magerr, sdss_mag, sdss_magerr, output_filename,
                         zp_median, zp_std, 
                         sdss_filtername, odi_filtername,
                         zp_upper1sigma=None, zp_lower1sigma=None,
                         zp_distribfull=None, zp_distribclipped=None,
                         title=None
                         ):


    stdout_write("Creating the photometric calibration plot ...")

    fig, ax = matplotlib.pyplot.subplots()

    zp_raw = sdss_mag - odi_mag

    zp_clipped, clipped = three_sigma_clip(zp_raw, return_mask=True)

    delta = 0.3 if zp_std < 0.3 else zp_std
    close_to_median = (zp_raw > zp_median - 3 * delta) & (zp_raw < zp_median + 3 * delta)

    zp_max, zp_min = zp_median+3*delta, zp_median-3*delta

    # Determine the min and max sdss magnitudes
    sdss_min = numpy.min(sdss_mag - sdss_magerr) - 1
    sdss_max = numpy.max(sdss_mag + sdss_magerr)
    # now round them to the nearest integer
    sdss_minint = int(math.floor(sdss_min))
    sdss_maxint = int(math.ceil(sdss_max))

    # Prepare some helpers for horizontal lines
    # There has to be a more elegant way to do this
    x_values = numpy.linspace(sdss_minint, sdss_maxint, 10)
    y_values = numpy.ones(shape=x_values.shape)

    #
    # Draw a grey-shaded region outlining the 1-sigma range
    #
    ax.barh(zp_median-zp_std, (sdss_maxint-sdss_minint), height=2*zp_std, 
            left=sdss_minint, label=u"1 sigma range", 
            color="#a0a0a0", edgecolor='#a0a0a0')


    #
    # Compute and plot a histogram showing the distribution of ZPs
    #

    # Prepare a histogram to illustrate the distribution of ZP values
    binwidth = (0.2 * zp_std)
    nbins = (zp_max - zp_min) / binwidth
    count, edges = numpy.histogram(zp_raw, bins=nbins, range=[zp_min, zp_max], normed=True)
    # Normalize histogram
    count = count * (zp_std*math.sqrt(2*math.pi)) / numpy.sum(count) / binwidth

    # also add a gaussian with the determined distribution overlayed on the histogram
    gauss_y = numpy.linspace(zp_min, zp_max, 1000)
    gauss_x = float(sdss_minint) + numpy.exp(-(gauss_y-zp_median)**2/(2*zp_std**2)) #/(zp_std*math.sqrt(2*math.pi))
    #/(zp_std*math.sqrt(2*math.pi))*binwidth*numpy.sum(clipped)

    ax.barh(edges[:-1], count*1.5, height=binwidth, left=sdss_minint, 
                           edgecolor='green', color='green',
                           label='distribution')

    matplotlib.pyplot.plot(gauss_x, gauss_y, linewidth=1, ls='-', color='black')
    #print gauss_y, gauss_x


    #
    # Plot the actual measurments and values deemed to be outliers
    #
    ax.plot(sdss_mag[clipped==False], zp_raw[clipped==False], "ro", fillstyle='none', label='outliers')
    ax.plot(sdss_mag[clipped], zp_raw[clipped], "bo", linewidth=0, label='valid')

    #
    # Overplot a white line to illustrate the median value
    #
    matplotlib.pyplot.plot(x_values+1, y_values*zp_median, linewidth=1, ls='-', color='white')

    ax.grid(True)
    ax.legend(loc='upper left', borderaxespad=1)

    #matplotlib.pyplot.plot(x_values+1, y_values*zp_median, linewidth=2, ls='-', color='blue')

    if (title == None):
        title_string = u"Photometric zeropoint: %.3f \u00b1 %.3f mag" % (zp_median, zp_std)
    else:
        title_string = u"%s\nPhotometric zeropoint: %.3f \u00b1 %.3f mag" % (title, zp_median, zp_std)

    matplotlib.pyplot.title(title_string)
    matplotlib.pyplot.xlabel("SDSS magnitude in %s" % (sdss_filtername), labelpad=7)
    matplotlib.pyplot.ylabel("zeropoint (sdss [%s] - odi [%s])" % (sdss_filtername, odi_filtername), labelpad=20)
    matplotlib.pyplot.xlim((sdss_minint, sdss_maxint))
    matplotlib.pyplot.ylim((zp_min, zp_max))
#    matplotlib.pyplot.axes().set_aspect('equal')

    if (output_filename == None):
        fig.show()
    else:
        fig.savefig(output_filename)

    matplotlib.pyplot.close()
    stdout_write(" done!\n")



if __name__ == "__main__":

    matched_cat = numpy.loadtxt("odi+2mass.matched")

    filename = "output"

#    make_plots(matched_cat, filename)
    wcsdiag_scatter(matched_cat, filename+".wcs1.png")
    wcsdiag_shift(matched_cat, filename+".wcs2.png")
#    wcsdiag_shift(matched_cat, None)
