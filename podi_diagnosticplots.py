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



#####################################################################################
#
#
# This creates the WCS scatter plots
#
#
#####################################################################################



def plot_wcsdiag_scatter(d_ra, d_dec, filename, extension_list, 
                         title="WCS Scatter",
                         ota_stats = None, ota_global_stats = None):

    #fig = matplotlib.pyplot.figure()
    fig, ax = matplotlib.pyplot.subplots()

    count, xedges, yedges = numpy.histogram2d(d_ra*3600., d_dec*3600.,
                                              bins=[60,60], range=[[-3,3], [-3,3]])
    img = ax.imshow(count.T, 
                                   extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]), 
                                   origin='lower', 
                                   cmap=cmap_bluewhite)
    # interpolation='nearest', 
    fig.colorbar(img)

    if (not ota_global_stats == None):
        x = ota_global_stats['MEDIAN-RA']
        y = ota_global_stats['MEDIAN-DEC']
        width = ota_global_stats['RMS-RA']
        height = ota_global_stats['RMS-DEC']
        ellipse = matplotlib.patches.Ellipse(xy=(x,y), width=width, height=height, 
                                            edgecolor='r', fc='None', lw=2)
        ax.add_patch(ellipse)

        global_text = """\
Overall WCS:
offset: %+0.3f'' / %+0.3f''
R.M.S. %0.3f'' / %.3f''
%.3f'' (combined)""" % (ota_global_stats['MEDIAN-RA'], ota_global_stats['MEDIAN-DEC'],
            ota_global_stats['RMS-RA'], ota_global_stats['RMS-DEC'], ota_global_stats['RMS'])
        ax.text(2.9, -2.9, global_text,
        horizontalalignment='right',
        verticalalignment='bottom',
                 fontsize=10, backgroundcolor='white')

    if (not ota_stats == None):
        x = ota_stats['MEDIAN-RA']
        y = ota_stats['MEDIAN-DEC']
        width = ota_stats['RMS-RA']
        height = ota_stats['RMS-DEC']
        ellipse = matplotlib.patches.Ellipse(xy=(x,y), width=width, height=height, 
                                            edgecolor='b', fc='None', lw=2)
        ax.add_patch(ellipse)

        local_text = """\
This OTA:
offset: %+0.3f'' / %+0.3f (DEC)
R.M.S. %0.3f / %.3f (DEC) 
%.3f'' (combined)""" % (ota_stats['MEDIAN-RA'], ota_stats['MEDIAN-DEC'],
            ota_stats['RMS-RA'], ota_stats['RMS-DEC'], ota_stats['RMS'])
        ax.text(-2.9, -2.9, local_text,
                 horizontalalignment='left',
                 verticalalignment='bottom',
                 fontsize=10, backgroundcolor='white')


    ax.plot(d_ra*3600., d_dec*3600., "b,", linewidth=0)
    ax.set_title(title)
    ax.set_xlabel("error RA ['']")
    ax.set_ylabel("error DEC ['']")
    ax.set_xlim((-3,3))
    ax.set_ylim((-3,3))
    ax.grid(True)
    ax.set_aspect('equal')

    if (filename == None):
        fig.show()
    else:
        for ext in extension_list:
            fig.savefig(filename+"."+ext)

    matplotlib.pyplot.close()

    return


def wcsdiag_scatter(matched_cat, filename, options=None, ota_wcs_stats=None,
                    also_plot_singleOTAs=True
                    ):

    stdout_write("Creating the WCS scatter plot ...")

    matches_zeroed = matched_cat
    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    d_ra  = good_matches[:,0] - good_matches[:,2]
    d_dec = good_matches[:,1] - good_matches[:,3]
    ota = good_matches[:,10]
    print ota

    if (options == None):
        extension_list = ('png')
    else:
        extension_list = options['plotformat']

    # Create one plot for the full focalplane
    ota_global_stats = None if ota_wcs_stats == None else ota_wcs_stats['full']
    title = "WCS Scatter - full focal plane"
    plot_wcsdiag_scatter(d_ra, d_dec, filename, extension_list, 
                         title=title,
                         ota_stats=None, ota_global_stats=ota_global_stats)

    # Now break down the plots by OTA
    if (also_plot_singleOTAs):
        list_of_otas = available_ota_coords
        if (options['central_only']):
            list_of_otas = central_array_ota_coords

        for (otax, otay) in list_of_otas:
            this_ota = otax * 10 + otay
            in_this_ota = (ota == this_ota)

            extname = "OTA%02d.SCI" % (this_ota)
            ota_stats = None if ota_wcs_stats == None else ota_wcs_stats[extname]
            title = "WSC Scatter - OTA %02d" % (this_ota)
                
            ota_stats = None if ota_wcs_stats == None else ota_wcs_stats[extname]
            print extname, ota_stats

            ota_plotfile = "%s_OTA%02d" % (filename, this_ota)
            plot_wcsdiag_scatter(d_ra[in_this_ota], d_dec[in_this_ota], ota_plotfile, extension_list,
                                 title=title,
                                 ota_stats=ota_stats, ota_global_stats=ota_global_stats)

    # Create some plots for WCS diagnosis
    stdout_write(" done!\n")




















#####################################################################################
#
#
# WCS shift vector plots
#
#
#####################################################################################


def plot_wcsdiag_shift(matched_cat, filename, extension_list=('png'), 
                       ota_outlines=None, 
                       title=None,
                       ):

    fig, ax = matplotlib.pyplot.subplots()
    
    d_ra  = matched_cat[:,0] - matched_cat[:,2]
    d_dec = matched_cat[:,1] - matched_cat[:,3]
    ramin, ramax = numpy.min(matched_cat[:,0]), numpy.max(matched_cat[:,0])
    decmin, decmax = numpy.min(matched_cat[:,1]), numpy.max(matched_cat[:,1])

    dimension = numpy.min([ramax-ramin, decmax-decmin])
    vector_scaling = 10 * dimension/100 * 3600. # size of 1 arcsec in percent of screen size

    ax.plot(matched_cat[:,0], matched_cat[:,1],
            color='red', marker='o', linestyle='None',
            markeredgecolor='none', markersize=4, 
            )
    Q = ax.quiver(matched_cat[:,0], matched_cat[:,1], 
                  d_ra*vector_scaling, d_dec*vector_scaling,
                  linewidth=0, edgecolor='red',
                  angles='xy', scale_units='xy', pivot='tail', zorder=99, 
                  scale=1,
                  headwidth=2, headlength=2, headaxislength=1.8,
                  )
    # Determine min and max values

    ax.set_title("WCS misalignment")
    ax.set_xlim((ramin-0.02, ramax+0.02))
    ax.set_ylim((decmin-0.02, decmax+0.02))
    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")

    # draw some arrow to mark how long the other arrows are
    #arrow_x = ramin  + 0.05 * (ramax - ramin)
    #arrow_y = decmax - 0.05 * (decmax - decmin)
    arrow_length = 1 / 3600. # arcsec
    # ax.plot(arrow_x, arrow_y, "ro")
    # ax.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
    #           angles='xy', scale_units='xy', scale=1, pivot='middle',
    #           headwidth=0)

    print "Adding quiverkey at pos"#, arrow_x, arrow_y
    # qk = ax.quiverkey(Q, arrow_x, arrow_y, arrow_length, 
    #                   label="1''", labelpos='S',
    #                   coordinates='data'
    #                   )
    qk = ax.quiverkey(Q, 0.05, 0.95, arrow_length*vector_scaling, 
                      label="1''", labelpos='S',
                      coordinates='axes'
                      )

    # ax.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
    #           angles='xy', scale_units='xy', scale=1, pivot='middle',
    #           headwidth=0)

    if (not ota_outlines == None):
        corners = numpy.array(ota_outlines)
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
        ax.add_collection(coll)

    # add a line
    # x,y = numpy.array([[arrow_x-arrow_length*vector_scaling, arrow_x+arrow_length*vector_scaling], [arrow_y, arrow_y]])
    # ax.plot(x,y, linewidth=3, color='black')

    # add label saying "2''"
    # ax.text(arrow_x, arrow_y-2*vector_scaling/3600., "%d''" % (2*arrow_length*3600), 
    #                        horizontalalignment='center')

    if (filename == None):
        fig.show()
    else:
        for ext in extension_list:
            fig.savefig(filename+"."+ext)

    matplotlib.pyplot.close()
    stdout_write(" done!\n")

    return



def wcsdiag_shift(matched_cat, filename, options=None, ota_outlines=None, 
                  ota_wcs_stats=None,
                  also_plot_singleOTAs=True
                  ):

    stdout_write("Creating the WCS offset/shift plot ...")

    matches_zeroed = matched_cat
    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    d_ra  = good_matches[:,0] - good_matches[:,2]
    d_dec = good_matches[:,1] - good_matches[:,3]
    ota = good_matches[:,10]
    
    if (options == None):
        extension_list = ('png')
    else:
        extension_list = options['plotformat']

    # Create one plot for the full focalplane
    title = "WCS Scatter - full focal plane"
    plot_wcsdiag_shift(good_matches, filename, extension_list, 
                       ota_outlines=None,
                       title=title,
                       )

    # Now break down the plots by OTA
    if (also_plot_singleOTAs):
        list_of_otas = available_ota_coords
        if (options['central_only']):
            list_of_otas = central_array_ota_coords

        for (otax, otay) in list_of_otas:
            this_ota = otax * 10 + otay
            in_this_ota = (ota == this_ota)

            extname = "OTA%02d.SCI" % (this_ota)
            title = "WSC Scatter - OTA %02d" % (this_ota)
                
            ota_plotfile = "%s_OTA%02d" % (filename, this_ota)
            plot_wcsdiag_shift(good_matches[in_this_ota], 
                               filename=ota_plotfile, 
                               extension_list=extension_list, 
                               ota_outlines=None,
                               title=title,
                               )

#     fig, ax = matplotlib.pyplot.subplots()
    
#     #matched_zeroed = matched_cat
#     #matched_zeroed[:,0:2] -= wcs_shift_refinement
#     #matplotlib.pyplot.plot(matched_zeroed[:,0], matched_zeroed[:,1], ",", color=(1,1,1))
#     #matching_radius_arcsec = 3. / 3600.
#     #valid = (numpy.fabs(matched_zeroed[:,0]-matched_zeroed[:,2]) < matching_radius_arcsec) & \
#     #    (numpy.fabs(matched_zeroed[:,1]-matched_zeroed[:,3]) < matching_radius_arcsec)
#     #matched_zeroed = matched_zeroed[valid]
    
#     ramin, ramax = numpy.min(matches_zeroed[:,0]), numpy.max(matches_zeroed[:,0])
#     decmin, decmax = numpy.min(matches_zeroed[:,1]), numpy.max(matches_zeroed[:,1])

#     dimension = numpy.min([ramax-ramin, decmax-decmin])
#     vector_scaling = 2 * dimension/100 * 3600. # size of 1 arcsec in percent of screen size
#     #vector_scaling = 0.02 * 3600.
#     matplotlib.pyplot.quiver(good_matches[:,0], good_matches[:,1], 
#                              d_ra*vector_scaling, d_dec*vector_scaling,
#                              linewidth=0, angles='xy', scale_units='xy', scale=1, pivot='middle')
#     # Determine min and max values
#     #ramin, ramax = numpy.min(good_matches[:,0]), numpy.max(good_matches[:,0])
#     #decmin, decmax = numpy.min(good_matches[:,1]), numpy.max(good_matches[:,1])
#     matplotlib.pyplot.title("WCS misalignment")
#     matplotlib.pyplot.xlim((ramin-0.02, ramax+0.02))
#     matplotlib.pyplot.ylim((decmin-0.02, decmax+0.02))
#     matplotlib.pyplot.xlabel("RA [degrees]")
#     matplotlib.pyplot.ylabel("DEC [degrees]")
#     # draw some arrow to mark how long the other arrows are
#     arrow_x = ramin  + 0.05 * (ramax - ramin)
#     arrow_y = decmax - 0.05 * (decmax - decmin)
#     arrow_length = 1 / 3600. # arcsec
#     matplotlib.pyplot.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
#                              angles='xy', scale_units='xy', scale=1, pivot='middle',
#                              headwidth=0)

#     if (not ota_outlines == None):
#         corners = numpy.array(ota_outlines)
#         coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
#         ax.add_collection(coll)

#     # add a line
#     x,y = numpy.array([[arrow_x-arrow_length*vector_scaling, arrow_x+arrow_length*vector_scaling], [arrow_y, arrow_y]])
#     matplotlib.pyplot.plot(x,y, linewidth=3, color='black')

#     # add label saying "2''"
#     matplotlib.pyplot.text(arrow_x, arrow_y-2*vector_scaling/3600., "%d''" % (2*arrow_length*3600), 
#                            horizontalalignment='center')
# #    matplotlib.pyplot.text(arrow_x, arrow_y, "%d''" % (2*arrow_length*3600), 
# #                           horizontalalignment='center', verticalalignment='top')

#     if (filename == None):
#         fig.show()
#     else:
#         if (not options == None):
#             for ext in options['plotformat']:
#                 if (ext != ''):
#                     fig.savefig(filename+"."+ext)
#         else:
#             fig.savefig(filename+".png")

#     matplotlib.pyplot.close()
#     stdout_write(" done!\n")

    return 0























def photocalib_zeropoint(odi_mag, odi_magerr, sdss_mag, sdss_magerr, output_filename,
                         zp_median, zp_std, 
                         sdss_filtername, odi_filtername,
                         zp_upper1sigma=None, zp_lower1sigma=None,
                         zp_distribfull=None, zp_distribclipped=None,
                         title=None,
                         options=None,
                         also_plot_singleOTAs=False,
                         ):


    stdout_write("Creating the photometric calibration plot ...")

    fig, ax = matplotlib.pyplot.subplots()

    zp_raw = sdss_mag - odi_mag
    zp_err = numpy.hypot(sdss_magerr, odi_magerr)

    if (zp_raw.shape[0] < 5):
        zp_clipped = zp_raw
        clipped = numpy.isfinite(zp_clipped)
    else:
        zp_clipped, clipped = three_sigma_clip(zp_raw, return_mask=True)

    delta = 0.3 if zp_std < 0.3 else zp_std
    close_to_median = (zp_raw > zp_median - 3 * delta) & (zp_raw < zp_median + 3 * delta)
    
    
    #zp_max, zp_min = zp_median+3*delta, zp_median-3*delta
    zp_max = zp_median+5*zp_std+0.3
    zp_min = zp_median-5*zp_std-0.3

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
    ax.errorbar(sdss_mag[clipped==False], zp_raw[clipped==False], 
                xerr=sdss_magerr[clipped==False], 
                yerr=zp_err[clipped==False], 
                capsize=0,
                fmt="r+", fillstyle='none', label='outliers')

    ax.errorbar(sdss_mag[clipped], zp_raw[clipped], 
                xerr=sdss_magerr[clipped], 
                yerr=zp_err[clipped], 
                capsize=0,
                fmt="b+", linewidth=1, label='valid')

    #
    # Overplot a white line to illustrate the median value
    #
    matplotlib.pyplot.plot(x_values+1, y_values*zp_median, linewidth=1, ls='-', color='white')

    ax.grid(True)
    ax.legend(loc='upper left', borderaxespad=1)

    photzp_text = u"ZP = %.3f \u00b1 %.3f mag" % (zp_median, zp_std)
    ax.text(0.96, 0.05, photzp_text, fontsize=15,
            horizontalalignment='right',
            verticalalignment='bottom',
            transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', linewidth=0))

    title_string = title if (title != None) else ""
    matplotlib.pyplot.title(title_string)
    matplotlib.pyplot.xlabel("SDSS magnitude in %s" % (sdss_filtername), labelpad=7)
    matplotlib.pyplot.ylabel("zeropoint (sdss [%s] - odi [%s])" % (sdss_filtername, odi_filtername), labelpad=20)
    matplotlib.pyplot.xlim((sdss_minint, sdss_maxint))
    matplotlib.pyplot.ylim((zp_min, zp_max))
#    matplotlib.pyplot.axes().set_aspect('equal')

    if (output_filename == None):
        fig.show()
    else:
        if (not options == None):
            for ext in options['plotformat']:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext)
                    print "plot saved as ",output_filename+"."+ext
        else:
            fig.savefig(output_filename+".png")

    matplotlib.pyplot.close()
    stdout_write(" done!\n")
























#####################################################################################
#
#
# Following are the routines to create the zeropoint maps
#
#
#####################################################################################

def plot_zeropoint_map(ra, dec, zp, ota_outlines, output_filename, options, zp_range):

    fig, ax = matplotlib.pyplot.subplots()
    zp_min, zp_max = zp_range

    if (ra.shape[0] < 5):
        zp_clipped = zp
        clipped = numpy.isfinite(zp)
    else:
        zp_clipped, clipped = three_sigma_clip(zp, return_mask=True)

    zp_median_ota = numpy.median(zp_clipped)
    zp_std_ota = numpy.std(zp_clipped)

    #print ra_ota
    #print zp_ota

    sc = matplotlib.pyplot.scatter(ra, dec, c=zp, alpha=0.75, 
                                   vmin=zp_min, vmax=zp_max, edgecolor='none',
                                   s=9, cmap=matplotlib.pyplot.cm.get_cmap('spectral'))
        #matplotlib.pyplot.colorbar(sc, cm)

    if (not ota_outlines == None):
        corners = numpy.array(ota_outlines)
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
        ax.add_collection(coll)

    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    ax.autoscale_view()

    colorbar = matplotlib.pyplot.colorbar(cmap=matplotlib.pyplot.cm.get_cmap('spectral'))
    colorbar.set_label("phot. zeropoint")

    if (output_filename == None):
        fig.show()
    else:
        if (not options == None):
            for ext in options['plotformat']:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext)
                    print "saving plot to",output_filename+"."+ext
        else:
            fig.savefig(output_filename+".png")

    matplotlib.pyplot.close()
    return


def photocalib_zeropoint_map(odi_mag, sdss_mag, ota, ra, dec, output_filename,
                             sdss_filtername, odi_filtername,
                             title=None,
                             ota_outlines=None,
                             options=None,
                             also_plot_singleOTAs=True
                             ):


    stdout_write("Creating the photometric calibration per OTA plot ...")

    zp_raw = sdss_mag - odi_mag
    zp_min = scipy.stats.scoreatpercentile(zp_raw,  5)
    zp_max = scipy.stats.scoreatpercentile(zp_raw, 95)

    # Create one plot for the full focal plane, using boxes to outlines OTAs
    zp_range = (zp_min, zp_max)
    plot_zeropoint_map(ra, dec, zp_raw, ota_outlines, output_filename, options, zp_range)
    print ota_outlines

    # If requested, do the same for the individual OTAs
    if (also_plot_singleOTAs):
        list_of_otas = available_ota_coords
        if (options['central_only']):
            list_of_otas = central_array_ota_coords

        for (otax, otay) in list_of_otas:
            this_ota = otax * 10 + otay
            in_this_ota = ota == this_ota
            zp_ota = zp_raw[in_this_ota]
            ra_ota = ra[in_this_ota]
            dec_ota = dec[in_this_ota]

            ota_plotfile = "%s_OTA%02d" % (output_filename, this_ota)
            plot_zeropoint_map(ra_ota, dec_ota, zp_ota, None, ota_plotfile, options, zp_range)
    
    stdout_write(" done!\n")

    return










#####################################################################################
#
#
# Following are the routines to create the PSF size maps
# These plots show the psf size for each source color-coded at its 
# position in the field of view
#
#
#####################################################################################

def plot_psfsize_map(ra, dec, fwhm, output_filename, 
                     fwhm_range,
                     title=None,
                     ota_outlines=None,
                     options=None,
                     ):
                       
    fig, ax = matplotlib.pyplot.subplots()

    fwhm_min, fwhm_max = fwhm_range

    sc = matplotlib.pyplot.scatter(ra, dec, c=fwhm, alpha=0.75, 
                                   vmin=fwhm_min, vmax=fwhm_max, edgecolor='none',
                                   s=9, cmap=matplotlib.pyplot.cm.get_cmap('spectral'))

    if (not ota_outlines == None):
        corners = numpy.array(ota_outlines)
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', 
                                                     linestyle='-', linewidth=0.5)
        ax.add_collection(coll)

    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    ax.autoscale_view()

    colorbar = matplotlib.pyplot.colorbar(cmap=matplotlib.pyplot.cm.get_cmap('spectral'))
    colorbar.set_label("FWHM [arcsec]")

    if (output_filename == None):
        fig.show()
    else:
        if (not options == None):
            for ext in options['plotformat']:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext)
                    print "saving plot to",output_filename+"."+ext
        else:
            fig.savefig(output_filename+".png")

    matplotlib.pyplot.close()
  

  
def diagplot_psfsize_map(ra, dec, fwhm, ota, output_filename,
                         title=None,
                         ota_outlines=None,
                         options=None,
                         also_plot_singleOTAs=True,
                         ):


    stdout_write("Creating the PSF size as heatmap plot ...")
    fwhm_scale_arcsec = 3600.0 * 0.11

    fwhm *= fwhm_scale_arcsec
    fwhm_min = scipy.stats.scoreatpercentile(fwhm,  5)
    fwhm_max = scipy.stats.scoreatpercentile(fwhm, 95)
    fwhm_clipped, clipmask = three_sigma_clip(fwhm, return_mask=True)

    fwhm_range = (fwhm_min, fwhm_max)

    # Create one plot for the full focal plane, using boxes to outlines OTAs
    plot_psfsize_map(ra, dec, fwhm, output_filename, 
                     fwhm_range=fwhm_range,
                     title=title,
                     ota_outlines=ota_outlines,
                     options=options)

    print ota_outlines

    # If requested, do the same for the individual OTAs
    if (also_plot_singleOTAs):
        list_of_otas = available_ota_coords
        if (options['central_only']):
            list_of_otas = central_array_ota_coords

        for (otax, otay) in list_of_otas:
            this_ota = otax * 10 + otay

            in_this_ota = (ota == this_ota)
            fwhm_ota = fwhm[in_this_ota]
            ra_ota = ra[in_this_ota]
            dec_ota = dec[in_this_ota]

            ota_plotfile = "%s_OTA%02d" % (output_filename, this_ota)
            plot_psfsize_map(ra_ota, dec_ota, fwhm_ota, 
                             output_filename=ota_plotfile, 
                             fwhm_range=fwhm_range,
                             title=title,
                             ota_outlines=None,
                             options=options)

    stdout_write(" done!\n")

    return





#####################################################################################
#
#
# Stand-alone routines are below.
#
#
#####################################################################################





if __name__ == "__main__":

    if (cmdline_arg_isset("-psfplot")):
        fitsfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(fitsfile)
        catalogfile = fitsfile+".src.cat"
        source_cat = numpy.loadtxt(catalogfile)

        from podi_collectcells import read_options_from_commandline
        options = read_options_from_commandline(None)

        flags = source_cat[:,7]
        valid_flags = flags == 0

        ra = source_cat[:,0][valid_flags]
        dec = source_cat[:,1][valid_flags]
        fwhm = source_cat[:,5][valid_flags]
        ota = source_cat[:,8][valid_flags]
        print fwhm

        ota_outlines = derive_ota_outlines(hdulist)

        output_filename = fitsfile[:-5]+".seeing"
        output_filename = "test.seeing"
        diagplot_psfsize_map(ra, dec, fwhm, ota, output_filename,
                             title="PSF size",
                             ota_outlines=ota_outlines,
                             options=options,
                             also_plot_singleOTAs=True
                             )
        sys.exit(0)

    if (cmdline_arg_isset("-zeropoint_map")):
        fitsfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(fitsfile)
        filtername = hdulist[0].header['FILTER']
        print filtername

        from podi_collectcells import read_options_from_commandline
        options = read_options_from_commandline(None)

        ota_outlines = derive_ota_outlines(hdulist)

        import podi_photcalib
        source_cat = numpy.loadtxt(fitsfile+".src.cat")
        podi_photcalib.photcalib(source_cat, "test.fits",
                                 filtername, exptime=1,
                                 diagplots=True,
                                 calib_directory=None,
                                 overwrite_cat=None,
                                 plottitle="standalone",
                                 otalist=ota_outlines,
                                 options=options)

        sys.exit(0)


   



    matched_cat = numpy.loadtxt("odi+2mass.matched")

    filename = "output"

#    make_plots(matched_cat, filename)
    wcsdiag_scatter(matched_cat, filename+".wcs1")
    wcsdiag_shift(matched_cat, filename+".wcs2")
#    wcsdiag_shift(matched_cat, None)


    
