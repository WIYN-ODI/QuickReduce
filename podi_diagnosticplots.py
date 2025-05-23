#!/usr/bin/env python3
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

Module containing all code pertaining to creating the diagnostic plots.

"""

import sys
import os
import astropy.io.fits as pyfits
import numpy
import math

import matplotlib
import matplotlib.colors
# from podi_plotting import *

gain_correct_frames = False

import matplotlib.lines

from podi_definitions import *
from podi_commandline import *
import podi_sitesetup as sitesetup
import podi_focalplanelayout

import queue
import multiprocessing
import podi_logging, logging

import matplotlib.patches

import pickle
import time
from multiprocessing.queues import _ForkingPickler
#####################################################################################
#
#
# This creates the WCS scatter plots
#
#
#####################################################################################



def plot_wcsdiag_scatter(d_ra, d_dec, filename, extension_list,
                         title="WCS Scatter",
                         high_s2n=None,
                         ota_stats = None, ota_global_stats = None):
    """

    Function generating the WCS scatter plot from the given data.

    """
    
    logger = logging.getLogger("DiagPlot_WCSScatter")
    import matplotlib.pyplot

    #fig = matplotlib.pyplot.figure()
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    count, xedges, yedges = numpy.histogram2d(d_ra*3600., d_dec*3600.,
                                              bins=[60,60], range=[[-3,3], [-3,3]])
    # img = ax.imshow(count.T,
    #                                extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),
    #                                origin='lower',
    #                                cmap=cmap_bluewhite)
    # # interpolation='nearest',
    # fig.colorbar(img, label='point density')

    max_dimension = 3.
    n_sigma = 7
    if (ota_global_stats is not None):
        max_dimension = numpy.max([n_sigma*ota_global_stats['RMS-RA-CLIP'],
                                   n_sigma*ota_global_stats['RMS-DEC-CLIP'],
                                   0.5,])
        max_dimension = 3. if max_dimension > 3 else max_dimension
        #print "max:",max_dimension

    #
    # Draw the rings in the center to outline the RMS of the OTA and focal-plane
    #
    if (ota_global_stats is not None):
        x = ota_global_stats['MEDIAN-RA']
        y = ota_global_stats['MEDIAN-DEC']
        width = ota_global_stats['RMS-RA']
        height = ota_global_stats['RMS-DEC']
        ellipse = matplotlib.patches.Ellipse(xy=(x,y), width=width, height=height, 
                                             edgecolor='r', fc='None', lw=1)
        ax.add_patch(ellipse)

        global_text = """\
Overall WCS (N=%(STARCOUNT)d):
offset: %(MEDIAN-RA)+0.3f'' / %(MEDIAN-DEC)+0.3f''
R.M.S. %(RMS-RA-CLIP)0.3f'' / %(RMS-DEC-CLIP)0.3f''
%(RMS-CLIP).3f'' (combined)\
""" % ota_global_stats
        ax.text(-0.97*max_dimension, -0.97*max_dimension, global_text,
        horizontalalignment='right',
        verticalalignment='bottom',
                 fontsize=10, backgroundcolor='white')

    if (ota_stats is not None):
        x = ota_stats['MEDIAN-RA']
        y = ota_stats['MEDIAN-DEC']
        width = ota_stats['RMS-RA']
        height = ota_stats['RMS-DEC']
        ellipse = matplotlib.patches.Ellipse(xy=(x,y), width=width, height=height, 
                                             edgecolor='b', fc='None', lw=1)
        ax.add_patch(ellipse)

        local_text = """\
This OTA (N=%(STARCOUNT)d):
offset: %(MEDIAN-RA)+0.3f'' / %(MEDIAN-DEC)+0.3f''
R.M.S. %(RMS-RA-CLIP)0.3f'' / %(RMS-DEC-CLIP)0.3f''
%(RMS-CLIP).3f'' (combined)\
""" % ota_stats
        ax.text(0.97*max_dimension, -0.97*max_dimension, local_text,
                 horizontalalignment='left',
                 verticalalignment='bottom',
                 fontsize=10, backgroundcolor='white')

    #
    # Add some histograms to the borders to illustrate the distribution
    # Only do so if there are at least 5 stars
    #
    histogram_scale = 0.3*max_dimension
    select = [
        (numpy.ones(d_ra.shape, dtype=numpy.bool), '#808080')
    ]
    if (high_s2n is not None):
        select.append((high_s2n, 'blue'))

    if (d_ra.shape[0] > 5):
        from scipy.stats import gaussian_kde
        x = numpy.linspace(-3,3,600)

        for (sample, color) in select:
            try:
                density_ra = gaussian_kde(d_ra[sample]*3600.)
                density_ra.covariance_factor = lambda : .1
                density_ra._compute_covariance()
                peak_ra = numpy.max(density_ra(x))
                ax.plot(x,max_dimension-density_ra(x)/peak_ra*histogram_scale, "-", color=color)#'black')

                density_dec = gaussian_kde(d_dec[sample]*3600.)
                density_dec.covariance_factor = lambda : .1
                density_dec._compute_covariance()
                peak_dec = numpy.max(density_dec(x))
                ax.plot(density_dec(x)/peak_dec*histogram_scale-max_dimension, x, "-", color=color)#'black')
            except ValueError:
                pass


    if (high_s2n is not None):
        # ax.plot(d_ra[high_s2n] * 3600., d_dec[high_s2n] * 3600., "b,", linewidth=0)
        ax.scatter(d_ra[~high_s2n] * 3600., d_dec[~high_s2n] * 3600.,
                   c='#aaaaaa', marker='.', linewidth=0, s=1)
        ax.scatter(d_ra[high_s2n] * 3600., d_dec[high_s2n] * 3600.,
                   c='blue', s=4, marker='.', linewidth=0)
    else:
        ax.plot(d_ra*3600., d_dec*3600., "b,", linewidth=0)
    ax.set_title(title)
    ax.set_xlabel("error RA * cos(DEC) [arcsec]")
    ax.set_ylabel("error DEC [arcsec]")
    ax.set_xlim((max_dimension,-1.*max_dimension))
    ax.set_ylim((-1.*max_dimension,max_dimension))
    ax.grid(True)
    ax.set_aspect('equal')

    if (filename is None):
        fig.show()
    else:
        for ext in extension_list:
            fig.set_size_inches(8,6)
            logger.debug("saving file: %s.%s" % (filename, ext))
            fig.savefig(filename+"."+ext, dpi=100,bbox_inches='tight')

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()

    return


def wcsdiag_scatter(matched_radec_odi, 
                    matched_radec_2mass, 
                    matched_ota,
                    matched_odierror,
                    filename, options=None, ota_wcs_stats=None,
                    also_plot_singleOTAs=True,
                    title_info = None,
                    ):
    """

    Function preparing the data for WCS scatter plot and forwarding all plotting
    tasks to `plot_wcsdiag_scatter`

    """

    logger = logging.getLogger("DiagPlot_WCSScatter")
    logger.info("Creating the WCS scatter plot ...")

    # Eliminate all matches with somehow illegal coordinates
    good_matches = numpy.isfinite(matched_radec_odi[:,0]) & numpy.isfinite(matched_radec_2mass[:,0])
    matched_radec_odi   = matched_radec_odi[good_matches]
    matched_radec_2mass = matched_radec_2mass[good_matches]
    matched_ota         = matched_ota[good_matches]

#    matches_zeroed = matched_cat
#    good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    cos_declination = numpy.cos(numpy.radians(0.5 * (matched_radec_odi[:,1] + matched_radec_2mass[:,1])))
    d_ra  = (matched_radec_odi[:,0] - matched_radec_2mass[:,0]) * cos_declination
    d_dec = matched_radec_odi[:,1] - matched_radec_2mass[:,1]
    ota = matched_ota

    if (options is None):
        extension_list = ['png']
    else:
        extension_list = options['plotformat']

    # Create one plot for the full focalplane
    ota_global_stats = None if ota_wcs_stats is None else ota_wcs_stats['full']

    title = "WCS Scatter - full focal plane"
    if (title_info is not None):
        logger.debug("Received information for a more descriptive plot title for WCS scatter")
        try:
            title = "WCS Scatter (reference: %(ASTRMCAT)s)\n%(OBSID)s (focal-plane) \n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec" % title_info
            logger.debug(title)
        except:
            pass

    processes = []
    plot_args= {"d_ra": d_ra, 
                "d_dec": d_dec, 
                "filename": filename, 
                "extension_list": extension_list, 
                "ota_stats": None,
                "ota_global_stats": ota_global_stats,
                "title": title,
                "high_s2n": (matched_odierror<0.02),
            }
    plot_wcsdiag_scatter(**plot_args)
    # p = multiprocessing.Process(target=plot_wcsdiag_scatter, kwargs=plot_args)
    # p.start()
    # processes.append(p)
    # plot_wcsdiag_scatter(d_ra, d_dec, filename, extension_list, 
    #                      title=title,
    #                      ota_stats=None, ota_global_stats=ota_global_stats)

    # Now break down the plots by OTA
    if (also_plot_singleOTAs):

        list_of_otas = set(ota)
        for this_ota in list_of_otas:
            in_this_ota = (ota == this_ota)
            if (numpy.sum(in_this_ota) <= 0):
                continue

            extname = "OTA%02d.SCI" % (this_ota)

            ota_stats = None
            if (ota_wcs_stats is not None):
                if (extname in ota_wcs_stats):
                    ota_stats = ota_wcs_stats[extname]

            title = "WSC Scatter - OTA %02d" % (this_ota)
            if (title_info is not None):
                title_info['OTA'] = this_ota
                try:
                    title = "WCS Scatter (reference: %(ASTRMCAT)s)\n%(OBSID)s - OTA %(OTA)02d\n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec" % title_info
                except:
                    pass
               
            logger.debug(extname+" -> "+str(ota_stats))

            ota_plotfile = create_qa_otaplot_filename(filename, this_ota, options['structure_qa_subdirs'])
            # ota_plotfile = "%s_OTA%02d" % (filename, this_ota)
            # if (options['structure_qa_subdirs']):
            #     ota_plotfile = "%s/OTA%02d" % (filename, this_ota)

            plot_args= {"d_ra": d_ra[in_this_ota], 
                        "d_dec": d_dec[in_this_ota],
                        "high_s2n": (matched_odierror[in_this_ota] < 0.02),
                        "filename": ota_plotfile, 
                        "extension_list": extension_list, 
                        "title": title,
                        "ota_stats": ota_stats,
                        "ota_global_stats": ota_global_stats,
                    }
            plot_wcsdiag_scatter(**plot_args)
            # p = multiprocessing.Process(target=plot_wcsdiag_scatter, kwargs=plot_args)
            # p.start()
            # processes.append(p)

            # plot_wcsdiag_scatter(d_ra[in_this_ota], d_dec[in_this_ota], ota_plotfile, extension_list,
            #                      title=title,
            #                      ota_stats=ota_stats, ota_global_stats=ota_global_stats)

    for p in processes:
        p.join()

    # Create some plots for WCS diagnosis
    logger.debug("done!")




















#####################################################################################
#
#
# WCS shift vector plots
#
#
#####################################################################################


def plot_wcsdiag_shift(radec, d_radec,
                       # matched_cat, 
                       filename, extension_list=['png'], 
                       ota_outlines=None, 
                       title=None,
                       ):
    """

    Function generating the WCS shift plot from the given data.

    """

    logger = logging.getLogger("DiagPlot_WCSShift")
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.ticklabel_format(useOffset=False)

#    d_ra  = matched_cat[:,0] - matched_cat[:,2]
#    d_dec = matched_cat[:,1] - matched_cat[:,3]
#    ramin, ramax = numpy.min(matched_cat[:,0]), numpy.max(matched_cat[:,0])
#    decmin, decmax = numpy.min(matched_cat[:,1]), numpy.max(matched_cat[:,1])

    around_zero = False
    if ((numpy.max(radec[:,0]) - numpy.min(radec[:,0])) > 180):
        
        # This means we most likely have to deal with coordinates around 0
        radec[:,0][radec[:,0] > 180] -= 360.
        around_zero = True

    ramin, ramax = numpy.min(radec[:,0]), numpy.max(radec[:,0])
    decmin, decmax = numpy.min(radec[:,1]), numpy.max(radec[:,1])

    dimension = numpy.min([ramax-ramin, decmax-decmin])
    vector_scaling = 10 * dimension/100 * 3600. # size of 1 arcsec in percent of screen size

    ax.plot(radec[:,0], radec[:,1],
            color='red', marker='o', linestyle='None',
            markeredgecolor='none', markersize=4, 
            )
    Q = ax.quiver(radec[:,0], radec[:,1], 
                  d_radec[:,0]*vector_scaling, d_radec[:,1]*vector_scaling,
                  angles='xy', scale_units='xy', pivot='tail', zorder=99, 
                  scale=1, 
                  headwidth=2, headlength=2, headaxislength=1.8,
                  width=1e-4, linewidth=1, edgecolor='#000000', color='#000000',
                  )
    # Determine min and max values

    ax.set_title(title)
    ax.set_xlim((ramin-0.02, ramax+0.02))
    ax.set_ylim((decmin-0.02, decmax+0.02))
    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    # Invert the x-axis to have north up and east to the left
    ax.set_xlim(ax.get_xlim()[::-1])

    # draw some arrow to mark how long the other arrows are
    #arrow_x = ramin  + 0.05 * (ramax - ramin)
    #arrow_y = decmax - 0.05 * (decmax - decmin)
    arrow_length = 1 / 3600. # arcsec
    # ax.plot(arrow_x, arrow_y, "ro")
    # ax.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
    #           angles='xy', scale_units='xy', scale=1, pivot='middle',
    #           headwidth=0)

    logger.debug("Adding quiverkey at pos") #, arrow_x, arrow_y
    # qk = ax.quiverkey(Q, arrow_x, arrow_y, arrow_length, 
    #                   label="1''", labelpos='S',
    #                   coordinates='data'
    #                   )
    # qk = ax.quiverkey(Q, 0.05, 0.95, arrow_length*vector_scaling, 
    #                   label="1''", labelpos='S',
    #                   coordinates='axes'
    #                   )

    # ax.quiver(arrow_x, arrow_y, arrow_length*vector_scaling, 0, linewidth=0,
    #           angles='xy', scale_units='xy', scale=1, pivot='middle',
    #           headwidth=0)

    if (ota_outlines is not None):
        corners = numpy.array(ota_outlines)
        if (around_zero):
            corners[:,:,0][corners[:,:,0] > 180.] -= 360.
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
        ax.add_collection(coll)

    # add a line
    # x,y = numpy.array([[arrow_x-arrow_length*vector_scaling, arrow_x+arrow_length*vector_scaling], [arrow_y, arrow_y]])
    # ax.plot(x,y, linewidth=3, color='black')

    # add label saying "2''"
    # ax.text(arrow_x, arrow_y-2*vector_scaling/3600., "%d''" % (2*arrow_length*3600), 
    #                        horizontalalignment='center')

    if (filename is None):
        fig.show()
    else:
        for ext in extension_list:
            logger.debug("saving file: %s.%s" % (filename, ext))
            fig.set_size_inches(8,6)
            fig.savefig(filename+"."+ext, dpi=100,bbox_inches='tight')

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()
    logger.debug("done!")

    return



def wcsdiag_shift(matched_radec_odi, 
                  matched_radec_2mass, 
                  matched_ota,
                  filename, 
                  options=None, 
                  ota_outlines=None, 
                  ota_wcs_stats=None,
                  also_plot_singleOTAs=True,
                  title_info = None,
                  ):
    """

    Function preparing the data for WCS shift plot and forwarding all plotting
    tasks to `plot_wcsdiag_shift`

    """

    logger = logging.getLogger("DiagPlot_WCSShift")
    logger.info("Creating the WCS offset/shift plot ...")

    # Eliminate all matches with somehow illegal coordinates
    good_matches = numpy.isfinite(matched_radec_odi[:,0]) & numpy.isfinite(matched_radec_2mass[:,0])
    matched_radec_odi   = matched_radec_odi[good_matches].copy()
    matched_radec_2mass = matched_radec_2mass[good_matches].copy()
    matched_ota         = matched_ota[good_matches].copy()

    around_zero = False
    if (((numpy.max(matched_radec_odi[:,0]) - numpy.min(matched_radec_odi[:,0])) > 180) or
        ((numpy.max(matched_radec_2mass[:,0]) - numpy.min(matched_radec_2mass[:,0])) > 180)):
        # This means we most likely have to deal with coordinates around 0
        matched_radec_odi[:,0][matched_radec_odi[:,0] > 180] -= 360.
        matched_radec_2mass[:,0][matched_radec_2mass[:,0] > 180] -= 360.
        around_zero = True

    cos_declination = numpy.cos(numpy.radians(0.5 * (matched_radec_odi[:,1] + matched_radec_2mass[:,1])))
    d_radec = matched_radec_odi - matched_radec_2mass
    d_radec[:,0] *= cos_declination
    ota = matched_ota

#    d_ra  = matched_radec_odi[:,0] - matched_radec_2mass[:,0]
#    d_dec = matched_radec_odi[:,1] - matched_radec_2mass[:,1]


    # matches_zeroed = matched_cat
    # good_matches = matches_zeroed[matches_zeroed[:,2] >= 0]
    # d_ra  = good_matches[:,0] - good_matches[:,2]
    # d_dec = good_matches[:,1] - good_matches[:,3]
    # ota = good_matches[:,10]
    
    if (options is None):
        extension_list = ['png']
    else:
        extension_list = options['plotformat']

    processes = []

    # Create one plot for the full focalplane
    title = "WCS errors - full focal plane"
    if (title_info is not None):
        logger.debug("Received information for a more descriptive plot title for WCS scatter")
        try:
            title = "WCS Errors - %(OBSID)s (focal-plane) \n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec" % title_info
            logger.debug(title)
        except:
            pass

    plot_args = {"radec": matched_radec_2mass,
                 "d_radec": d_radec,
                 "filename": filename, 
                 "extension_list": extension_list, 
                 "ota_outlines": None,
                 "title": title,
             }
    plot_wcsdiag_shift(**plot_args)
    # p = multiprocessing.Process(target=plot_wcsdiag_shift, kwargs=plot_args)
    # p.start()
    # processes.append(p)
 
    # Now break down the plots by OTA
    if (also_plot_singleOTAs):
        list_of_otas = set(ota)
        for this_ota in list_of_otas:
            in_this_ota = (ota == this_ota)
            if (numpy.sum(in_this_ota) <= 0):
                continue

            extname = "OTA%02d.SCI" % (this_ota)
            title = "WSC Error - OTA %02d" % (this_ota)
            if (title_info is not None):
                title_info['OTA'] = this_ota
                try:
                    title = "WCS Errors - %(OBSID)s OTA %(OTA)02d\n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec" % title_info
                except:
                    pass
                
            ota_radec = matched_radec_2mass[in_this_ota]
            ota_d_radec = d_radec[in_this_ota]

            ota_plotfile = create_qa_otaplot_filename(filename, this_ota, options['structure_qa_subdirs'])
            # ota_plotfile = "%s_OTA%02d" % (filename, this_ota)
            # if (options['structure_qa_subdirs']):
            #     ota_plotfile = "%s/OTA%02d" % (filename, this_ota)

            plot_args = {"radec": ota_radec,
                         "d_radec": ota_d_radec,
                         "filename": ota_plotfile, 
                         "extension_list": extension_list, 
                         "ota_outlines": None,
                         "title": title,
                     }
            # p = multiprocessing.Process(target=plot_wcsdiag_shift, kwargs=plot_args)
            # p.start()
            # processes.append(p)
            plot_wcsdiag_shift(**plot_args)

    for p in processes:
        p.join()

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























def photocalib_zeropoint(output_filename,
                         sdss_filtername, odi_filtername,
                         zp_distribfull=None, zp_distribclipped=None,
                         title=None,
                         options=None,
                         also_plot_singleOTAs=False,
                         details=None,
                         plot_formats=None,
                         ):

    """

    Generate the photometric zeropoint plot

    """

    logger = logging.getLogger("DiagPlot_PhotZP")
    logger.info("Creating the photometric calibration plot ...")

    if (details['photref_col_mag'] < 0):
        logger.debug("Error: No valid reference photometry found!")
        return

    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    #small_errors = details['odi_sdss_matched_smallerrors']

    use_for_calibration = details['use_for_calibration_mask']
    small_errors = details['odi_sdss_matched'][use_for_calibration]

    sdss_mag = small_errors[:, details['photref_col_mag']]
    sdss_magerr = small_errors[:, details['photref_col_err']]

    # zp = small_errors[:, photref_col_mag] - small_errors[:, SXcolumn[col_mag]+2]
    # zp_err = numpy.hypot(small_errors[:, photref_col_err], small_errors[:, SXcolumn[col_magerr]+2])

    # zp_raw = sdss_mag - odi_mag
    # zp_err = numpy.hypot(sdss_magerr, odi_magerr)

    # if (zp_raw.shape[0] < 5):
    #     zp_clipped = zp_raw
    #     clipped = numpy.isfinite(zp_clipped)
    # else:
    #     zp_clipped, clipped = three_sigma_clip(zp_raw, return_mask=True)

    zp_median = details['median']
    zp_std = details['std']
    delta = 0.3 if zp_std < 0.3 else zp_std
    # close_to_median = (zp_raw > zp_median - 3 * delta) & (zp_raw < zp_median + 3 * delta)
    

    #zp_max, zp_min = zp_median+3*delta, zp_median-3*delta
    zp_min = zp_median-sitesetup.diagplot__zeropoint_ZPrange[0] if \
             (sitesetup.diagplot__zeropoint_ZPrange[0] > 0) else zp_median-5*zp_std-0.3
    zp_max = zp_median+sitesetup.diagplot__zeropoint_ZPrange[1] if \
             (sitesetup.diagplot__zeropoint_ZPrange[1] > 0) else zp_median+5*zp_std+0.3

    # Determine the min and max sdss magnitudes
    sdss_min = numpy.min(details['odi_sdss_matched'][:,details['photref_col_mag']] - 
                         details['odi_sdss_matched'][:,details['photref_col_err']]) - 1
    sdss_max = numpy.max(details['odi_sdss_matched'][:,details['photref_col_mag']] +
                         details['odi_sdss_matched'][:,details['photref_col_err']])
    # now round them to the nearest integer
    sdss_minint = sitesetup.diagplot__zeropoint_magrange[0] if \
                  (sitesetup.diagplot__zeropoint_magrange[0] > 0) else int(math.floor(sdss_min))
    sdss_maxint = sitesetup.diagplot__zeropoint_magrange[1] if \
                  (sitesetup.diagplot__zeropoint_magrange[1] > 0) else int(math.ceil(sdss_max))

    # Prepare some helpers for horizontal lines
    # There has to be a more elegant way to do this
    x_values = numpy.linspace(sdss_minint, sdss_maxint, 10)
    y_values = numpy.ones(shape=x_values.shape)

    #
    # Draw a grey-shaded region outlining the 1-sigma range
    # Plot it as a horizontal bar histogram with a single bin.
    #
    #ax.barh(zp_median-zp_std, (sdss_maxint-sdss_minint), height=2*zp_std,
    # ax.barh(y=0, # CHECK
    #         bottom=zp_median, height=2*zp_std,
    #         left=sdss_minint, width=(sdss_maxint - sdss_minint),
    #         label=u"1$\sigma$ range",
    #         color="#a0a0a0", edgecolor='#a0a0a0')
    ax.add_patch(matplotlib.patches.Rectangle(
        xy=(sdss_minint, zp_median-zp_std),
        width=sdss_maxint-sdss_minint,
        height=2*zp_std,
        facecolor="#a0a0a0", edgecolor='#a0a0a0',
    ))

    #
    # Compute and plot a histogram showing the distribution of ZPs
    #

    # Prepare a histogram to illustrate the distribution of ZP values
    binwidth = numpy.max([(0.2 * zp_std), 0.025])
    nbins = int(math.ceil((zp_max - zp_min) / binwidth))
    small_errors = details['odi_sdss_matched_smallerrors']
    zp_good = small_errors[:, details['photref_col_mag']] \
              - small_errors[:, details['odi_col_mag']]
    count, edges = numpy.histogram(zp_good, bins=nbins, range=[zp_min, zp_max], density=True)
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
    # Count the number of sources in each of the categories
    #
    n_ref = numpy.sum(use_for_calibration) # details['odi_sdss_matched_ref'].shape[0]
    n_outlier = details['odi_sdss_matched_outlier'].shape[0]
    n_all = details['odi_sdss_matched'].shape[0]
    n_ignored = n_all-(n_ref+n_outlier)


    #
    # Plot the actual measurments and values deemed to be outliers
    #
    # full_cat = details['odi_sdss_matched']
    # full_sdss_mag = full_cat[:, details['photref_col_mag']]
    # full_odi_mag = full_cat[:, details['odi_col_mag']]
    # full_sdss_err = full_cat[:, details['photref_col_err']]
    # full_odi_err = full_cat[:, details['odi_col_err']]
    # full_zp = full_sdss_mag - full_odi_mag
    # within_std = numpy.fabs(full_zp - zp_median) <= zp_std
    # outside_std = numpy.fabs(full_zp - zp_median) > zp_std


    #
    # Plot sources with large (i.e.exceeding the limit specified
    # in the sitesetup configuration) errors
    #
    # largeerr_cat = details['odi_sdss_matched_largeerrors']
    # largeerr_sdss_mag = largeerr_cat[:, details['photref_col_mag']]
    # largeerr_odi_mag = largeerr_cat[:, details['odi_col_mag']]
    # largeerr_sdss_err = largeerr_cat[:, details['photref_col_err']]
    # largeerr_odi_err = largeerr_cat[:, details['odi_col_err']]
    # largeerr_zp = largeerr_sdss_mag - largeerr_odi_mag
    # ax.errorbar(largeerr_sdss_mag,
    #             largeerr_zp,
    #             xerr=largeerr_sdss_err,
    #             yerr=numpy.hypot(largeerr_sdss_err, largeerr_odi_err),
    #             capsize=0,
    #             fmt='.', ms=0, color='#606060',
    #             label='ignored (#: %d)' % n_ignored)


    # smallerr_cat = details['odi_sdss_matched_smallerrors']
    # smallerr_sdss_mag = smallerr_cat[:, details['photref_col_mag']]
    # smallerr_odi_mag = smallerr_cat[:, details['odi_col_mag']]
    # smallerr_sdss_err = smallerr_cat[:, details['photref_col_err']]
    # smallerr_odi_err = smallerr_cat[:, details['odi_col_err']]
    # smallerr_zp = smallerr_sdss_mag - smallerr_odi_mag

    #
    # And for completeness, also plot sources deemed outliers
    #
    outlier_cat = details['odi_sdss_matched_outlier']
    ax.errorbar(outlier_cat[:, details['photref_col_mag']],
                outlier_cat[:, details['photref_col_mag']]-outlier_cat[:, details['odi_col_mag']],
                xerr=outlier_cat[:, details['photref_col_err']],
                yerr=numpy.hypot(outlier_cat[:, details['photref_col_err']],
                                 outlier_cat[:, details['odi_col_err']]),
                capsize=0,
                fmt="+", ms=5, color='red', linewidth=1,
                label='outlier (#: %d)' % n_outlier)

    #
    # Plot all valid calibrator sources
    #
    ref_cat = details['odi_sdss_matched_ref']
    ax.errorbar(ref_cat[:, details['photref_col_mag']],
                ref_cat[:, details['photref_col_mag']]-ref_cat[:, details['odi_col_mag']],
                xerr=ref_cat[:, details['photref_col_err']],
                yerr=numpy.hypot(ref_cat[:, details['photref_col_err']],
                                 ref_cat[:, details['odi_col_err']]),
                capsize=0,
                fmt="+", ms=5, color='blue', linewidth=1,
                label='valid calibrator (#: %d)' % n_ref)

    #
    # Overplot a white line to illustrate the median value
    #
    matplotlib.pyplot.plot(x_values+1, y_values*zp_median, linewidth=1, ls='-',
                           color='white', zorder=99)
    # ax.hline(y=zp_median, color='white', zorder=99)

    #
    # Add some additional data about the restricted fitting range and
    # the slope of phot. ZP with SDSS magnitude
    #
    if (details is not None):
        # Add some information about the fitted slope
        if (details['zp_magnitude_slope'] is not None):
            fit, err = details['zp_magnitude_slope']
            full_cat = details['odi_sdss_matched']
            minx = numpy.min(full_cat[:, details['photref_col_mag']])
            maxx = numpy.max(full_cat[:, details['photref_col_mag']])
            slopefit_x = numpy.linspace(minx-0.1*(maxx-minx), maxx+0.1*(maxx-minx), 100)
            slopefit_y = fit[0] + fit[1] * slopefit_x
            # ax.plot(slopefit_x, slopefit_y, "k-", label="fit ZP(%s-ODI)" % (details['catalog']))

    ax.grid(True)
    # ax.legend(loc='upper left', borderaxespad=0.5, prop={'size':9})

    rms = details['rms']
    sem = details['sem']
    photzp_text = u"ZP = %.3f \u00b1 %.3f mag\nr.m.s/s.e.m. = %.3f/%.3f mag" % (
        zp_median, zp_std, rms, sem,
    )
    ax.text(0.96, 0.05, photzp_text, fontsize=15,
            horizontalalignment='right',
            verticalalignment='bottom',
            transform=ax.transAxes,
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', linewidth=0))

    title_string = title if (title != None) else ""
    matplotlib.pyplot.title(title_string)
    ref_mag_name = "reference"
    ref_filter = "none"
    if (details is not None and 'catalog' in details):
        ref_mag_name = details['catalog']
        ref_filter = details['reference_filter']

    matplotlib.pyplot.xlabel("%s magnitude in %s" % (ref_mag_name, ref_filter), labelpad=7)
    matplotlib.pyplot.ylabel("zeropoint (%s [%s] - odi [%s])" % (ref_mag_name, ref_filter, odi_filtername), labelpad=20)
    matplotlib.pyplot.xlim((sdss_minint, sdss_maxint))
    matplotlib.pyplot.ylim((zp_min, zp_max))
#    matplotlib.pyplot.axes().set_aspect('equal')

    if (output_filename is None):
        fig.show()
    else:
        fig.set_size_inches(8,6)
        if (options is not None):
            for ext in options['plotformat']:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext, bbox_inches='tight')
                    logger.info("plot saved as %s.%s" % (output_filename, ext))
        else:
            fig.savefig(output_filename+".png", dpi=100,bbox_inches='tight')

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()
    logger.debug("done!")
























#####################################################################################
#
#
# Following are the routines to create the zeropoint maps
#
#
#####################################################################################

def plot_zeropoint_map(details, ota_select, ota_outlines, output_filename, zp_range, options=None, plot_formats=None):
    """

    Plot the photometric zeropoint, color-coded, as function of position

    """

    logger = logging.getLogger("DiagPlot_ZPmap")

    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    ax.ticklabel_format(useOffset=False)
    zp_min, zp_max = zp_range

    full_catalog = details['odi_sdss_matched_smallerrors'].copy()
    if (ota_select is not None):
        in_this_ota = full_catalog[:, SXcolumn['ota']+2] == ota_select
        full_catalog = full_catalog[in_this_ota]

    ra = full_catalog[:, 0]
    dec = full_catalog[:, 1]
    zp = full_catalog[:, details['photref_col_mag']] - full_catalog[:, details['odi_col_mag']]

    around_zero = False
    min_ra, max_ra = numpy.min(ra), numpy.max(ra)
    if ((max_ra - min_ra) > 180 or (max_ra > 0 and min_ra < 0) or max_ra > 360):
        # This means we most likely have to deal with coordinates around 0
        ra[ra > 180] -= 360.
        around_zero = True

    # if (ra.shape[0] < 5):
    #     zp_clipped = zp
    #     clipped = numpy.isfinite(zp)
    # else:
    #     zp_clipped, clipped = three_sigma_clip(zp, return_mask=True)

    zp_median_ota = numpy.median(zp)
    zp_std_ota = numpy.std(zp)

    #print ra_ota
    #print zp_ota

    colormap_name = "nipy_spectral"
    sc = matplotlib.pyplot.scatter(ra, dec, c=zp, alpha=0.75, 
                                   vmin=zp_min, vmax=zp_max, edgecolor='none',
                                   s=9, cmap=matplotlib.pyplot.cm.get_cmap(colormap_name))
        #matplotlib.pyplot.colorbar(sc, cm)

    if (ota_outlines is not None):
        corners = numpy.array(ota_outlines)
        if (around_zero):
            corners[:,:,0][corners[:,:,0] > 180.] -= 360.
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
        ax.add_collection(coll)

    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    ax.autoscale_view()
    # Invert the x-axis to have north up and east to the left
    ax.set_xlim(ax.get_xlim()[::-1])

    colorbar = matplotlib.pyplot.colorbar(sc) #cmap=matplotlib.pyplot.cm.get_cmap(colormap_name))
    colorbar.set_label("phot. zeropoint")

    if (output_filename is None):
        fig.show()
    else:
        fig.set_size_inches(8,6)
        if (plot_formats is not None):
            for ext in plot_formats:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext, dpi=100, bbox_inches='tight')
                    logger.debug("saving plot to %s.%s" % (output_filename, ext))
        else:
            fig.savefig(output_filename+".png", dpi=100)

    matplotlib.pyplot.close()
    return


def photocalib_zeropoint_map(details,
                             output_filename,
                             sdss_filtername, odi_filtername,
                             title=None,
                             ota_outlines=None,
                             options=None,
                             also_plot_singleOTAs=True,
                             allow_parallel=True
                             ):

    """

    Prepare all necessary data and generate all photometric zeropoint maps.

    """
    logger = logging.getLogger("DiagPlot_ZPmap")
    logger.info("Creating the photometric calibration per OTA plot ...")

    processes = []

    # print(details)

    zp_median = details['median'] 
    zp_raw = details['odi_sdss_matched'][:, details['photref_col_mag']] \
             - details['odi_sdss_matched'][:, details['odi_col_mag']]
    zp_min = scipy.stats.scoreatpercentile(zp_raw,  5) if \
             (sitesetup.diagplot__zeropointmap_range[0] <= 0) else \
             zp_median - sitesetup.diagplot__zeropointmap_range[0]
    zp_max = scipy.stats.scoreatpercentile(zp_raw, 95) if \
             (sitesetup.diagplot__zeropointmap_range[1] <= 0) else \
             zp_median + sitesetup.diagplot__zeropointmap_range[1]

    # Create one plot for the full focal plane, using boxes to outlines OTAs
    zp_range = (zp_min, zp_max)
    kwargs = {"details": details,
              "ota_select": None,
              "ota_outlines": ota_outlines,
              "output_filename": output_filename,
              #"options": options,
              "zp_range": zp_range,
              'plot_formats': options['plotformat']
    }
    # if (not allow_parallel):
    plot_zeropoint_map(**kwargs) #details, None, ota_outlines, output_filename, options, zp_range)
    # else:
    #     # print(kwargs)
    #     p = multiprocessing.Process(target=plot_zeropoint_map,
    #                                 kwargs=kwargs)
    #     p.daemon = True
    #     p.start()
    #     processes.append(p)

#    return

    logger.debug(ota_outlines)

    # If requested, do the same for the individual OTAs
    if (also_plot_singleOTAs):
        list_of_otas = set(details['odi_sdss_matched_smallerrors'][:, SXcolumn['ota']+2])
        for this_ota in list_of_otas:
            in_this_ota = (details['odi_sdss_matched_smallerrors'][:, SXcolumn['ota']+2] == this_ota)
            if (numpy.sum(in_this_ota) <= 0):
                continue

            ota_plotfile = create_qa_otaplot_filename(output_filename, this_ota, options['structure_qa_subdirs'])


            kwargs = {"details": details, 
                      "ota_select": this_ota,
                      "ota_outlines": None,
                      "output_filename": ota_plotfile,
                      #"options": options,
                      "zp_range": zp_range,
                      'plot_formats': options['plotformat']
            }
            # if (not allow_parallel or True):
            plot_zeropoint_map(**kwargs) #details, ota_select=this_ota, ota_outlines, output_filename, options, zp_range)
                #ra_ota, dec_ota, zp_ota, None, ota_plotfile, options, zp_range)

    for p in processes:
        p.join()

    logger.debug("done!")

    # import time
    # time.sleep(2)
    # print("done with zpmap")
    # time.sleep(2)

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
                     plot_formats=None,
                     ):
    """

    Generate a PSF size map, with points representing sources. The color of the 
    point encodes the image quality / PSF FWHM.

    """
    logger = logging.getLogger("DiagPlot_PSFSize")
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    ax.ticklabel_format(useOffset=False)

    fwhm_min, fwhm_max = fwhm_range

    colormap_name = "nipy_spectral" #"hsv", "jet" #"rainbow" #"Dark2"

    around_zero = False
    if ((numpy.max(ra) - numpy.min(ra)) > 180):
        # This means we most likely have to deal with coordinates around 0
        ra[ra > 180] -= 360.
        around_zero = True

    sc = matplotlib.pyplot.scatter(ra, dec, c=fwhm, alpha=0.75, 
                                   vmin=fwhm_min, vmax=fwhm_max, edgecolor='none',
                                   s=9, cmap=matplotlib.pyplot.cm.get_cmap(colormap_name))

    if (ota_outlines is not None):
        corners = numpy.array(ota_outlines)
        if (around_zero):
            corners[:,:,0][corners[:,:,0] > 180.] -= 360.

        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', 
                                                     linestyle='-', linewidth=0.5)
        ax.add_collection(coll)

    ax.set_title(title)
    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    ax.autoscale_view()
    # Invert the x-axis to have north up and east to the left
    ax.set_xlim(ax.get_xlim()[::-1])

    colorbar = matplotlib.pyplot.colorbar(sc) #cmap=matplotlib.pyplot.cm.get_cmap(colormap_name)) #spectral
    colorbar.set_label("FWHM [arcsec]")

    if (output_filename is None):
        fig.show()
    else:
        fig.set_size_inches(8,6)
        if (plot_formats is not None):
            for ext in plot_formats:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext, dpi=100,bbox_inches='tight')
                    logger.debug("saving plot to %s.%s" % (output_filename, ext))
        else:
            fig.savefig(output_filename+".png", dpi=100)

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()
  

  
def diagplot_psfsize_map(ra, dec, fwhm, ota, output_filename,
                         title=None,
                         ota_outlines=None,
                         options=None,
                         also_plot_singleOTAs=True,
                         title_info=None,
                         ):

    """

    Generate all PSF FWHM plots for the entire flocalplane and, optioanlly, for
    each OTA.

    """

    logger = logging.getLogger("DiagPlot_PSFSize")
    logger.info("Creating the PSF size as heatmap plot ...")

    fwhm_scale_arcsec = 3600.0 

    fwhm *= fwhm_scale_arcsec
    fwhm_min = scipy.stats.scoreatpercentile(fwhm,  5)
    fwhm_max = scipy.stats.scoreatpercentile(fwhm, 95)
    fwhm_clipped, clipmask = three_sigma_clip(fwhm, return_mask=True)

    fwhm_range = (fwhm_min, fwhm_max)
    fwhm_range = (0.4, 2.5)

    processes = []

    title = "PSF map"
    if (title_info is not None):
        try:
            title = "FWHM map - %(OBSID)s (focal plane)\n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec)" % title_info
        except:
            pass

    # Create one plot for the full focal plane, using boxes to outlines OTAs
    plot_args = {"ra": ra, 
                 "dec": dec, 
                 "fwhm": fwhm, 
                 "output_filename": output_filename, 
                 "fwhm_range": fwhm_range,
                 "title": title,
                 "ota_outlines": ota_outlines, 
                 #"options": options,
                 "plot_formats": options['plotformat']

    }
    plot_psfsize_map(**plot_args)
    # p = multiprocessing.Process(target=plot_psfsize_map, kwargs=plot_args)
    # p.daemon = True
    # p.start()
    # processes.append(p)
   
    # plot_psfsize_map(ra, dec, fwhm, output_filename, 
    #                  fwhm_range=fwhm_range,
    #                  title=title,
    #                  ota_outlines=ota_outlines,
    #                  options=options)

    logger.debug("OTA-outlines:")
    logger.debug(ota_outlines)

    # If requested, do the same for the individual OTAs
    if (also_plot_singleOTAs):
        list_of_otas = set(ota)
        for this_ota in list_of_otas:

            in_this_ota = (ota == this_ota)
            if (numpy.sum(in_this_ota) <= 0):
                continue
            fwhm_ota = fwhm[in_this_ota]
            ra_ota = ra[in_this_ota]
            dec_ota = dec[in_this_ota]

            title = "PSF map, OTA %02d" % (this_ota)
            if (title_info is not None):
                title_info['OTA'] = this_ota
                try:
                    title = "FWHM map - %(OBSID)s OTA %(OTA)02d\n(%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec)" % title_info
                except:
                    pass

            ota_plotfile = create_qa_otaplot_filename(output_filename, this_ota, options['structure_qa_subdirs'])

            plot_args = {"ra": ra_ota, 
                         "dec": dec_ota, 
                         "fwhm": fwhm_ota, 
                         "output_filename": ota_plotfile, 
                         "fwhm_range": fwhm_range,
                         "title": title,
                         "ota_outlines": None,
                         #"options": options,
                         "plot_formats": options['plotformat']
            }
            # print(plot_args)
            # _obj = _ForkingPickler.dumps(plot_args)
            # print(_obj)
            #
            # pickle.dumps(plot_args)
            # # continue
            # p = multiprocessing.Process(target=plot_psfsize_map, kwargs=dict(plot_args))
            # p.daemon = True
            # p.start()
            # processes.append(p)

            plot_psfsize_map(**plot_args)

            #time.sleep(1)


    for p in processes:
        p.join()

    logger.debug("done!")

    return




from matplotlib.collections import LineCollection
elongation_limit = 3.5


def plot_psfshape_map(ra, dec, elongation, angle, fwhm, 
                      output_filename='psfshape_test',
                      options=None,
                      title='test',
                      ota_outlines=None,
                      show_round_stars=True,
                      plot_formats=None
                      ):

    logger = logging.getLogger("DiagPlot_PSFShape")
    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.ticklabel_format(useOffset=False)

    around_zero = False
    if ((numpy.max(ra) - numpy.min(ra)) > 180):
        
        # This means we most likely have to deal with coordinates around 0
        ra[ra > 180] -= 360.
        around_zero = True

    ramin, ramax = numpy.min(ra), numpy.max(ra)
    decmin, decmax = numpy.min(dec), numpy.max(dec)

    dimension = numpy.min([ramax-ramin, decmax-decmin])

    # Convert axis ratio to axis excess (i.e. round stars have then 0 elongation)
    elongation -= 1.0

    #
    # Plot all round stars as small grey datapoints, mostly for reference 
    # and orientation
    #
    round_stars = elongation < 0.25
    elongated = elongation > 0.1
    logger.debug("Plotting round stars")
    if (numpy.sum(round_stars) > 0 and show_round_stars):
        ax.plot(ra[round_stars], dec[round_stars],
                color='grey', marker='o', linestyle='None',
                markeredgecolor='none', markersize=3, alpha=0.5
            )


    #
    # Compute the arrow positions
    #
    # scaling is 10x, i.e. stars with 1:2 axis ration (elongation=1) 
    # get drawn with 10'' long streaks one either side of the star
    #

    scale = 10./3600. 
    dx = numpy.cos(numpy.radians(angle))*elongation * scale / numpy.cos(numpy.radians(dec))
    dy = numpy.sin(numpy.radians(angle))*elongation * scale
    # print dx[:25]
    # print dy[:25]

    # convert FWHM from degrees to arcseconds
    # fwhm *= 3600.
    # only use stars with reasonable seeing - this is for color-mapping only
    seeing_range = [0.4, 2.5]
    rel_seeing = fwhm - seeing_range[0] / (seeing_range[1] - seeing_range[0])
    arrows = []
    for i in range(ra.shape[0]):
        if (not elongated[i] or fwhm[i] < seeing_range[0]):
            continue
        position = [[ra[i]-dx[i], dec[i]-dy[i]],
                    [ra[i]+dx[i], dec[i]+dy[i]]]
        arrows.append(position)

    # print elongated[:25]
    # print arrows[:25]
    # print arrows[elongated]

    # Draw all streaks
    if (numpy.sum(elongated) > 0):
        logger.debug("Plotting elongated stars")
        colormap = matplotlib.pyplot.cm.get_cmap("nipy_spectral")
        lines = LineCollection(arrows,
                           alpha=0.8, 
                           array=fwhm[elongated],
                           cmap=colormap, 
                           norm=matplotlib.pyplot.Normalize(0.5,2.5),
                           zorder=99)
        ax.add_collection(lines)

    # Determine min and max values

    ax.set_title(title)
    ax_dx = 0.1*(ramax - ramin)
    ax_dy = 0.1*(decmax - decmin)
    ax.set_xlim((ramin-ax_dx, ramax+ax_dx))
    ax.set_ylim((decmin-ax_dy, decmax+ax_dy))
    ax.set_xlabel("RA [degrees]")
    ax.set_ylabel("DEC [degrees]")
    # Invert the x-axis to have north up and east to the left
    ax.set_xlim(ax.get_xlim()[::-1])

    # Draw OTA outlines if requested
    if (ota_outlines is not None):
        corners = numpy.array(ota_outlines)
        if (around_zero):
            corners[:,:,0][corners[:,:,0] > 180.] -= 360.
        coll = matplotlib.collections.PolyCollection(corners,facecolor='none',edgecolor='#808080', linestyle='-')
        ax.add_collection(coll)

    # Save plot to file
    if (output_filename is None):
        matplotlib.pyplot.show()
        fig.show()
    else:
        fig.set_size_inches(8,6)
        if (plot_formats is not None):
            for ext in plot_formats:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext, dpi=100,bbox_inches='tight')
                    logger.debug("saving plot to %s.%s" % (output_filename, ext))
        else:
            fig.savefig(output_filename+".png", dpi=100,bbox_inches='tight')

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()
    logger.debug("done!")

    return


def  diagplot_psfshape_map(ra, dec, elongation, angle, fwhm, ota, 
                           output_filename='psfshape_test',
                           options=None,
                           title='test',
                           title_info=None,
                           ota_outlines = None,
                           also_plot_singleOTAs=False,
                           plot_formats=None,
):



    """

    Generate all PSF Shape plots for the entire focalplane and, optionally, for
    each OTA.

    """

    logger = logging.getLogger("DiagPlot_PSFShape")
    logger.info("Creating the PSF-shape diagnostic plot ...")

    processes = []

    title = "PSF Shape"
    if (title_info is not None):
        try:
            title = "PSF Shape - %(OBSID)s (focal plane)\n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec)" % title_info
        except:
            pass

    elongation[elongation > elongation_limit] = elongation_limit

    # Create one plot for the full focal plane, using boxes to outlines OTAs
    plot_args = {"ra": ra, 
                 "dec": dec, 
                 "elongation": elongation,
                 "angle": angle,
                 "fwhm": fwhm, 
                 "output_filename": output_filename, 
                 "title": title,
                 "ota_outlines": ota_outlines, 
                 #"options": options,
                 "show_round_stars": False,
                 "plot_formats": options['plotformat']
    }
    plot_psfshape_map(**plot_args)
    # p = multiprocessing.Process(target=plot_psfshape_map, kwargs=plot_args)
    # p.start()
    # processes.append(p)
   
    # plot_psfsize_map(ra, dec, fwhm, output_filename, 
    #                  fwhm_range=fwhm_range,
    #                  title=title,
    #                  ota_outlines=ota_outlines,
    #                  options=options)

    logger.debug("OTA-outlines:")
    logger.debug(ota_outlines)

    # If requested, do the same for the individual OTAs
    if (also_plot_singleOTAs):
        list_of_otas = set(ota)
        for this_ota in list_of_otas:
            in_this_ota = (ota == this_ota)
            if (numpy.sum(in_this_ota) <= 0):
                continue

            ra_ota = ra[in_this_ota]
            dec_ota = dec[in_this_ota]
            elongation_ota = elongation[in_this_ota]
            angle_ota = angle[in_this_ota]
            fwhm_ota = fwhm[in_this_ota]
            title = "PSF map, OTA %02d" % (this_ota)

            if (title_info is not None):
                title_info['OTA'] = this_ota
                try:
                    title = "PSF Shape - %(OBSID)s OTA %(OTA)02d\n%(OBJECT)s - %(FILTER)s - %(EXPTIME)dsec)" % title_info
                except:
                    pass

            ota_plotfile = create_qa_otaplot_filename(output_filename, this_ota, options['structure_qa_subdirs'])

            plot_args = {"ra": ra_ota, 
                         "dec": dec_ota, 
                         "elongation": elongation_ota,
                         "angle": angle_ota,
                         "fwhm": fwhm_ota,
                         "output_filename": ota_plotfile, 
                         "title": title,
                         "ota_outlines": ota_outlines, 
                         #"ota_outlines": None,
                         #"options": options,
                         "plot_formats": options['plotformat']
            }
            # p = multiprocessing.Process(target=plot_psfshape_map, kwargs=plot_args)
            # p.start()
            # processes.append(p)
            plot_psfshape_map(**plot_args)

    for p in processes:
        p.join()

    logger.debug("done!")

    return


def diagplot_photflat(extnames, data, one_sigma=None,
                      title=None,
                      output_filename="photflat",
                      options=None,
                      n_sigma=3,
                      force_symmetric=False,
                      plot_formats=None,
                      ):

    logger = logging.getLogger("DiagPlot_PhotFlat")

    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    px_y, px_x = numpy.indices((9,9), dtype=float)
    px_x *= 512.
    px_y *= 512.
    ota_px_size = 1
    rel_ota_size = 0.95

    fluxfac = data[:, :-1, :-1] #numpy.power(10., 0.4*data[:, :-1, :-1])

    if (one_sigma is None or not force_symmetric):
        user_sigma_plus = 50. + 50.*scipy.special.erf(float(n_sigma)/numpy.sqrt(2))
        user_sigma_minus = 100. - user_sigma_plus
        print(user_sigma_plus, user_sigma_minus)

        stats = numpy.nanpercentile(fluxfac, [16,84,50,5,95, user_sigma_minus, user_sigma_plus])
        print(stats)
        one_sigma = 0.5*(stats[1]-stats[0])
        print(one_sigma, 0.25*(stats[4]-stats[3]))
        median = stats[2]
        print(median)
        one_sigma += numpy.fabs(median-1.)
        logger.info("Using automatic intensity scaling")
    if (force_symmetric):
        logger.info("Photflat-scaling: 1.0 +/- %.3f" % (one_sigma))
        vmin = 1. - n_sigma*one_sigma
        vmax = 1. + n_sigma*one_sigma
    else:
        vmin = stats[-2]
        vmax = stats[-1]
        logger.info("photflat-scaling: %f .. %f" % (vmin, vmax))


    for ota, extname in enumerate(extnames): #range(data.shape[0]):
        #ox = ota_x[ota]
        #oy = ota_y[ota]
        ox,oy = int(extname[3]), int(extname[4])
        # print ox, oy, data[ota]
        # fluxdata = numpy.power(10., 0.4*data[ota,:-1,:-1])
        fluxdata = fluxfac[ota, :, :] #numpy.power(10., 0.4*data[ota,:-1,:-1])

        this_x = px_x + ox*ota_px_size
        this_y = px_y + oy*ota_px_size

        # print this_x

        data_prepped = numpy.flip(numpy.ndarray.transpose(fluxdata), 0)

        ims = ax.imshow(data_prepped,
            #aspect='equal',
            interpolation='bilinear',
            extent=((ox-0.5*rel_ota_size)*ota_px_size, (ox+0.5*rel_ota_size)*ota_px_size,
                    (oy-0.5*rel_ota_size)*ota_px_size, (oy+0.5*rel_ota_size)*ota_px_size),
            cmap=matplotlib.pyplot.cm.get_cmap('nipy_spectral'), #matplotlib.cm.gray,
            vmin=vmin, vmax=vmax,
        )

    cbar = fig.colorbar(ims, orientation='vertical')
    print("intensity range: %.3f ... %.3f" % (vmin, vmax))

    if (title is None):
        title = "Photometric flatfield"
    ax.set_title(title)

    ax.set_xlim((0.8-0.5*rel_ota_size,5.2+0.5*rel_ota_size))
    ax.set_ylim((0.8-0.5*rel_ota_size,6.2+0.5*rel_ota_size))

    if (output_filename is None):
        matplotlib.pyplot.show()
        fig.show()
    else:
        fig.set_size_inches(8,6)
        if (plot_formats is not None):
            for ext in plot_formats:
                if (ext != ''):
                    fig.savefig(output_filename+"."+ext, dpi=100,bbox_inches='tight')
                    logger.debug("saving plot to %s.%s" % (output_filename, ext))
        else:
            fig.savefig(output_filename+".png", dpi=100,bbox_inches='tight')



def crossout(x, y, size):
    # [
    #     [ota_x, ota_y],
    #     [ota_x + ota_cellsize, ota_y],
    #     [ota_x + ota_cellsize, ota_y + ota_cellsize],
    #     [ota_x, ota_y + ota_cellsize],
    #     [ota_x, ota_y],
    #     [ota_x + ota_cellsize, ota_y + ota_cellsize],
    #     [ota_x + ota_cellsize, ota_y],
    #     [ota_x, ota_y + ota_cellsize],
    # ])

    coords = [
        [x, y],
        [x + size, y],
        [x + size, y + size],
        [x, y + size],
        [x, y],
        [x + size, y + size],
        [None, None],
        [x + size, y],
        [x, y + size],
        [None, None],
    ]
    return coords

def plot_cellbycell_stats(
        hdulist,
        title=None,
        vmin=None, vmax=None,
        plotfile=None,
        stats=None,
        numberformat="%.3f",
        showlabels=True,
        units=None,
        crossout_missing=True,
        plotformat=None,
):
    """

    Function generating the WCS shift plot from the given data.

    """

    logger = logging.getLogger("DiagPlot_CellbyCellStats")

    import matplotlib.pyplot
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)
    ax.ticklabel_format(useOffset=False)

    if (title is not None):
        ax.set_title(title, fontsize=10.)

    if (plotformat is None):
        plotformat = ['pdf']

    fpl = podi_focalplanelayout.FocalPlaneLayout(hdulist)

    if (stats is None):
        stats = [('', numpy.nanmean),
                 ('', numpy.std),
                 ]

    # prepare boxes and text
    cellsize = 0.12
    ota_cellsize = 8*cellsize
    all_corners = []
    all_colors = []
    all_labels = []
    all_label_coords = []
    # get stats by ota

    lines_to_plot = []

    for ota_x, ota_y in itertools.product(range(8), repeat=2):
        extname = "OTA%d%d.SCI" % (ota_x, ota_y)

        try:
            data = hdulist[extname].data
        except:
            logger.debug("There's no OTA %s" % (extname))
            lines_to_plot.append(crossout(ota_x, ota_y, ota_cellsize))
            continue
        logger.debug("Working on cell-by-cell stats from OTA %s" % (extname))

        for cell_x, cell_y in itertools.product(range(8), repeat=2):

            # figure out if the cell is one of the broken ones
            cell = (cell_x, cell_y)
            if (cell in fpl.broken_cells[hdulist[extname].header['OTA_ID']]):
                lines_to_plot.append(crossout(ota_x+cell_x*cellsize, ota_y+(7-cell_y)*cellsize, cellsize))
                continue

            labels = []
            colorvalue = None
            # all is good, compute statistics
            for statname, statfct in stats:
                # print statname, statfct

                cellarea = cell2ota__get_target_region(
                    x=cell_x,
                    y=cell_y,
                    binning=fpl.get_hardware_binning(),
                    trimcell=15)
                x1, x2, y1, y2 = cellarea
                value = statfct(data[y1:y2, x1:x2])

                if (statname != ''):
                    label = "%s: %s" % (statname, numberformat % (value))
                else:
                    label = numberformat % value
                labels.append(label)

                if (colorvalue is None):
                    colorvalue = value

            # save the color value for this cell
            all_colors.append(colorvalue)

            # compute the corners of each cell so we can plot each cell as little rectangular polygon later
            x1 = ota_x + cell_x * cellsize
            y1 = ota_y + (7 - cell_y) * cellsize
            corners = [[x1, y1], [x1 + cellsize, y1], [x1 + cellsize, y1 + cellsize], [x1, y1 + cellsize]]
            all_corners.append(corners)

            # save what text to put where
            all_labels.append("\n".join(labels))
            all_label_coords.append([x1, y1])

    if (crossout_missing):
        all_lines = numpy.array(lines_to_plot)
        # print all_lines
        # print all_lines.shape
        # print all_lines[0, :, :]

        ax.plot(all_lines[:, :, 0].flatten(),
                all_lines[:, :, 1].flatten(),
                c='#bbbbbb',
                linewidth=0.2)

    #
    #
    #

    cmap = matplotlib.pyplot.cm.get_cmap('nipy_spectral')

    corners = numpy.array(all_corners)
    ax.set_xlim([0, 8])
    ax.set_ylim([0, 8])

    nl_min = numpy.min(all_colors[numpy.isfinite(all_colors)]) if vmin is None else vmin
    nl_max = numpy.max(all_colors[numpy.isfinite(all_colors)]) if vmax is None else vmax

    colorvalues = cmap(
        (numpy.array(all_colors) - nl_min) / (nl_max - nl_min))  # [matplotlib.pyplot.cm.jet(x) for x in all_intensity]
    norm_colors = (numpy.array(all_colors) - nl_min) / (nl_max - nl_min)

    coll = matplotlib.collections.PolyCollection(corners,  # facecolor='#505050', #
                                                 facecolor=colorvalues,
                                                 edgecolor='black', linestyle='-', linewidth=0.1,
                                                 cmap=matplotlib.pyplot.cm.get_cmap('nipy_spectral'),
                                                 )

    img = matplotlib.pyplot.imshow([[1e9], [1e9]], vmin=nl_min, vmax=nl_max, cmap=cmap, extent=(0, 0, 0, 0),
                                   origin='lower')

    if (units is None):
        fig.colorbar(img, format=numberformat)
    else:
        fig.colorbar(img, format=numberformat, label=units)

    #
    # also show the values as text labels in each cell
    #
    if (showlabels):
        for cell in range(len(all_labels)):
            x, y = all_label_coords[cell]
            label = all_labels[cell]
            color = "white" if norm_colors[cell] < 0.15 else 'black'
            ax.text(x+0.5*cellsize, y+0.5*cellsize, label,
                    ha='center', va='center', color=color,
                    fontsize=1.1)
            # print x,y,label

    ax.set_xlim(-0.1, 8.1)
    ax.set_ylim(-0.1, 8.1)
    ax.add_collection(coll)

    logger.debug("Saving plot as %s" % (plotfile))
    fig.set_size_inches(8, 6)
    for ft in plotformat:
        plot_fn = "%s.%s" % (plotfile, ft)
        try:
            fig.savefig(plot_fn, dpi=300, bbox_inches='tight')
        except (IOError, ValueError) as err:
            logger.error("Unable to create %s: %s" % (plot_fn, str(err)))

    matplotlib.pyplot.close(fig)
    matplotlib.pyplot.close()
    logger.debug("done!")

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
        print(fwhm)

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

    elif (cmdline_arg_isset("-zeropoint_map")):
        fitsfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(fitsfile)
        filtername = hdulist[0].header['FILTER']
        print(filtername)

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


    elif (cmdline_arg_isset("-wcstest")):
        fitsfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(fitsfile)
        
        matched = hdulist['WCSCAL.CAT']
        
        table = matched.data
        n_src = matched.header['NAXIS2']

        ra = table['ODI_RA'].reshape((n_src,1))
        dec = table['ODI_DEC'].reshape((n_src,1))
        matched_radec_odi = numpy.append(ra, dec, axis=1)
        print(matched_radec_odi.shape)

        ra = table['REF_RA'].reshape((n_src,1))
        dec = table['REF_DEC'].reshape((n_src,1))
        matched_radec_2mass = numpy.append(ra, dec, axis=1)
 
        matched_ota = table['OTA']

        import podi_collectcells
        options = podi_collectcells.read_options_from_commandline()

        wcsdiag_scatter(matched_radec_odi, 
                        matched_radec_2mass, 
                        matched_ota,
                        "wcs_test", options=options, ota_wcs_stats=None,
                        also_plot_singleOTAs=True,
                        title_info = None,
                    )

    elif (cmdline_arg_isset("-psfshape")):

        inputframe = get_clean_cmdline()[1]
        plotname = get_clean_cmdline()[2]

        print(inputframe,"-->",plotname)
        

        from podi_collectcells import read_options_from_commandline
        options = read_options_from_commandline(None)

        hdulist = pyfits.open(inputframe)
        hdulist.info()
        try:
            odi_cat = hdulist['CAT.ODI'].data
        except:
            print("no source catalog found")
            sys.exit(0)

        ota_outlines = derive_ota_outlines(hdulist)

        ota = odi_cat.field('OTA')
        flags = odi_cat.field('FLAGS')

        in_ota = (flags == 0) #& (ota == 32)

        ra = odi_cat.field('RA')[in_ota]
        dec = odi_cat.field('DEC')[in_ota]
        elongation = odi_cat.field('ELONGATION')[in_ota]
        angle = odi_cat.field('THETAWIN_IMAGE')[in_ota]
        ota = odi_cat.field('OTA')[in_ota]
        fwhm = odi_cat.field('FWHM_WORLD')[in_ota]

        diagplot_psfshape_map(ra, dec, elongation, angle, fwhm, ota, 
                              output_filename=plotname,
                              options=options,
                              also_plot_singleOTAs=True,
                              
                              ota_outlines=ota_outlines,
                      )

    else:
        matched_cat = numpy.loadtxt("odi+2mass.matched")

        filename = "output"

    #    make_plots(matched_cat, filename)
        wcsdiag_scatter(matched_cat, filename+".wcs1")
        wcsdiag_shift(matched_cat, filename+".wcs2")
    #    wcsdiag_shift(matched_cat, None)







