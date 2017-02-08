#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits, os, sys
import numpy
import itertools

import matplotlib
import matplotlib.pyplot as plot
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

qr_dir = "/work/podi_prep56"
sys.path.append(qr_dir)

import podi_swarpstack
from podi_commandline import *
from podi_definitions import *
import podi_logging
import podi_sitesetup as sitesetup

from scipy.stats import gaussian_kde

limit_fwhm_max = 3.5


def get_guidephotom_filelist(directory, filebase):
    
    filelist = []
    for l,n in itertools.product(['A', 'B', 'C', 'F', 'G', 'H'], range(1,5)):
        
        filename = "%s/expVideo/%s.odi%s%d.photom.fits" % (
            directory, filebase, l, n)
        #print filename

        if (os.path.isfile(filename)):
            filelist.append(filename)

    return filelist




def draw_guidestarplot(filelist, title=None, plot_filename=None):

    #
    # Start plot
    #
    fig = plot.figure()
    fig.set_facecolor('#e0e0e0')
    fig.canvas.set_window_title(
        title if not title == None else "Guide Star details")
    fig.suptitle(
        title if not title == None else "Guide Star details")
    #
    # Set all plot dimensions
    #
    plot_width = 0.73
    ax_flux = plot.axes([0.10, 0.57, 0.73, 0.36])
    ax_fwhm = plot.axes([0.10, 0.20, 0.73, 0.36])

    hist_flux = plot.axes([0.84, 0.57, 0.14, 0.36])
    hist_fwhm = plot.axes([0.84, 0.20, 0.14, 0.36])

    # 
    # Set all other plot-specific properties
    #

    null_fmt = plot.NullFormatter()

    # Top panel
    ax_flux.grid(True)
    #ax_flux.set_xlabel("time [seconds]")
    ax_flux.set_ylabel("normalized flux")
    ax_flux.xaxis.set_major_formatter(null_fmt)
    ax_fwhm.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))
    
    # Bottom panel
    ax_fwhm.grid(True)
    ax_fwhm.set_xlabel("time [seconds]")
    ax_fwhm.set_ylabel("FWHM [arcsec]")
    ax_fwhm.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.1f"))

    hist_fwhm.xaxis.set_major_formatter(null_fmt)
    hist_fwhm.yaxis.set_major_formatter(null_fmt)
    hist_flux.xaxis.set_major_formatter(null_fmt)
    hist_flux.yaxis.set_major_formatter(null_fmt)

    # from matplotlib.patches import Rectangle
    # currentAxis = fig.canvas
    # currentAxis.add_patch(Rectangle((0, 0), 1, 0.10, facecolor="white", fill=True))

    global_min_flux = 1e9
    shutter_delay = 1.25
    pixelscale = 0.118

    max_time = -1.

    all_flux = []
    all_fwhm = []

    colors = ['blue', 
              'red', 
              'green',
              'purple', 
              'orange', 
              'dimgrey',
              'sienna', ]

    text_positions = numpy.array([
        [0.01, 0.05],
        [0.01, 0.01],
        [0.50, 0.05],
        [0.50, 0.01],
        [0.99, 0.05],
        [0.99, 0.01],
        ])
    text_align = ['left', 'left', 'center', 'center', 'right', 'right']

    guidestats = {}
    if (len(filelist) <= 0):
        return guidestats

    for idx, infile in enumerate(filelist):
        
        hdulist = pyfits.open(infile)

        data = hdulist['VideoPhotometry'].data
        timestamp = data['TimeStamp']
        flux = data["ApertureBasedFluxDN"]
        
        fwhm_x = data["FWHMX"]
        fwhm_y = data["FWHMY"]
        fwhm = numpy.hypot(fwhm_x, fwhm_y)

        # correct timestamp to seconds since exposure start
        timestamp = (timestamp - numpy.min(timestamp)) / 1000.
        
        shutter_open = (timestamp > shutter_delay) & \
                       (timestamp < (numpy.max(timestamp) - shutter_delay))

        flux = flux[shutter_open]
        fwhm = fwhm[shutter_open] * pixelscale
        timestamp = timestamp[shutter_open]

        min_flux, max_flux = numpy.min(flux), numpy.max(flux)
        global_min_flux = numpy.min([global_min_flux, min_flux/max_flux])
        max_time = numpy.max([max_time, numpy.max(timestamp)])

        min_fwhm, max_fwhm = numpy.min(fwhm), numpy.max(fwhm)

        #ax_fwhm.plot(timestamp, fwhm)
        ax_fwhm.scatter(timestamp, fwhm, 
                        linewidth=0, alpha=0.3,
                        c=colors[idx])

        #ax_flux.plot(timestamp, flux/max_flux)
        ax_flux.scatter(timestamp, flux/max_flux,
                        linewidth=0, alpha=0.3,
                        c=colors[idx])

        #
        # Plot some smooth curves through the data points
        #
        spline_t = numpy.linspace(timestamp[1], timestamp[-2], timestamp[-1]/5.)

        # smooth_flux = scipy.interpolate.UnivariateSpline(
        #     x=timestamp, y=flux/max_flux, w=None, bbox=[None, None], k=3, 
        #     s=0.3)
        smooth_flux = scipy.interpolate.LSQUnivariateSpline(
            x=timestamp, y=flux/max_flux, t=spline_t, k=3)
        ax_flux.plot(timestamp, smooth_flux(timestamp),
                     color=colors[idx])

        smooth_fwhm = scipy.interpolate.UnivariateSpline(
            x=timestamp, y=fwhm, k=3) #, w=None, bbox=[None, None], s=0.8*fwhm.shape[0])

        smooth_fwhm = scipy.interpolate.LSQUnivariateSpline(
            x=timestamp, y=fwhm, 
            t=spline_t,
            )
        #ax_fwhm.scatter(spline_t, smooth_fwhm(spline_t))
        #ax_fwhm.plot(spline_t, smooth_fwhm(spline_t), linewidth=5)
        ax_fwhm.plot(timestamp, smooth_fwhm(timestamp), 
                     #linewidth=4,
                     color=colors[idx])


        #print fwhm

        #
        # Make histograms
        #

        # h_flux_count, h_flux = numpy.histogram(flux/max_flux, bins=15, range=[min_flux/max_flux, 1.0])
        # h_flux_pos = 0.5*(h_flux[:-1] + h_flux[1:])
        # print h_flux
        # hist_flux.plot(h_flux_count, h_flux_pos,
        #              color=colors[idx])

        # h_fwhm_count, h_fwhm = numpy.histogram(fwhm, bins=15, range=[min_fwhm, max_fwhm])
        # h_fwhm_pos = 0.5*(h_fwhm[:-1] + h_fwhm[1:])
        # print h_fwhm
        #hist_fwhm.plot(h_fwhm_count, h_fwhm_pos,
        #             color=colors[idx])

        flux_range = (max_flux - min_flux)/max_flux
        x_flux = numpy.linspace(min_flux/max_flux-0.05*flux_range, 1.+0.05*flux_range, 100)
        density_flux = gaussian_kde(flux/max_flux)
        #density_flux.covariance_factor = lambda : .1
        density_flux._compute_covariance()
        peak_flux = 1. #numpy.max(density_flux(x_flux))
        hist_flux.plot(density_flux(x_flux)/peak_flux, x_flux, "-", color=colors[idx])

        fwhm_range = max_fwhm - min_fwhm
        x_fwhm = numpy.linspace(min_fwhm-0.05*fwhm_range, max_fwhm+0.05*fwhm_range, 100)
        density_fwhm = gaussian_kde(fwhm)
        #density_fwhm.covariance_factor = lambda : .1
        density_fwhm._compute_covariance()
        peak_fwhm = 1. #numpy.max(density_fwhm(x_fwhm))
        hist_fwhm.plot(density_fwhm(x_fwhm)/peak_fwhm, x_fwhm, "-", color=colors[idx])


        #
        #
        #
        all_flux.append(flux/max_flux)
        all_fwhm.append(fwhm)

        s2p = lambda x : 0.5*(1+scipy.special.erf(x/numpy.sqrt(2)))*100
        #print s2p([-3,-1,1,3])
        try:
            sigmas = scipy.stats.scoreatpercentile(flux/max_flux, s2p([-3,-1,1,3]))
            sigma1 = sigmas[2] - sigmas[1]
            sigma3 = sigmas[3] - sigmas[0]
        except:
            sigma1, sigma3 = -1., -1.
            pass
        txt = u'GS %d: 1\u03C3=%.3f - 3\u03C3=%.3f' % (
            idx+1, sigma1, sigma3)
        fig.text(text_positions[idx,0], text_positions[idx,1], 
                 txt,
                 horizontalalignment=text_align[idx])

        #hist_fwhm.hist(flux/max_flux, bins=15, orientation='horizontal')
        #axHisty.hist(y, bins=bins, orientation='horizontal')

        one_return = {
            'fwhm_min': min_fwhm,
            'fwhm_max': max_fwhm,
            'flux_min': min_flux,
            'flux_max': max_flux,
            'flux_1sigma': sigma1,
            'flux_3sigma': sigma3,
            'n_guide_samples': fwhm.shape[0]
        }
        guidestats[infile] = one_return

    #
    # Determine and set all axis limits
    #
    all_flux = numpy.array(all_flux)
    all_fwhm = numpy.array(all_fwhm)
    #print all_fwhm.shape
    #numpy.savetxt("fwhm", all_fwhm.reshape((-1,1)))

    try:
        _min_flux, _max_flux = numpy.min(all_flux), 1.0 # normalized !
    except:
        _min_flux, _max_flux = 0.0, 1.0
    ax_flux.set_ylim((_min_flux-0.05, 1.05))
    hist_flux.set_ylim((_min_flux-0.05, 1.05))

    try:
        _min_fwhm, _max_fwhm = numpy.min(all_fwhm), numpy.max(all_fwhm)
    except:
        _min_fwhm, _max_fwhm = 0.0, limit_fwhm_max
    #print _min_fwhm, _max_fwhm
    if (_max_fwhm > limit_fwhm_max): _max_fwhm = limit_fwhm_max
    _range_fwhm = _max_fwhm - _min_fwhm
    ax_fwhm.set_ylim((_min_fwhm-0.05*_range_fwhm, _max_fwhm+0.05*_range_fwhm))
    hist_fwhm.set_ylim((_min_fwhm-0.05*_range_fwhm, _max_fwhm+0.05*_range_fwhm))

    d_time = numpy.max([2, 0.03*max_time])

    ax_flux.set_xlim((0-d_time, max_time+d_time))
    ax_fwhm.set_xlim((0-d_time, max_time+d_time))

    hist_flux.set_ylim((global_min_flux-0.05, 1.05))

    if (not plot_filename == None):
        fig.set_size_inches(8,6)
        fig.savefig(plot_filename, dpi=100,
                    facecolor=fig.get_facecolor(), edgecolor='none')
    else:
        fig.show()
        plot.show()
        pass


    # catfile = infile[:-5]+".cat"

    # sex_cmd = """

    # %(sex)s -c %(qr_base)s/config/wcsfix.sex
    # -PARAMETERS_NAME %(qr_base)s/config/wcsfix.sexparam
    # -CATALOG_NAME %(catfile)s
    # %(infile)s """ % {
    #     'sex': sitesetup.sextractor,
    #     'qr_base': qr_dir,
    #     'catfile': catfile,
    #     'infile': infile,
    # }
    # print " ".join(sex_cmd.split())
    # os.system(" ".join(sex_cmd.split()))

    return guidestats

    
    
if __name__ == "__main__":

    filelist = sys.argv[1:]
    
    guidestats = draw_guidestarplot(filelist, title="some title", plot_filename="guide.png")
    print guidestats
