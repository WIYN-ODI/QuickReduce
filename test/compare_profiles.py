#!/usr/bin/env python


import pyfits
import numpy
import math
import scipy
import os, sys
import scipy.interpolate

_dir, _ = os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append(_dir+"/../")

import matplotlib.pyplot

from podi_commandline import *
from podi_definitions import clobberfile
from profile import *

from meanprofile import compute_mean_profile


if __name__ == "__main__":

    mot_x, mot_y = 0, 0
    mxy = cmdline_arg_set_or_default("-nonsidereal", None)
    if (not mxy == None):
        items = mxy.split(",")
        mot_x = float(items[0])
        mot_y = float(items[1])

    aperture_size = float(cmdline_arg_set_or_default("-apsize", 3.))

    dz2_limit = float(cmdline_arg_set_or_default("-dz2", 0.5))
    max_radius = float(cmdline_arg_set_or_default("-maxr", 5.0))

    # Compute mean profile for all stars
    ref_frame = get_clean_cmdline()[1]
    coords = []
    for idx, p in enumerate(get_clean_cmdline()[2:]):
        if (p == ":"):
            break
        coords.append(p)
        print p, idx

    #
    # Compute mean profile
    #
    ref_data = compute_mean_profile(ref_frame,
                                    coords, 
                                    mode='profile',
                                    mot_xy=(mot_x, mot_y),
                                    dz2_limit=dz2_limit,
                                )
    
    print idx

    print "\n\n\n\n\n"
    print get_clean_cmdline()[idx+3]
    print "\n\n\n\n\n"

    frame2 = get_clean_cmdline()[idx+3]
    print frame2
    coords2 = []
    for idx2, p in enumerate(get_clean_cmdline()[idx+4:]):
        print p
        if (p == ":"):
            break
        coords2.append(p)

    data = compute_mean_profile(frame2,
                                coords2, 
                                mode='radial',
                                mot_xy=(mot_x, mot_y),
                                dz2_limit=dz2_limit,
                                )

    numpy.savetxt("data.1", ref_data)
    numpy.savetxt("data.2", data)

    print coords
    print coords2

    #
    # Now we have both PSF profiles.
    # fit a spline curve to the reference PSF
    #
    ref_data_mirrored = numpy.append((ref_data * numpy.array([-1.,1.,0.0]))[::-1],
                                     ref_data,
                                     axis=0)
    numpy.savetxt("data.1m", ref_data_mirrored)

    #psf_fit = scipy.interpolate.interp1d(ref_data[:,0], ref_data[:,1])
    # psf_fit = scipy.interpolate.interp1d(ref_data_mirrored[:,0], ref_data_mirrored[:,1])
    # psf_fit = scipy.interpolate.InterpolatedUnivariateSpline(ref_data[:,0], ref_data[:,1])
    base_points = numpy.linspace(-4.8,4.8,200)
    #base_points = numpy.linspace(0.02,4.8,200)
    print numpy.min(ref_data_mirrored[:,0]), numpy.max(ref_data_mirrored[:,0])
    print numpy.min(base_points), numpy.max(base_points)
    psf_fit = scipy.interpolate.LSQUnivariateSpline(
        ref_data_mirrored[:,0], ref_data_mirrored[:,1],
        t=base_points,
        k=3)
    # psf_fit = scipy.interpolate.LSQUnivariateSpline(
    #     ref_data[:,0], ref_data[:,1],
    #     t=base_points,
    #     k=3)
    
    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    ax.scatter(ref_data[:,0], ref_data[:,1])
    
    plot_x = numpy.linspace(numpy.min(ref_data_mirrored[:,0]), numpy.max(ref_data_mirrored[:,0]), 1000)
    ax.plot(plot_x, psf_fit(plot_x), "r-", linewidth=3)
    ax.scatter(base_points, psf_fit(base_points), marker="o", color="green", s=50)

    fig.show()

    #
    # Now in a second step, compute synthetic PSF at 
    # all points in the analysis PSF
    #
    data_mirrored = numpy.append((data*numpy.array([-1.,1.]))[::-1],
                                 data, axis=0)
    data_psf = scipy.interpolate.LSQUnivariateSpline(
        data_mirrored[:,0], data_mirrored[:,1],
        t=base_points,
        k=3)
    
    #
    # Normalize both PSFs to the same integrated intensity
    #
    aperture_x = numpy.linspace(-aperture_size, aperture_size, 500)
    ref_flux = numpy.sum(psf_fit(aperture_x))
    data_flux = numpy.sum(data_psf(aperture_x))

    model_psf = psf_fit(data[:,0])
    scale_match = data_flux/ref_flux
    difference = data[:,1] - model_psf*scale_match
    print "scaling reference psf by %7.5f" % (scale_match)

    #
    # Make a second plot showing the data, reference psf, and difference
    #
    fig2 = matplotlib.pyplot.figure()
    ax2 = fig2.add_subplot(111)
    
    numpy.savetxt("data.diff", 
                  numpy.append(data, difference.reshape((-1,1)), axis=1))

    ax2.scatter(data[:,0], data[:,1], marker=".")
    ax2.scatter(data[:,0], difference-0.2, marker="x")
    ax2.plot(plot_x, psf_fit(plot_x)*scale_match, "r-", linewidth=3)
    ax2.set_xlim((0, numpy.max(base_points)))
    ax2.axhline(color='#80ff80', linewidth=3)
    ax2.axhline(y=-0.2, color='#8080ff', linewidth=3)
    ax2.set_xlabel("Radius [arcsec]")
    ax2.set_ylabel("Intensity [flux-normalized counts]")
    fig2.show()

    matplotlib.pyplot.show()

    sys.exit(0)




    if (radial_mode):

        data = compute_mean_profile(infilename, 
                                    get_clean_cmdline()[2:], 
                                    mode='radial',
                                    mot_xy=(mot_x, mot_y),
                                    dz2_limit=None,
                                )

        if (not savefile == None):
            f = open(savefile, "a")
            numpy.savetxt(f, data)
            print >>f, "\n\n\n\n\n"
            f.close()

        ax.set_xlabel("Radius [arcsec]")
        ax.set_ylabel("Intensity [flux-normalized counts]")
        ax.scatter(data[:,0], data[:,1])
        ax.set_xlim((0, max_radius))

    if (profile_mode):

        data = compute_mean_profile(infilename, 
                                    get_clean_cmdline()[2:], 
                                    mode='profile',
                                    mot_xy=(mot_x, mot_y),
                                    dz2_limit=dz2_limit,
                                )

        # all_dz = numpy.array((0))
        # all_dz2 = numpy.array((0))
        # all_data = numpy.array((0))

        # for fxy in get_clean_cmdline()[2:]:
        #     items = fxy.split(",")
        #     print fxy, items
        #     fx = float(items[0]) - 1.
        #     fy = float(items[1]) - 1.

        #     distance_z, distance_z2, cutout = get_profile(
        #         data, center_x=fx, center_y=fy, 
        #         mx=mot_x,  my=mot_y, width=15, 
        #         mode='crosscut', 
        #         normalize=normalize,
        #         dz2_limit=dz2_limit,
        #     )
        #     all_dz = numpy.append(all_dz, distance_z)
        #     all_dz2 = numpy.append(all_dz2, distance_z2)
        #     all_data = numpy.append(all_data, cutout)

        # si = numpy.argsort(all_dz)
        # all_dz = all_dz[si]
        # all_dz2 = all_dz2[si]
        # all_data = all_data[si]

        if (not savefile == None):
            f = open(savefile, "a")
            # out = numpy.zeros((all_dz.shape[0], 3))
            # out[:,0] = all_dz  * arcsec_per_pixel
            # out[:,1] = all_data
            # out[:,2] = all_dz2  * arcsec_per_pixel
            numpy.savetxt(f, data)
            print >>f, "\n\n\n\n\n"
            f.close()


        # to_plot = (all_dz2<2)
        # ax.scatter((all_dz*arcsec_per_pixel)[to_plot], 
        #            all_data[to_plot])
        ax.scatter(data[:,0], data[:,1])
        ax.set_xlim((0, max_radius))
        print max_radius

        # to_plot = (all_dz < 2)
        # ax.scatter((all_dz2*arcsec_per_pixel)[to_plot], 
        #            all_data[to_plot], 
        #            c='red', marker='x')
        ax.scatter(data[:,2], data[:,1],
                   c='red', marker='x')
        
    fig.show()
    matplotlib.pyplot.show()
