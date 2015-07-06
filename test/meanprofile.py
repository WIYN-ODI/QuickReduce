#!/usr/bin/env python


import pyfits
import numpy
import math
import scipy
import os, sys

_dir, _ = os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append(_dir+"/../")

import matplotlib.pyplot

from podi_commandline import *
from podi_definitions import clobberfile
from profile import *


def compute_mean_profile(filename, 
                         fxy_list, 
                         mode='radial',
                         mot_xy=(0.,0.),
                         dz2_limit=None,
                         max_radius=5.0,
                         save_tmp=True
                         ):

    
    mot_x, mot_y = mot_xy

    hdu = pyfits.open(filename)
    arcsec_per_pixel = hdu[0].header['CD2_2'] * 3600
    data = hdu[0].data.T

    # Convert zeopoint from mag/pixel to mag/arcsec
    pixelarea = arcsec_per_pixel**2
    width = int(math.ceil(max_radius / arcsec_per_pixel))

    if (not dz2_limit == None):
        dz2_limit /= arcsec_per_pixel

    if (mode == 'radial'):

        all_r = numpy.array([])
        all_data = numpy.array([])

        for fxy in fxy_list:
            if (type(fxy) == numpy.ndarray):
                fx = fxy[0] - 1.
                fy = fxy[1] - 1.
            else:
                items = fxy.split(",")
                print fxy, items
                fx = float(items[0]) - 1.
                fy = float(items[1]) - 1.

            r, cutout, nf = get_profile(data, center_x=fx, center_y=fy, 
                                    mx=mot_x,  my=mot_y, width=width, 
                                    mode='radial', 
                                    normalize=True,
                                )
            all_r = numpy.append(all_r,r)
            all_data = numpy.append(all_data,cutout)

            if (save_tmp):
                numpy.savetxt("radial__%.2f__%.2f.dump" % (fx,fy), 
                              numpy.append(r.reshape(-1,1), cutout.reshape(-1,1), axis=1))

        # Sort by radius
        si = numpy.argsort(all_r)
        all_r = all_r[si]
        all_data = all_data[si]

        # combine radii and data and save to file
        out = numpy.zeros((all_r.shape[0], 2))
        out[:,0] = all_r * arcsec_per_pixel
        out[:,1] = all_data
        return out

    elif (mode == 'profile'):

        all_dz = numpy.array([])
        all_dz2 = numpy.array([])
        all_data = numpy.array([])

        for idx, fxy in enumerate(fxy_list):
            items = fxy.split(",")
            print fxy, items
            fx = float(items[0]) - 1.
            fy = float(items[1]) - 1.

            distance_z, distance_z2, cutout, nf = get_profile(
                data, center_x=fx, center_y=fy, 
                mx=mot_x,  my=mot_y, width=width, 
                mode='crosscut', 
                normalize=True,
                dz2_limit=dz2_limit,
            )
            all_dz = numpy.append(all_dz, distance_z)
            all_dz2 = numpy.append(all_dz2, distance_z2)
            all_data = numpy.append(all_data, cutout)

            numpy.savetxt("profile__%.2f__%.2f.dump" % (fx,fy), 
                          numpy.append(distance_z.reshape(-1,1), cutout.reshape(-1,1), axis=1))
            numpy.savetxt("profile__%d__%.2f__%.2f.dump_2" % (idx+1,fx,fy), 
                          numpy.append(all_dz.reshape(-1,1), all_data.reshape(-1,1), axis=1))

        si = numpy.argsort(all_dz)
        all_dz = all_dz[si]
        all_dz2 = all_dz2[si]
        all_data = all_data[si]

        out = numpy.zeros((all_dz.shape[0], 3))
        out[:,0] = all_dz  * arcsec_per_pixel
        out[:,1] = all_data
        out[:,2] = all_dz2  * arcsec_per_pixel

        return out

    else:
        print "mode unknown:", mode
        

    return

if __name__ == "__main__":

    radial_mode = cmdline_arg_isset("-radial")
    profile_mode = cmdline_arg_isset("-profile")
    if (not radial_mode and not profile_mode):
        print "Need to specify either -radial or -profile"
        sys.exit(0)

    savefile = cmdline_arg_set_or_default("-save", None)
    if (not savefile == None):
        clobberfile(savefile)

    infilename = get_clean_cmdline()[1]

    mot_x, mot_y = 0, 0
    mxy = cmdline_arg_set_or_default("-nonsidereal", None)
    if (not mxy == None):
        items = mxy.split(",")
        mot_x = float(items[0])
        mot_y = float(items[1])


    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    dz2_limit = float(cmdline_arg_set_or_default("-dz2", 0.5))
    max_radius = float(cmdline_arg_set_or_default("-maxr", 5.0))

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
