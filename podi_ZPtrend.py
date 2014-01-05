#!/usr/bin/env python
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

Create a night-log for a given list of input frame. 

This log lists, for each frame 
* filename
* type of observation: bias/dark/flat/science 
* binning
* filter name
* exposure time
* Object/target name as specified during the observation (this is not necessarily
  the real name of the target, rather what the user entered)
* pointing coordinates: Ra/Dec

The output of this file is compatible with the file-list requirements of, for 
example, podi_makecalibrations.

"""

import sys
import os
import pyfits
import numpy
import math

import matplotlib
import matplotlib.pyplot
from podi_observingplots import *


if __name__ == "__main__":

    
    print "Reading data"
    direntry, arrays = read_data_from_files(sys.argv[1:])
    obstype, exptime, filtername, photzp, photzpe, mjd, dateobs = arrays

    print mjd

    print "Plotting"
    fig, ax = matplotlib.pyplot.subplots()
    #tfig, tax = matplotlib.pyplot.subplots()

    cc = matplotlib.colors.ColorConverter()

    # This is the MJD of 01/01/0001
    mjd_zeropoint = 1721424.500000 - 2400000.5 + (7./24.0)

    def dzp_to_transparency(d_zp):
        return 100.*numpy.power(10., 0.4*d_zp)

    for thisfilter in set(filtername):

        print thisfilter

        legendname = None
        color='grey'
        if (thisfilter in known_filters):
            legendname = thisfilter
            ref_zp, color = known_filters[thisfilter]

        # Select all datapoints for this filter
        data = []
        colors = []
        for i in range(len(filtername)):
            if (not filtername[i] == thisfilter):
                continue
                
            if (thisfilter in known_filters):
                ref_zp, color = known_filters[thisfilter]
                d_zp = photzp[i] - ref_zp
                zperr = photzpe[i]
                this_color = color
            else:
                d_zp = -999

            if (d_zp > 0.5 or d_zp < -5):
                d_zp = 0
                zperr = 0
                this_color = "grey"

            this_data = [mjd[i], exptime[i], photzp[i], zperr, d_zp]
            data.append(this_data)
            colors.append(this_color) #cc.to_rgba(this_color))
            #colors.append(this_color) #cc.to_rgba(this_color))

        data = numpy.array(data)
        timestamp = data[:,0] + 0.5 * data[:,1] / 86400.  - mjd_zeropoint# mjdobs + 0.5*exptime
        
        print colors
        print list(colors)
        print "\n\n\n\n"
        ax.errorbar(x=timestamp, y=data[:,4], xerr = 0.5*data[:,1]/86400., yerr=data[:,3],
#                    color=[colors[i] for i in range(len(colors))],
                    c=color,
                    marker="o",
                    fmt="o",
                    label=legendname)
        # ax.plot(x=timestamp, y=data[:,4],
        #         color=colors,
        #         marker="o",
        #         label=legendname)

        # tax.errorbar(x=timestamp, y=dzp_to_transparency(data[:,4]), xerr = 0.5*data[:,1]/86400., 
        #              yerr=data[:,3]*dzp_to_transparency(data[:,4]),
        #              marker="o",
        #              fmt="o",
        #              label=legendname)

    # matplotlib date format object
    hfmt = matplotlib.dates.DateFormatter('%m/%d/%Y\n%H:%M MST')
    #locator = matplotlib.dates.AutoDateLocator()
    #hfmt = matplotlib.dates.AutoDateFormatter(locator)

    ax.set_ylabel("ZP difference Ref-ODI / Throughput loss /\nCloud attenuation [mag]")
    #ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator())
    ax.xaxis.set_major_formatter(hfmt)
    ax.set_ylim(bottom = 0)
    matplotlib.pyplot.xticks(rotation='vertical')
    #ax.get_xaxis.set_xticks(rotation='vertical')
    matplotlib.pyplot.subplots_adjust(bottom=.22, top=0.95, right=0.9)

    time_start = numpy.min(mjd)-mjd_zeropoint
    time_end = numpy.max(numpy.array(mjd)+numpy.array(exptime)/86400)-mjd_zeropoint
    frac_diff = 0.05 * (time_end - time_start)
    time_start, time_end = time_start - frac_diff, time_end + frac_diff
    matplotlib.pyplot.hlines(0,time_start, time_end)


    ax.set_xlim((time_start, time_end))

    ax.set_ylim((-5,0.3))
    ax.legend(loc='best', borderaxespad=1)

    #
    # Now add the polygons to show the shutter-open efficiency.
    #
    # First, compute all times for each of the frames
    top_level = 0.3
    height = 0.2
    efficiency_plot = []
    efficiency_colors = []
    for filename in direntry:
        
        this_file = direntry[filename]

        # compute all times
        # but first apply the MJD zeropoint to convert times to the matplotlib
        # format
        mjdobs = this_file['MJD-OBS'] - mjd_zeropoint
        init = mjdobs - seconds2mjd(10.)
        start = mjdobs
        end = mjdobs + seconds2mjd(this_file['EXPMEAS'])
        complete = end + seconds2mjd(25.)
        this_poly = [[init, top_level],
                     [start, top_level-height],
                     [end, top_level-height],
                     [complete, top_level]
                     ]

        this_color = 'grey'
        if (this_file['FILTER'] in known_filters):
            zp,col = known_filters[this_file['FILTER']]
            this_color = col
            
        efficiency_plot.append(this_poly)
        efficiency_colors.append(this_color)

    # and then plots all the polygons
    coll = matplotlib.collections.PolyCollection(efficiency_plot,
                                                 facecolor=efficiency_colors,
                                                 #edgecolor='#808080', 
                                                 edgecolor=efficiency_colors,
                                                 #edgecolor='none',
                                                 linestyle='-')
    ax.add_collection(coll)

    # Set output size to 900x500 pixels
    fig.set_size_inches(9,5)
    fig.savefig("photzp_trend.png",dpi=100)

    matplotlib.pyplot.show()


