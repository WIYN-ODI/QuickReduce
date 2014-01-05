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
import datetime

import matplotlib
import matplotlib.pyplot

from podi_observingplots import *
from podi_definitions import *



if __name__ == "__main__":

    print "Reading data"
    direntry, arrays = read_data_from_files(get_clean_cmdline()[1:])
    obstype, exptime, filtername, photzp, photzpe, mjd, dateobs, airmass = arrays

    # Now create the plots


    # This is the MJD of 01/01/0001
    mjd_zeropoint = 1721424.500000 - 2400000.5 + (7./24.0)
    # Find out the start and end times of the data block
    time_start = numpy.min(mjd)-mjd_zeropoint
    time_end = numpy.max(numpy.array(mjd)+numpy.array(exptime)/86400)-mjd_zeropoint

    print time_start*24, time_end*24.
    hour_start = math.floor(time_start*24.)
    hour_end = math.ceil(time_end*24.)
    n_hours = int(hour_end - hour_start)
    print n_hours, hour_start, hour_end

    
    print "setting up plot"
    fig = matplotlib.pyplot.figure()

    all_axes = []
    for cur_hour in range(n_hours):

        subplot_id = n_hours * 100 + 10 + (cur_hour+1)
        if (len(all_axes) == 0):
            ax = fig.add_subplot(n_hours, 1, cur_hour+1)
        else:
            ax = fig.add_subplot(n_hours, 1, cur_hour+1, sharey=all_axes[0])

        fiveminutes = matplotlib.dates.MinuteLocator(interval=5)
        minutes = matplotlib.dates.MinuteLocator(interval=1)
        hfmt = matplotlib.dates.DateFormatter(':%M')
        ax.xaxis.set_major_locator(fiveminutes)
        ax.xaxis.set_minor_locator(minutes)
        ax.xaxis.set_major_formatter(hfmt)
        ax.set_xlim(((hour_start+cur_hour)/24., (hour_start+cur_hour+1)/24.-1e-6))
        ax.axes.yaxis.set_ticks([])

        hour_of_the_day = int(hour_start + cur_hour) % 24
        #print hour_of_the_day

        date = datetime.date.fromordinal(int(math.floor((hour_start + cur_hour)/24.))).strftime("%d/%m/%y")
        
        if (n_hours > 7):
            ylabel = "%s-%02dh" % (date, hour_of_the_day)
        else:
            ylabel = "%s\n%02dh MST" % (date, hour_of_the_day)

        ax.set_ylabel(ylabel,
                      rotation="horizontal",
                      verticalalignment="center",
                      horizontalalignment="right")
        all_axes.append(ax)

    all_axes[0].set_ylim((0,1))
#    all_axes[0].set_title("Observing efficiency")
    fig.suptitle("Observing efficiency", fontsize=20)

    for i in range(len(all_axes)-1): #ax in all_axes[:-1]:
        #ax.axes.get_xaxis().set_visible(False)
        all_axes[i].axes.xaxis.set_ticklabels([])

    cc = matplotlib.colors.ColorConverter()


    def dzp_to_transparency(d_zp):
        return 100.*numpy.power(10., 0.4*d_zp)

    matplotlib.pyplot.subplots_adjust(bottom=.05, 
                                      top=0.92, 
                                      right=0.98, 
                                      left=0.13 if n_hours <= 7 else 0.18)

    #
    # Now add the polygons to show the shutter-open efficiency.
    #
    # First, compute all times for each of the frames
    top_level = 1
    height = 1
    efficiency_plot = []
    efficiency_colors = []

    poly_for_axes = [[]] * n_hours
    polyc_for_axes = [[]] * n_hours
    print poly_for_axes

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
        
        # Determine which hour we need
        hour_slot = int( math.floor(init*24.) - hour_start )
        print hour_start, init*24, hour_slot
        hour_slot=0

        this_color = 'grey'
        if (this_file['FILTER'] in known_filters):
            zp,amt,col = known_filters[this_file['FILTER']]
            this_color = col
            
        if (math.floor(init*24.) < math.floor(complete*24.)):
            # This block spans multiple hours
            hour_break = math.floor(complete*24.)/24.
            poly_start = [[init, top_level],
                         [start if start < hour_break else hour_break , top_level-height],
                         [end if end < hour_break else hour_break, top_level-height],
                         [complete if end < hour_break else hour_break, top_level]
            ]
            poly_end = [[init if init > hour_break else hour_break, top_level],
                        [start if start > hour_break else hour_break , top_level-height],
                        [end if end > hour_break else hour_break, top_level-height],
                        [complete if end > hour_break else hour_break, top_level]
            ]

            poly_for_axes[hour_slot].append(poly_start)
            polyc_for_axes[hour_slot].append(this_color)

            poly_for_axes[hour_slot+1].append(poly_end)
            polyc_for_axes[hour_slot+1].append(this_color)

            pass
        else:
            this_poly = [[init, top_level],
                         [start, top_level-height],
                         [end, top_level-height],
                         [complete, top_level]
            ]
            poly_for_axes[hour_slot].append(this_poly)
            polyc_for_axes[hour_slot].append(this_color)


    for i in range(len(all_axes)):

        # and then plots all the polygons
        coll = matplotlib.collections.PolyCollection(poly_for_axes[i],
                                                     facecolor=polyc_for_axes[i],
                                                     #edgecolor='#808080', 
                                                     edgecolor=polyc_for_axes[i],
                                                     #edgecolor='none',
                                                     linestyle='-')
        all_axes[i].add_collection(coll)


    #tax.set_ylim((1,200))
    #tax.set_yscale('log')
    #tax.legend(loc='best', borderaxespad=1)

    # Set output size to 900x500 pixels
    fig.set_size_inches(9,5)
    output_filename = cmdline_arg_set_or_default("-output", "shutter_open.png")
    print "Saving output to file",output_filename
    fig.savefig(output_filename, dpi=100)

    #tfig.savefig("transparency_trend.png")

    if (cmdline_arg_isset("-show")):
        matplotlib.pyplot.show()


