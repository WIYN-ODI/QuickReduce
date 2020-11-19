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
import astropy.io.fits as pyfits
import numpy
import math

import matplotlib
import matplotlib.pyplot
from podi_observingplots import *
from podi_definitions import *
from podi_commandline import *

dzp_limit_max = 0.3
dzp_limit_min = -3.


def create_zptrend_plot(files, output_filename, show_plot=False, zpoffset=0.):

    print "Reading data"
    direntry, arrays = read_data_from_files(files)
    obstype, exptime, filtername, photzp, photzpe, mjd, dateobs, airmass = arrays

    print "Plotting"
    fig, ax = matplotlib.pyplot.subplots()
    # tfig, tax = matplotlib.pyplot.subplots()

    cc = matplotlib.colors.ColorConverter()

    # This is the MJD of 01/01/0001
    mjd_zeropoint = 1721424.500000 - 2400000.5 + (7. / 24.0)

    def dzp_to_transparency(d_zp):
        return 100. * numpy.power(10., 0.4 * d_zp)

    all_dzps = numpy.array([0.2, -0.3])
    all_dzp_errs = numpy.array([0., 0., ])

    # print all_dzps
    for thisfilter in set(filtername):

        print thisfilter

        legendname = None
        color = 'grey'
        if (thisfilter in known_filters):
            legendname = thisfilter
            ref_zp, airmassterm, color = known_filters[thisfilter]

        # Select all datapoints for this filter
        data = []
        colors = []
        for i in range(len(filtername)):
            if (not filtername[i] == thisfilter):
                continue

            if (thisfilter in known_filters):
                ref_zp, airmassterm, color = known_filters[thisfilter]
                d_zp = photzp[i] - ref_zp + airmassterm * (airmass[i] - 1)
                zperr = photzpe[i]
                this_color = color
            else:
                d_zp = -999

            d_zpx = d_zp

            if (d_zp > 0.5 or d_zp < -5):
                d_zp = 0
                zperr = 0
                this_color = "grey"

            this_data = [mjd[i], exptime[i], photzp[i], zperr, d_zp, d_zpx]
            data.append(this_data)
            colors.append(this_color)  # cc.to_rgba(this_color))
            # colors.append(this_color) #cc.to_rgba(this_color))

        data = numpy.array(data)
        timestamp = data[:, 0] + 0.5 * data[:, 1] / 86400. - mjd_zeropoint  # mjdobs + 0.5*exptime
        all_dzps = numpy.append(all_dzps, data[:, 5])
        all_dzp_errs = numpy.append(all_dzp_errs, data[:, 3])
        # print all_dzps

        # print colors
        # print list(colors)
        # print "\n\n\n\n"
        ax.errorbar(x=timestamp, y=data[:, 4]-zpoffset, xerr=0.5 * data[:, 1] / 86400., yerr=data[:, 3],
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
    # locator = matplotlib.dates.AutoDateLocator()
    # hfmt = matplotlib.dates.AutoDateFormatter(locator)

    ax.set_ylabel("ZP difference Ref-ODI / Throughput loss /\nCloud attenuation [mag]")
    # ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator())
    ax.xaxis.set_major_formatter(hfmt)
    ax.set_ylim(bottom=0)
    matplotlib.pyplot.xticks(rotation='vertical')
    # ax.get_xaxis.set_xticks(rotation='vertical')
    matplotlib.pyplot.subplots_adjust(bottom=.22, top=0.97, right=0.96)

    time_start = numpy.min(mjd) - mjd_zeropoint
    time_end = numpy.max(numpy.array(mjd) + numpy.array(exptime) / 86400) - mjd_zeropoint

    if ((time_end - time_start) > 8):
        round_off = 24.  # 1 hour
    else:
        round_off = 48.

    # Round times to closest hours
    time_start = math.floor(time_start * round_off) / round_off
    time_end = math.ceil(time_end * round_off) / round_off
    matplotlib.pyplot.hlines(0, time_start, time_end)

    if ((time_end - time_start) * 24. > 10):
        # if time-range exceeds 10 hours, plot labels every 2 hours
        major = matplotlib.dates.HourLocator(interval=2)
    else:
        # otherwise every hour
        major = matplotlib.dates.HourLocator(interval=1)
    # Add minor ticks every 15 minutes
    minor = matplotlib.dates.MinuteLocator(interval=15)

    ax.xaxis.set_major_locator(major)
    ax.xaxis.set_minor_locator(minor)

    ax.set_xlim((time_start, time_end))

    # Determine the min and max y-range
    good_d_zp = (all_dzps < 50) & (all_dzps > -50)
    # print all_dzps, numpy.array(all_dzps)

    min_dzp = numpy.min((all_dzps - all_dzp_errs)[good_d_zp])
    max_dzp = numpy.max((all_dzps + all_dzp_errs)[good_d_zp])
    min_dzp = min_dzp if min_dzp > dzp_limit_min else dzp_limit_min
    max_dzp = max_dzp if max_dzp < dzp_limit_max else dzp_limit_max

    ax.set_ylim((min_dzp - 0.1, max_dzp + 0.1))
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
                     [start, top_level - height],
                     [end, top_level - height],
                     [complete, top_level]
                     ]

        this_color = 'grey'
        if (this_file['FILTER'] in known_filters):
            zp, amt, col = known_filters[this_file['FILTER']]
            this_color = col

        efficiency_plot.append(this_poly)
        efficiency_colors.append(this_color)

    # and then plots all the polygons
    coll = matplotlib.collections.PolyCollection(efficiency_plot,
                                                 facecolor=efficiency_colors,
                                                 # edgecolor='#808080',
                                                 edgecolor=efficiency_colors,
                                                 # edgecolor='none',
                                                 linestyle='-')
    ax.add_collection(coll)

    # Set output size to 900x500 pixels
    fig.set_size_inches(9, 5)
    print "Saving output to file", output_filename
    fig.savefig(output_filename, dpi=100)

    if (show_plot):
        matplotlib.pyplot.show()


if __name__ == "__main__":

    output_filename = cmdline_arg_set_or_default("-output", "photzp_trend.png")
    show = cmdline_arg_isset("-show")
    zpoffset = float(cmdline_arg_set_or_default("-zpoffset", 0.0))

    watch = cmdline_arg_isset("-watch")
    old_filelist = []
    if (watch):

        watch_interval = float(cmdline_arg_set_or_default("-interval", 1.0))
        day_offset = int(cmdline_arg_set_or_default("-offset", 0.))

        # in minutes

        # select all of tonights files, starting from noon
        import datetime, glob, time
        now = datetime.datetime.today() + datetime.timedelta(days=day_offset)
        print now

        if (now.hour < 12):
            # it's morning
            start_day = now - datetime.timedelta(days=1)
            end_day = now
        else:
            # it's afternoon
            start_day = now
            end_day = now + datetime.timedelta(days=1)

        print start_day, end_day

        afternoon = start_day.strftime("%Y%m%d")+"T1[2-9]*.fits"
        night_firsthalf = start_day.strftime("%Y%m%d")+"T2[0-4]*.fits"
        night_secondhalf = end_day.strftime("%Y%m%d")+"T0[0-9]*.fits"
        wildcards = "%s %s %s" % (afternoon, night_firsthalf, night_secondhalf)
        print wildcards

        while (True):
            files = glob.glob(afternoon) + glob.glob(night_firsthalf) + glob.glob(night_secondhalf)

            if (len(files) > len(old_filelist)):
                print files
                create_zptrend_plot(files, output_filename, zpoffset=zpoffset)

            old_filelist = files
            time.sleep(watch_interval*60)

    else:
        files = get_clean_cmdline()[1:]
        create_zptrend_plot(files, output_filename, zpoffset=zpoffset)