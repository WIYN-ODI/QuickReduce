#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import astropy.io.fits as pyfits
import numpy
import math
import ephem
import datetime
import glob
import time

from podi_commandline import *

import podi_almanac

ephem2mjd = +2400000.5 - 2415020
# print(ephem2mjd)

if __name__ == "__main__":

    # select all of tonights files, starting from noon
    now = datetime.datetime.today() #+ datetime.timedelta(days=day_offset)
    print(now)

    if (now.hour < 12):
        # it's morning
        start_day = now - datetime.timedelta(days=1)
        end_day = now
    else:
        # it's afternoon
        start_day = now
        end_day = now + datetime.timedelta(days=1)

    # print(start_day, end_day)

    wiyn = ephem.Observer()
    wiyn.lat = podi_almanac.wiyn_lat
    wiyn.lon = podi_almanac.wiyn_lon
    wiyn.elevation = podi_almanac.wiyn_elevation
    wiyn.pressure = 0
    wiyn.horizon = '-18'
    wiyn.date = start_day.strftime("%Y/%m/%d 19:00")

    sunset = wiyn.next_setting(ephem.Sun(), use_center=True)
    sunrise = wiyn.next_rising(ephem.Sun(), use_center=True)
    print("SUNSET: ", sunset)
    print("SUNRISE:", sunrise)

    night_start = sunset - 15 * ephem.minute - ephem2mjd
    night_end = sunrise + 15 * ephem.minute - ephem2mjd
    print(night_start, night_end)
    night_length = night_end - night_start

    afternoon = start_day.strftime("%Y%m%d")+"T1[2-9]*.fits"
    night_firsthalf = start_day.strftime("%Y%m%d")+"T2[0-4]*.fits"
    night_secondhalf = end_day.strftime("%Y%m%d")+"T0[0-9]*.fits"
    wildcards = "%s %s %s" % (afternoon, night_firsthalf, night_secondhalf)
    print(wildcards)

    files = glob.glob(afternoon) + glob.glob(night_firsthalf) + glob.glob(night_secondhalf)

    # Read all files, extract seeing measurements
    seeing_data = []
    for fn in files:
        print("Checking %s" % (fn))
        hdu = pyfits.open(fn)
        seeing_data.append([hdu[0].header['MJD-OBS'],
                            hdu[0].header['SEEING']])
        hdu.close()
    seeing_data = numpy.array(seeing_data)
    numpy.savetxt("seeing.data", seeing_data)
    seeing_data = seeing_data[numpy.isfinite(seeing_data[:,1])]

    for quarter in range(4):
        mjd_start = night_start + quarter*night_length/4.
        mjd_end = mjd_start + night_length/4.

        in_this_quarter = (seeing_data[:,0] >= mjd_start) & (seeing_data[:,0] <= mjd_end) & \
                          (seeing_data[:,1] > 0) & (seeing_data[:,1] < 5.0)

        seeing = numpy.nanmedian(seeing_data[in_this_quarter,1])
        seeing_min = numpy.min(seeing_data[in_this_quarter,1])
        seeing_max = numpy.max(seeing_data[in_this_quarter,1])
        print("Quarter %d: seeing = %.2f [%.2f..%.2f] (mjd=%.2f--%2f, %d **)" % (
            quarter+1, seeing, seeing_min, seeing_max,
            mjd_start, mjd_end, numpy.sum(in_this_quarter)))
