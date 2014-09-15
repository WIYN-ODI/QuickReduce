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
from podi_definitions import get_binning

import datetime
import ephem
import numpy


if __name__ == "__main__":

    for filename in sys.argv[1:]:

        if (not os.path.isfile(filename)):
            continue

        try:
            hdulist = pyfits.open(filename)
            hdr = hdulist[0].header
        except:
            continue

        timeobs = hdr['TIME-OBS']
        
        ra = hdr['RA']
        dec = hdr['DEC']

        date_format = "%Y-%m-%dT%H:%M:%S.%f"
        ephem_format = "%Y/%m/%d %H:%M:%S"

        time_format = "%H:%M:%S.%f"
        expmeas = hdulist[0].header['EXPMEAS'] if 'EXPMEAS' in hdulist[0].header \
            else hdulist[0].header['EXPTIME']
        date_obs = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], date_format)
        date_mid = date_obs + datetime.timedelta(seconds=(0.5*expmeas))
        date_end = date_obs + datetime.timedelta(seconds=expmeas)
        mjd = hdulist[0].header['MJD-OBS']

        target = ephem.Equatorial(hdulist[0].header['RA'], hdulist[0].header['DEC'])
        print target.ra, target.dec


        sun = ephem.Sun()
        moon = ephem.Moon()

        wiyn = ephem.Observer()
        wiyn.lat = '31.9575'   # N is positive
        wiyn.lon = '-111.6008' # E is positive
        wiyn.elevation = 2096.
        #print datetime.datetime.strftime(date_mid, ephem_format)
        wiyn.date = datetime.datetime.strftime(date_mid, ephem_format)

        wiyn_astrotwilight = ephem.Observer()
        wiyn_astrotwilight.lat = '31.9575'   # N is positive
        wiyn_astrotwilight.lon = '-111.6008' # E is positive
        wiyn_astrotwilight.elevation = 2096.
        wiyn_astrotwilight.horizon = '-18.'
        wiyn_astrotwilight.date = datetime.datetime.strftime(date_mid, ephem_format)

        print "current time:", wiyn.date

        print "next sunrise", wiyn.next_rising(sun)
        print "last sunset ", wiyn.previous_setting(sun)
        print "next transit", wiyn.next_transit(sun)

        sun.compute(wiyn)
        print "SUN: ",sun.alt, sun.az, sun.ra, sun.dec

        moon.compute(wiyn)
        print "MOON:", moon.alt, moon.az, moon.ra, moon.dec
        print "MOON: phase=",moon.moon_phase

        print "moon-sun:", ephem.separation(moon, sun)

        body = ephem.FixedBody(target) #ra=target.ra, dec=target.dec, epoch=target.epoch)
        # body2 = ephem.FixedBody(ra=target.ra, dec=0, epoch=target.epoch)
        # body2.compute(wiyn)

        body._ra = target.ra
        body._dec = target.dec
        body._epoch = target.epoch
        body.compute(wiyn)
        print "body:", body.ra, body.dec
        print "body:", body.alt, body.az, hdulist[0].header['ZD']

        body2 = ephem.FixedBody()#_ra=target.ra, _dec=0, _epoch=target.epoch)
        body2._ra = target.ra
        body2._dec = 0
        body2._epoch = target.epoch
        body2.compute(wiyn)

        print "object-moon", ephem.separation(body, moon)
#        print "object-moon", ephem.separation((target.lon, target.lat), moon)
        print "object-sun", ephem.separation(body, sun)
        print "sun-sun", ephem.separation(sun, sun)
        print "body-eq", ephem.separation(body, body2), hdulist[0].header['DEC']

        hdulist[0].header['SUN__RA'] = (numpy.degrees(sun.ra),
                                        "R.A. of Sun during obs [deg]")
        hdulist[0].header['SUN__DEC'] = (numpy.degrees(sun.dec),
                                         "declination of Sun during obs [deg]")
        hdulist[0].header['SUN__ALT'] = (numpy.degrees(sun.alt),
                                         "altitude of Sun during obs [deg]")
        hdulist[0].header['SUN__AZ'] = (numpy.degrees(sun.az),
                                        "azimuth of Sun during obs [deg]")

        hdulist[0].header['MOON_RA'] = (numpy.degrees(moon.ra),
                                        "R.A. of Moon during obs [deg]")
        hdulist[0].header['MOON_DEC'] = (numpy.degrees(moon.dec),
                                         "declination of Moon during obs [deg]")
        hdulist[0].header['MOON_ALT'] = (numpy.degrees(moon.alt),
                                         "altitude of Moon during obs [deg]")
        hdulist[0].header['MOON_AZ'] = (numpy.degrees(moon.az),
                                        "azimuth of Moon during obs [deg]")

        hdulist[0].header['SUN__D'] = (numpy.degrees(ephem.separation(body, sun)),
                                       "angle between target and sun [deg]")
        hdulist[0].header['MOON_D'] = (numpy.degrees(ephem.separation(body, moon)),
                                       "angle between target and moon [deg]")

        hdulist[0].header['MOONPHSE'] = (moon.moon_phase,
                                         "Moon phase")
        
        if (moon.alt > 0):
            hdulist[0].header['MOON_RSE'] = (24.*(wiyn.date-wiyn.previous_rising(moon)), 
                                             'time since moon rise [hours]')
            hdulist[0].header['MOON_SET'] = (24.*(wiyn.next_setting(moon)-wiyn.date), 
                                             'time until moon set [hours]')
        else:
            hdulist[0].header['MOON_RSE'] = (24.*(wiyn.next_rising(moon)-wiyn.date), 
                                             'time until moon rise [hours]')
            hdulist[0].header['MOON_SET'] = (24.*(wiyn.date-wiyn.previous_setting(moon)), 
                                             'time since moon set [hours]')
            
        hdulist[0].header['SUN_RISE'] = (24.*(wiyn.next_rising(sun)-wiyn.date), 
                                         'time until sunrise [hours]')
        hdulist[0].header['SUN_SET'] = (24.*(wiyn.date-wiyn.previous_setting(sun)), 
                                        'time since sunset [hours]')

        hdulist[0].header['SUN_RS18'] = (24.*(wiyn_astrotwilight.next_rising(sun)-wiyn.date), 
                                         'time until astronomical twilight [hours]')
        hdulist[0].header['SUN_ST18'] = (24.*(wiyn.date-wiyn_astrotwilight.previous_setting(sun)), 
                                         'time since astronomical twilight [hours]')
        print hdulist[0].header

        print wiyn.date+0.
