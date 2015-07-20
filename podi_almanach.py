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
from podi_definitions import get_binning, add_fits_header_title
from podi_commandline import *
import logging
import podi_logging

import datetime
import ephem
import numpy

wiyn_lat = '43.05' #'31.9575'   # N is positive
wiyn_lon = '-88.0' #'-111.6008' # E is positive
wiyn_elevation = 2096. # m


def add_ephem_data_to_header(hdr, time_overwrite):

    logger = logging.getLogger("SunMoonData")

    #
    # Read in the target coordinates
    #
    ra = hdr['RA']
    dec = hdr['DEC']
    target = ephem.Equatorial(hdr['RA'], hdr['DEC'])


    #
    # Convert the timing information from the header 
    # into a format that pyephem understands
    #
    date_format = "%Y-%m-%dT%H:%M:%S.%f"
    ephem_format = "%Y/%m/%d %H:%M:%S"
    time_format = "%H:%M:%S.%f"
    expmeas = hdr['EXPMEAS'] if 'EXPMEAS' in hdr \
        else hdr['EXPTIME']
    date_obs = datetime.datetime.strptime(hdr['DATE-OBS'], date_format)
    date_mid = date_obs + datetime.timedelta(seconds=(0.5*expmeas))
    date_end = date_obs + datetime.timedelta(seconds=expmeas)
    mjd = hdr['MJD-OBS']

    time_difference = 5*ephem.hour

    #
    # Create objects for sun, moon, and WIYN
    #
    sun = ephem.Sun()
    moon = ephem.Moon()

    wiyn = ephem.Observer()
    wiyn.lat = wiyn_lat
    wiyn.lon = wiyn_lon
    wiyn.elevation = wiyn_elevation
    wiyn.date = datetime.datetime.strftime(date_mid, ephem_format)

    if (not time_overwrite == None):
        wiyn.date = time_overwrite

    logger.debug("Current time (UTC): %s" % (wiyn.date))
    logger.debug("Current coordinates: %s %s" % (target.ra, target.dec))

    # #
    # # Compute some stats for astronomical twilight
    # #
    # wiyn_astrotwilight = ephem.Observer()
    # wiyn_astrotwilight.lat = wiyn_lat
    # wiyn_astrotwilight.lon = wiyn_lon
    # wiyn_astrotwilight.elevation = wiyn_elevation
    # wiyn_astrotwilight.horizon = '-18.'
    # # wiyn_astrotwilight.date = datetime.datetime.strftime(date_mid, ephem_format)
    # sun_astrotwi = ephem.Sun()
    # sun_astrotwi.compute(wiyn_astrotwilight)

    logger.debug("next sunrise: %s" % (ephem.localtime(wiyn.next_rising(ephem.Sun()))))
    logger.debug("next sunset : %s" % (ephem.localtime(wiyn.next_setting(ephem.Sun()))))
    logger.debug("last sunset : %s" % (ephem.localtime(wiyn.previous_setting(ephem.Sun()))))
    logger.debug("next transit: %s" % (ephem.localtime(wiyn.next_transit(ephem.Sun()))))

    # Compute the sun's position at WIYN
    sun.compute(wiyn)
    logger.debug("SUN:  alt=%-10s az=%-10s   RA=%-10s DEC=%-10s" % (
        sun.alt, sun.az, sun.ra, sun.dec)) #, " alt/az=", numpy.degrees(sun.alt)

    # and the Moon's position at WIYN
    moon.compute(wiyn)
    logger.debug("MOON: alt=%-10s az=%-10s   RA=%-10s DEC=%-10s  PHASE:%s" % (
        moon.alt, moon.az, moon.ra, moon.dec, moon.moon_phase)) 

    # Compute the position of the target on sky
    body = ephem.FixedBody(target)
    body._ra = target.ra
    body._dec = target.dec
    body._epoch = target.epoch
    body.compute(wiyn)

    # print "object-moon", ephem.separation(body, moon)
    # print "object-sun", ephem.separation(body, sun)
    # print "sun-sun", ephem.separation(sun, sun)
    # print "body-eq", ephem.separation(body, body2), hdr['DEC']

    logger.debug("Writing data to FITS header")

    #
    # Compute positions for sun and moon, and separation to the target 
    #

    # positions for the sun
    hdr['SUN__RA']  = (numpy.degrees(sun.ra),
                       "R.A. of Sun during obs [deg]")
    hdr['SUN__DEC'] = (numpy.degrees(sun.dec),
                       "declination of Sun during obs [deg]")
    hdr['SUN__ALT'] = (numpy.degrees(sun.alt),
                       "altitude of Sun during obs [deg]")
    hdr['SUN__AZ']  = (numpy.degrees(sun.az),
                       "azimuth of Sun during obs [deg]")

    # Positions for the moon
    hdr['MOON_RA']  = (numpy.degrees(moon.ra),
                       "R.A. of Moon during obs [deg]")
    hdr['MOON_DEC'] = (numpy.degrees(moon.dec),
                       "declination of Moon during obs [deg]")
    hdr['MOON_ALT'] = (numpy.degrees(moon.alt),
                       "altitude of Moon during obs [deg]")
    hdr['MOON_AZ']  = (numpy.degrees(moon.az),
                       "azimuth of Moon during obs [deg]")

    # separation between sun/moon and target
    hdr['SUNANGLE']   = (numpy.degrees(ephem.separation(body, sun)),
                       "angle between target and sun [deg]")
    hdr['MOONANGL']   = (numpy.degrees(ephem.separation(body, moon)),
                       "angle between target and moon [deg]")
    logger.debug("Distance object-moon (sun): %7.3f (%7.3f)" % (
        ephem.separation(body, moon), ephem.separation(body, sun)))


    # Moon phase: fraction of moon that's illuminated by the sun
    hdr['MOONPHSE'] = (moon.moon_phase,
                       "Moon phase")


    # print
    # print "last sun transit    ", wiyn.previous_transit(ephem.Sun())
    # print "last sun rise       ", wiyn.previous_rising(ephem.Sun())
    # print "last sun set        ", wiyn.previous_setting(ephem.Sun())
    # print "next sun transit    ", wiyn.next_transit(ephem.Sun())
    # print "next sun rise       ", wiyn.next_rising(ephem.Sun())
    # print "next sun set        ", wiyn.next_setting(ephem.Sun())
    # print "next sun rise (-18) ", wiyn_astrotwilight.next_rising(ephem.Sun())
    # print "next sun set (-18)  ", wiyn_astrotwilight.next_setting(ephem.Sun())
    # print "current             ", wiyn.date, float(wiyn.date)
    # print "time to sunset:     ", (wiyn.next_setting(ephem.Sun())-wiyn.date)*24.,"hours"
    # # print "time since sunrise: ", wiyn.date-wiyn.previous_rising(sun), (wiyn.date-wiyn.previous_rising(sun))*24.
    # # print "time to sunset:     ", (wiyn.next_setting(sun)-wiyn.date)*24.,"hours"
    # # print "time since sunrise: ", (wiyn.date-wiyn.previous_rising(sun))*24.,"hours"
    # print "sun altitude:       ", sun.alt, float(sun.alt), float(sun.az)
    # print

    #
    # Compute times since/until sunset and sunrise
    # Make sure to compute times to the closest sunrise and sunset, i.e. 
    # computations are different for the morning and afternoon
    #
    if (sun.alt > 0):
        # The sun is up, let's figure out if we are just beyond sunrise 
        # or just before sunset (i.e. it's morning or afternoon)
        time_since_last_noon = wiyn.previous_transit(ephem.Sun())
        time_to_next_noon = wiyn.next_transit(ephem.Sun())

        time_to_sunset = wiyn.next_setting(ephem.Sun())-wiyn.date
        time_since_sunrise = wiyn.date-wiyn.previous_rising(ephem.Sun())
        if (time_to_sunset < time_since_sunrise):
            # We are close to sunset
            hdr['SUN_RISE'] = (24.*(wiyn.next_rising(ephem.Sun())-wiyn.date), 
                                             'time until sunrise [hours]')
            hdr['SUN_SET'] = (24.*(wiyn.date-wiyn.next_setting(ephem.Sun())), 
                                            'time since sunset [hours]')
        else:
            # We are still before noon
            hdr['SUN_RISE'] = (24.*(wiyn.previous_rising(ephem.Sun())-wiyn.date), 
                                             'time until sunrise [hours]')
            hdr['SUN_SET'] = (24.*(wiyn.date-wiyn.previous_setting(ephem.Sun())), 
                                            'time since sunset [hours]')

    else:
        # Sun is below the horizon
        hdr['SUN_RISE'] = (24.*(wiyn.next_rising(ephem.Sun())-wiyn.date), 
                                         'time until sunrise [hours]')
        hdr['SUN_SET'] = (24.*(wiyn.previous_setting(ephem.Sun())-wiyn.date), 
                                        'time since sunset [hours]')


    #
    # Compute times since/until moon rise and set
    #
    if (moon.alt > 0):
        # Moon is up
        hdr['MOON_RSE'] = (24.*(wiyn.previous_rising(moon)-wiyn.date), 
                                         'time since moon rise [hours]')
        hdr['MOON_SET'] = (24.*(wiyn.next_setting(moon)-wiyn.date), 
                                         'time until moon set [hours]')

    else:
        # Moon is down
        hdr['MOON_RSE'] = (24.*(wiyn.next_rising(moon)-wiyn.date), 
                                         'time until moon rise [hours]')
        hdr['MOON_SET'] = (24.*(wiyn.previous_setting(moon)-wiyn.date), 
                                         'time since moon set [hours]')



    # time_since_ast_twilight = 24.*(wiyn.date-wiyn_astrotwilight.previous_setting(ephem.Sun()))
    # time_until_ast_twilight = 24.*(wiyn_astrotwilight.next_rising(ephem.Sun())-wiyn.date)
    # # if (numpy.degrees(sun.alt) > -18. and time_since_ast_twilight > 12.):
    # #     time_since_ast_twilight = 24.*(wiyn.date-wiyn_astrotwilight.next_setting(ephem.Sun()))
    # # if (numpy.degrees(sun.alt) > -18. and time_since_to_twilight > 16.):
    # #     time_until_ast_twilight = 24.*(wiyn_astrotwilight.previous_rising(ephem.Sun())-wiyn.date)

    # hdr['SUN_RS18'] = (time_until_ast_twilight,
    #                                  'time until astronomical twilight [hours]')
    # hdr['SUN_ST18'] = (time_since_ast_twilight,
    #                                  'time since astronomical twilight [hours]')

    #
    # Determine sky brightness
    #
    moon.compute(wiyn)
    sun.compute(wiyn)
    skybrite = None
    sun_alt = numpy.degrees(sun.alt)
    if (sun_alt < -18.):
        if (numpy.degrees(moon.alt) < -18):
            skybrite = 'dark'
        else:
            if (moon.moon_phase > 0.5):
                skybrite = "bright"
            else:
                skybrite = "grey"
    elif (sun_alt >= -18 and sun_alt < -12):
        skybrite = "astro.twilight"
    elif (sun_alt >= -12 and sun_alt < -6):
        skybrite = "naut.twilight"
    elif (sun_alt >= -6 and sun_alt < 0):
        skybrite = "civil.twilight"
    else:
        skybrite = 'daytime'
    hdr['SKYBRITE'] = (skybrite,
                                     "sky brightness quality")

    #
    # Compute lunar age, i.e. days since/until new moon
    #
    lunar_age = wiyn.date - ephem.previous_new_moon(wiyn.date)
    if (ephem.previous_full_moon(wiyn.date) > ephem.previous_new_moon(wiyn.date)):
        lunar_age = wiyn.date - ephem.next_new_moon(wiyn.date)
    hdr['LUNARAGE'] = (lunar_age,
                                     "days since last new moon")
    add_fits_header_title(hdr,
                          "Almanach data", 'SUN__RA')

    #for key in hdr:
    #    print "%-8s" % (key), " --> ", hdr[key]

    # hdr.toTxtFile("ephemdebug.head", clobber=True)

    #print wiyn.date+0.

if __name__ == "__main__":

    time_overwrite = None
    if (cmdline_arg_isset("-time")):
        time_overwrite = cmdline_arg_set_or_default("-time", "2014/14/20 14:00:00")

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)

    for filename in sys.argv[1:]:

        if (not os.path.isfile(filename)):
            continue

        try:
            hdulist = pyfits.open(filename)
            hdr = hdulist[0].header
        except:
            continue

        add_ephem_data_to_header(hdr, time_overwrite)



    podi_logging.shutdown_logging(options)
