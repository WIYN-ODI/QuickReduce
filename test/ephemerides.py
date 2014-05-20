#!/usr/bin/env python

import numpy
import scipy
import sys,os
import ephem
import math

import matplotlib
import matplotlib.pyplot
import datetime
import scipy.interpolate


def load_ephemerides(filename, plot=False):

    # Load file
    dat = open(filename, "r")
    lines = dat.readlines()
    # print lines[:2]

    ref_date_str = "2014-Jan-01 00:00"
    date_format = "%Y-%b-%d %H:%M"
    ref_date = datetime.datetime.strptime(ref_date_str, date_format)
    # print int(ref_date)
    ref_mjd = 56658.

    data = []
    for line in lines:
#        print line.split()
        items = line.split()

        date_time = items[0:2]

        coords = 'deg'
        if (coords == 'deg'):
            ra = float(items[2])
            dec = float(items[3])
            rate_ra = float(items[4])
            rate_dec = float(items[5])
        else:
            str_ra = items[2:5]
            str_dec = items[5:8]
            coord = ephem.Equatorial(" ".join(str_ra), " ".join(str_dec), epoch=ephem.J2000)
            ra = math.degrees(coord.ra)
            dec = math.degrees(coord.dec)
            rate_ra = float(items[8])
            rate_dec = float(items[9])

        mjd = 0

        # print date_time
        date_obs = datetime.datetime.strptime(" ".join(date_time), 
                                              date_format)

        delta_t = date_obs - ref_date
        mjd = (delta_t.total_seconds()/86400.) + ref_mjd

        data.append( [mjd, 
                      ra,
                      dec,
                      rate_ra, 
                      rate_dec]
                     )

    data = numpy.array(data)

    ra_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,1], kind='linear' )
    dec_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,2], kind='linear' )

    if (plot):
        fig = matplotlib.pyplot.figure()
        ax = fig.add_subplot(111)

        # #ax.plot(data[:,3], data[:,4])
        # # ax.plot(data[:,1], data[:,2])
        ax.plot(data[:,0], data[:,1])

        x = numpy.linspace(data[0,0], data[-1,0], 20)
        ax.scatter(x, dec_vs_mjd(x))

        fig.show()
        matplotlib.pyplot.show()

    return ra_vs_mjd, dec_vs_mjd, data


if __name__ == "__main__":
    filename = sys.argv[1]
    load_ephemerides(filename, plot=True)



    
