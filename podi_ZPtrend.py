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


try:
    import cPickle as pickle
except:
    import pickle

if __name__ == "__main__":

    
    obstype, exptime, filtername, photzp, photzpe, mjd, dateobs = [], [], [], [], [], [], []

    print "Reading data"

    # Open the file, read its content, and add to the existing filelist
    pickled_file = "index.pickle"
    direntry = {}
    if (os.path.isfile(pickled_file)):
        try:
            pickle_dict = open(pickled_file, "rb")
            print "Reading pickled file..."
            direntry = pickle.load(pickle_dict)
            close(pickle_dict)
        except:
            pass


    for filename in sys.argv[1:]:

        if (not os.path.isfile(filename)):
            continue

        if (filename in direntry):
            # We know about this file from the pickle
            (_obstype, _exptime, _filtername, _photzp, _photzpe, _mjdobs, _dateobs) = direntry[filename]

        else:
            "Getting info for file", filename
            try:
                hdulist = pyfits.open(filename)
                hdr = hdulist[0].header
            except:
                continue


            _obstype, _exptime, _filtername, _photzp, _photzpe, _mjd, _dateobs = "???", 0, "???", -999, -999, 0, "..."

            try:
                _obstype = hdr['OBSTYPE']
                _exptime = hdr['EXPTIME']
                _filtername = hdr['FILTER']
                _photzp = hdr['PHOTZP']
                _photzpe = hdr['PHOTZPE']
                _mjdobs = hdr['MJD-OBS']
                _dateobs = hdr['DATE-OBS']
                direntry[filename] = (_obstype, _exptime, _filtername, _photzp, _photzpe, _mjdobs, _dateobs)
            except:
                pass



        obstype.append(_obstype)
        exptime.append(_exptime)
        filtername.append(_filtername)
        photzp.append(_photzp)
        photzpe.append(_photzpe)
        mjd.append(_mjdobs)
        dateobs.append(_dateobs)


    # Now pickle the dir_entry for later use
    picklejar = "index.pickle"
    with open(picklejar, "wb") as pf:
        pickle.dump(direntry, pf)
        pf.close()


    # Now create the plots
    known_filters = {
        "odi_g": (26.1, "blue"),
        "odi_r": (26.1, "green"),
        "odi_i": (27.0, "orange"),
        "odi_z": (25.0, "red"),
    }


    fig, ax = matplotlib.pyplot.subplots()
    #tfig, tax = matplotlib.pyplot.subplots()

    cc = matplotlib.colors.ColorConverter()

    # This is the MJD of 01/01/0001
    mjd_zeropoint = 1721424.500000 - 2400000.5 + (7./24.0)

    def dzp_to_transparency(d_zp):
        return 100.*numpy.power(10., 0.4*d_zp)

    print "Plotting"
    for thisfilter in set(filtername):

        print thisfilter

        legendname = None
        if (thisfilter in known_filters):
            legendname = thisfilter

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
            colors.append(cc.to_rgba(this_color))
            #colors.append(this_color) #cc.to_rgba(this_color))

        data = numpy.array(data)
        timestamp = data[:,0] + 0.5 * data[:,1] / 86400.  - mjd_zeropoint# mjdobs + 0.5*exptime
        
        print colors
        print list(colors)
        print "\n\n\n\n"
        ax.errorbar(x=timestamp, y=data[:,4], xerr = 0.5*data[:,1]/86400., yerr=data[:,3],
#                    color=[colors[i] for i in range(len(colors))],
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

    ax.set_ylabel("ZP difference Ref-ODI [mag]")
    #ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator())
    ax.xaxis.set_major_formatter(hfmt)
    ax.set_ylim(bottom = 0)
    matplotlib.pyplot.xticks(rotation='vertical')
    #ax.get_xaxis.set_xticks(rotation='vertical')
    matplotlib.pyplot.subplots_adjust(bottom=.20, top=0.95, right=0.9)

    #tax.set_ylabel("Transparency")
    #tax.xaxis.set_major_formatter(hfmt)

    # ax2 = ax.twinx()
    # ax2.set_ylabel("Transparency")
    # matplotlib.pyplot.yticks((0,-1,-2.5), (100, 40, 10))

    #ax2.set_yscale('log')

    time_start = numpy.min(mjd)-mjd_zeropoint
    time_end = numpy.max(mjd)-mjd_zeropoint
    frac_diff = 0.05 * (time_end - time_start)
    time_start, time_end = time_start - frac_diff, time_end + frac_diff
    matplotlib.pyplot.hlines(0,time_start, time_end)


    ax.set_xlim((time_start, time_end))
    #tax.set_xlim((time_start, time_end))

    ax.set_ylim((-5,0.3))
    ax.legend(loc='best', borderaxespad=1)

    #tax.set_ylim((1,200))
    #tax.set_yscale('log')
    #tax.legend(loc='best', borderaxespad=1)

    fig.savefig("photzp_trend.png")
    #tfig.savefig("transparency_trend.png")

    matplotlib.pyplot.show()


