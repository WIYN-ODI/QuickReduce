#! /usr/bin/env python
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

import sys
import os
import pyfits
import numpy
from astLib import astWCS
import scipy
import scipy.stats
import matplotlib.pyplot


def odi2wcs(xy, wcs):
    wcs.updateFromHeader()
    radec = numpy.ones((xy.shape[0],2)) * 1e9
    for i in range(xy.shape[0]):
        x, y = xy[i,0], xy[i,1]
        radec[i,0], radec[i,1] = wcs.pix2wcs(x, y)
    return radec


loopcounter = 0

def ferr(p, wcs, cat):

    global loopcounter

    # Set all WCS headers

    wcs.header['CRPIX1'] = p[0]
    wcs.header['CRPIX2'] = p[1]

    wcs.header['CRVAL1'] = p[2]
    wcs.header['CRVAL2'] = p[3]

    wcs.header['CD1_1'] = p[4]
    wcs.header['CD1_2'] = p[5]
    wcs.header['CD2_1'] = p[6]
    wcs.header['CD2_2'] = p[7]

    wcs.header['PV1_0'] = p[8]
    wcs.header['PV1_1'] = p[9]
    wcs.header['PV1_2'] = p[10]
    wcs.header['PV1_4'] = p[11]
    wcs.header['PV1_5'] = p[12]
    wcs.header['PV1_6'] = p[13]

    wcs.header['PV2_0'] = p[14]
    wcs.header['PV2_1'] = p[15]
    wcs.header['PV2_2'] = p[16]
    wcs.header['PV2_4'] = p[17]
    wcs.header['PV2_5'] = p[18]
    wcs.header['PV2_6'] = p[19]

    wcs.header['PV1_7'] = p[20]
    wcs.header['PV1_8'] = p[21]
    wcs.header['PV1_9'] = p[22]
    wcs.header['PV1_10']= p[23]

    wcs.header['PV2_7'] = p[24]
    wcs.header['PV2_8'] = p[25]
    wcs.header['PV2_9'] = p[26]
    wcs.header['PV2_10']= p[27]

    #pi = 8
    #for i in range(7):
    #    pv_name = "PV1_%d" % i
    #    wcs.header[pv_name] = p[pi]
    #    pi += 1

    #for i in range(7):
    #    pv_name = "PV2_%d" % i
    #    wcs.header[pv_name] = p[pi]
    #    pi += 1

    wcs.updateFromHeader()

    # Now convert all x/y to Ra/Dec

    #print cat.shape

    radec_ref = cat[:,2:4]
    radec_comp = odi2wcs(cat[:,0:2], wcs)

    d2 = (radec_ref - radec_comp) / (0.11/3600.)
    d_radec = numpy.sqrt(numpy.sum(d2**2, axis=1))
    #print d_radec.shape

    # filter out all signals exceeding 3x the 1-sigma deviation
    lowlimit, highlimit = -50, +50

    matching_radius = 50
    low_ra, high_ra, low_dec, high_dec = -matching_radius,matching_radius,-matching_radius,matching_radius
    good_matches = True
    for i in range(0):

        #new_good_matches = (d_radec > lowlimit) & (d_radec < highlimit) & (numpy.isfinite(d_radec))
        #new_good_matches = (d2[:,0] > low_ra) & (d2[:,0] < high_ra) & \
        #    (d2[:,1] > low_dec) & (d2[:,1] < high_dec) & (numpy.isfinite(d_radec))
        new_good_matches = (d_radec < highlimit) & (numpy.isfinite(d_radec))

        if (numpy.sum(new_good_matches) > len(p)):
            good_matches = new_good_matches
        else:
            good_matches = True
            break

        med = numpy.median(d_radec[good_matches])
        onesigma = scipy.stats.scoreatpercentile(d_radec[good_matches], 68)
        
        #lowlimit = #med - 1 * onesigma
        highlimit = 3 * onesigma

        #med = numpy.median(d_radec[good_matches])
        #onesigma = 0.5 * (scipy.stats.scoreatpercentile(d_radec[good_matches], 84) - scipy.stats.scoreatpercentile(d_radec[good_matches], 16))
        #lowlimit = med - 1 * onesigma
        #highlimit = med + 1 * onesigma

        #med_ra = numpy.median(d2[:,0][good_matches])
        #onesigma_ra = 0.5 * (scipy.stats.scoreatpercentile(d2[:,0][good_matches], 84) - scipy.stats.scoreatpercentile(d2[:,0][good_matches], 16))
        #low_ra = med_ra - 1 * onesigma_ra
        #high_ra = med_ra + 1 * onesigma_ra

        #med_dec = numpy.median(d_radec[good_matches])
        #onesigma_dec = 0.5 * (scipy.stats.scoreatpercentile(d2[:,1][good_matches], 84) - scipy.stats.scoreatpercentile(d2[:,1][good_matches], 16))
        #low_dec = med_dec - 1 * onesigma_dec
        #high_dec = med_dec + 1 * onesigma_dec

    med = numpy.median(d_radec)
    onesigma = scipy.stats.scoreatpercentile(d_radec, 68)
    highlimit = 2 * onesigma

    good_matches = d_radec < highlimit
    bad_matches  = d_radec > highlimit

    near_x = d2[:,0][good_matches]
    near_y = d2[:,1][good_matches]

    far_x = d2[:,0][bad_matches]
    far_y = d2[:,1][bad_matches]

    fig = matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(far_x, far_y, "r.", near_x, near_y, "bo")

    #matplotlib.pyplot.xlim((low_ra,high_ra))
    #matplotlib.pyplot.ylim((low_dec,high_dec))
    matplotlib.pyplot.xlim((-highlimit, highlimit))
    matplotlib.pyplot.ylim((-highlimit, highlimit))
    extname = wcs.header['EXTNAME']
    figname = "fit_%s_%d.png" % (extname, loopcounter)
    fig.savefig(figname)

    #print d_radec[0:5]
    sys.stdout.write("\rFitting %s, iteration %d" % (extname, loopcounter))
    sys.stdout.flush()

    loopcounter += 1

    if (onesigma < 5):
        d_radec[:] = 0

    filtered = good_matches & (numpy.isfinite(d_radec))
    return d_radec[filtered]
    
def save_get_key(hdr, key):

    try:
        return hdr[key]
    except:
        pass
    return 0

def get_init_value(wcs):

    p = [
        save_get_key(wcs.header, "CRPIX1"),
        save_get_key(wcs.header, "CRPIX2"),
        save_get_key(wcs.header, "CRVAL1"),
        save_get_key(wcs.header, "CRVAL2"),
        save_get_key(wcs.header, "CD1_1"),
        save_get_key(wcs.header, "CD1_2"),
        save_get_key(wcs.header, "CD2_1"),
        save_get_key(wcs.header, "CD2_2"),
        save_get_key(wcs.header, "PV1_0"),
        save_get_key(wcs.header, "PV1_1"),
        save_get_key(wcs.header, "PV1_2"),
        #save_get_key(wcs.header, "PV1_3"),
        save_get_key(wcs.header, "PV1_4"),
        save_get_key(wcs.header, "PV1_5"),
        save_get_key(wcs.header, "PV1_6"),
        save_get_key(wcs.header, "PV2_0"),
        save_get_key(wcs.header, "PV2_1"),
        save_get_key(wcs.header, "PV2_2"),
        #save_get_key(wcs.header, "PV2_3"),
        save_get_key(wcs.header, "PV2_4"),
        save_get_key(wcs.header, "PV2_5"),
        save_get_key(wcs.header, "PV2_6"),
        \
        save_get_key(wcs.header, "PV1_7"),
        save_get_key(wcs.header, "PV1_8"),
        save_get_key(wcs.header, "PV1_9"),
        save_get_key(wcs.header, "PV1_10"),
        \
        save_get_key(wcs.header, "PV2_7"),
        save_get_key(wcs.header, "PV2_8"),
        save_get_key(wcs.header, "PV2_9"),
        save_get_key(wcs.header, "PV2_10"),
        ]
        
    return p



if __name__ == "__main__":
    global loopcounter

    filename = sys.argv[1]
    print filename

    hdulist = pyfits.open(filename)

    for extension in range(1, len(hdulist)):
        loopcounter = 0

        extname = hdulist[extension].header['EXTNAME']
        if (extname[0:3] != "OTA"):
            continue

        ota = int(extname[3:5])

        catfile = "matchcat_%02d.dat" % ota
        cat = numpy.loadtxt(catfile)
        
        wcs = astWCS.WCS(hdulist[extension].header, mode="pyfits")
        pinit = get_init_value(wcs)
        
        radec_ref = cat[:,2:4]
        radec_comp = odi2wcs(cat[:,0:2], wcs)
        d2 = (radec_ref - radec_comp) / (0.11/3600.)
        d_radec = numpy.sqrt(numpy.sum(d2**2, axis=1))
        cat = cat[d_radec < 30]

        errfunc = lambda p, wcs, cat: ferr(p, wcs, cat)

        #ferr(pinit, wcs, cat)

        out = scipy.optimize.leastsq(errfunc, pinit,
                                     args=(wcs, cat)) #, full_output=1)

    
        hdulist[extension].header = wcs.header

        print "done (%d)" % (loopcounter)
         #   d_crval1, d_crval2 = out[0][0]-pinit[0], out[0][1]-pinit[1]
         #   ota_list[ext].header.update("D_CRVAL1", d_crval1)
         #   ota_list[ext].header.update("D_CRVAL2", d_crval2)

    outfile = sys.argv[2]
    hdulist.writeto(outfile, clobber=True)
