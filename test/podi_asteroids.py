#!/usr/bin/env python

import pyfits
from astLib import astWCS

import numpy

import astropy
import astropy.io
import astropy.io.votable
import astropy.wcs

import ephem

import os
import sys

mjd0 = + 2400000.5

def get_minorobjects(epoch, ra, dec, sr):


    # now query skybot
    query = "http://vo.imcce.fr/webservices/skybot/conesearch?EPOCH=%(epoch)f&RA=%(ra)f&DEC=%(dec)f&SR=%(sr)f&VERB=3"


    import urllib
    params = urllib.urlencode({'EPOCH': epoch, 
                               'RA': ra,
                               'DEC': dec,
                               'SR': sr,
                               'VERB': 3,
                           })
    url = 'http://vo.imcce.fr/webservices/skybot/conesearch?'
    srccat = urllib.urlopen(url+'?%s' % params).read()

    # print srccat
    file = open("test.xml", "w")
    file.write(srccat)
    file.close()

    try:
        votab = astropy.io.votable.parse_single_table('test.xml')
    except IndexError:
        print "no table found"
        return None

    # print votab.array


    # print [x.name for x in votab.fields]  

    return votab



def create_ds9_region_file(votab, region_file, epoch, mjd_start, mjd_end, rel_nonsidereal):

    obj_ra = votab.array['RA']
    obj_dec = votab.array['DEC']
    names = votab.array['Name']

    exptime = (mjd_end - mjd_start) * 86400.

    pm_ra = votab.array['dracosdec'] - rel_nonsidereal[0]
    pm_dec = votab.array['ddec'] - rel_nonsidereal[1]

    # Now create a ds9 region file with all information

    # coord_j2000 = ephem.Equatorial(obj_ra, obj_dec, epoch=ephem.J2000)

    year_2014 = 2456658.500000
    
    ds9 = open(region_file, 'w')
    print >>ds9, 'global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
    print >>ds9, 'fk5'

    for obj in range(obj_ra.shape[0]):

        coord_j2000 = ephem.Equatorial(obj_ra[obj], obj_dec[obj], epoch=ephem.J2000)
        print "  %-13s : %s %s" % (names[obj], coord_j2000.ra, coord_j2000.dec)
        length = numpy.hypot(pm_ra[obj], pm_dec[obj]) * (exptime/3600.) / 3600.
        angle = numpy.degrees(numpy.arctan2(pm_ra[obj], pm_dec[obj])) + 90.

        offset_ra = (epoch - mjd_start)*24.*pm_ra[obj] if mjd_start > 0 else 0.
        offset_dec = (epoch - mjd_start)*24.*pm_dec[obj] if mjd_start > 0 else 0.

        
        dict = {
            'ra': numpy.degrees(coord_j2000.ra)-offset_ra,
            'dec': numpy.degrees(coord_j2000.dec)-offset_ra,
            'pm_ra': pm_ra[obj] * (exptime/3600.) / 3600.,
            'pm_dec': pm_dec[obj] * (exptime/3600.) / 3600.,
            'length': length,
            'angle': angle,
            'name': names[obj],
        }

        print >>ds9, '# vector(%(ra)f, %(dec)f, %(length)f, %(angle)f) vector=1' % dict
        print >>ds9, '# text(%(ra)f,%(dec)f) text={%(name)s}' % dict

    ds9.close()



def get_asteroid_list_from_fitsfile(fitsfile):

    hdulist = astropy.io.fits.open(fitsfile)

    wcs = astropy.wcs.WCS(hdulist[0].header)
    print wcs.wcs.name
    # print wcs.wcs.print_contents()


    # get corner coordinates
    naxis1, naxis2 = hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2']

    centerx, centery = naxis1/2, naxis2/2
    center = numpy.array([[centerx, centery]])
    print center

    center_coord = wcs.wcs_pix2world(center, 1)

    print center_coord

    epoch = hdulist[0].header['MJD-OBS'] + mjd0

    ra = center_coord[0,0]
    dec = center_coord[0,1]
    sr = 0.4

    votab = get_minorobjects(epoch, ra, dec, sr)

    hdulist.close()

    return votab

if __name__ == "__main__":



    fitsfile = sys.argv[1]

    votab = get_asteroid_list_from_fitsfile(fitsfile)

    hdulist = astropy.io.fits.open(fitsfile)
    rel_nonsidereal_ra = 0.
    rel_nonsidereal_dec = 0.
    if ('_NSID_RA' in hdulist[0].header):
        rel_nonsidereal_ra = hdulist[0].header['_NSID_RA']
    if ('_NSID_DE' in hdulist[0].header):
        rel_nonsidereal_dec = hdulist[0].header['_NSID_DE']


    epoch = hdulist[0].header['MJD-OBS'] + mjd0
    mjd_start = hdulist[0].header['MJD-STRT']+mjd0 if ('MJD-STRT' in hdulist[0].header) else 0
    mjd_end   = hdulist[0].header['MJD-END']+mjd0  if ('MJD-END'  in hdulist[0].header) else 0

    exptime = 3500.

    obj_ra = votab.array['RA']
    obj_dec = votab.array['DEC']
    names = votab.array['Name']

    pm_ra = votab.array['dracosdec'] - rel_nonsidereal_ra
    pm_dec = votab.array['ddec'] - rel_nonsidereal_dec

    # Now create a ds9 region file with all information

    # coord_j2000 = ephem.Equatorial(obj_ra, obj_dec, epoch=ephem.J2000)

    region_file = sys.argv[2]

    create_ds9_region_file(votab, region_file, epoch, mjd_start, mjd_end, [rel_nonsidereal_ra, rel_nonsidereal_dec])




