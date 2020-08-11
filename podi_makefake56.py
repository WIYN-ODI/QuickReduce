#!/usr/bin/env python3

import sys, os
import astropy.io.fits as pyfits
import datetime
import warnings
warnings.filterwarnings('ignore')
import podi_focalplanelayout

old_new = [(00,55),
           (16,25),
           (22,26),
           (23,14),
           (24,11),
           (32,12),
           (33,45),
           (34,13),
           (42,46),
           (43,36),
           (44,35),
           (55,16),
           (61,15),

           (33,21), (33,22), (33,23), (33,24),
           (33,31), (33,32), (33,33), (33,34),
           (33,41), (33,42), (33,43), (33,44),
           (33,51), (33,52), (33,53), (33,54),
           (33,56),
]

if __name__ == "__main__":

    outpath = sys.argv[1]

    # Now rename a bunch of files to reflect their position in the upgraded focalplane

    for dirname in sys.argv[2:]:
        #
        #
        #
        _, obsid = os.path.split(os.path.abspath(dirname))
        print "Working on ",obsid

        exptype = obsid[0]
        # print "exptype=",exptype

        for (old,new) in old_new:
            # print old,"-->",new, obsid

            _pos = obsid.index(".")
            dither_number = int(obsid[_pos+1:])
            # print dither_number

            infilename = "%s/%s.%02d.fits" % (obsid, obsid,old)
            hdulist = pyfits.open(infilename)

            hdulist[0].header['FPPOS'] = 'xy%02d' % (new)

            # chaneg all observing date forward by 3 years
            date_format = "%Y-%m-%dT%H:%M:%S.%f"
            date_obs = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], date_format)

            date_new = date_obs + datetime.timedelta(days=3*365)
            hdulist[0].header['DATE-OBS'] = date_new.strftime(date_format)
            # print date_new.strftime(date_format)

            #print hdulist[0].header['DATE-OBS'], date_obs, date_new

            hdulist[0].header['MJD-OBS'] += 3*365 # change time
            
            fpl = podi_focalplanelayout.FocalPlaneLayout(hdulist[0])
            
            obsid_timefmt = "%Y%m%dT%H%M%S"
            #localtime = date_new - timedelta(hours=7)
            #new_obsid = "%s.%d" % (date_new.strftime(obsid_timefmt), dither_number)
            #print obsid
            new_obsid = "%s%04d%s" % (obsid[0], int(obsid[1:5])+3, obsid[5:])
            #print new_obsid

            # Change the OTA-ID to simulate a more realistic pattern for broken cells
            hdulist[0].header["OTA_ID"] = "%d" % (fpl.get_otaid_from_position(new))
            hdulist[0].header['OBSID'] = new_obsid[1:]

            hdulist[0].header['FILENAME'] = "%s.%02d.fits" % (new_obsid, new)

            outdir = "%s/%s" % (outpath, new_obsid)
            if (not os.path.isdir(outdir)):
                os.mkdir(outdir)

            outfilename = "%s/%s.%02d.fits" % (outdir, new_obsid, new)
            if (not os.path.isfile(outfilename)):
                print "   writing output:",outfilename
                hdulist.writeto(outfilename, overwrite=True)
        
