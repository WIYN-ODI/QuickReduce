#!/usr/bin/env python3


import os
import sys
import numpy
import astropy.io.fits as pyfits
from astLib import astWCS
import ephem


def get_size(s):

    if s.endswith('"'):
        return float(s[:-1])/3600.
    elif s.endswith("'"):
        return float(s[:-1])/60.

    return float(s)

if __name__ == "__main__":

    region_fn = sys.argv[1]
    cutout_dir = sys.argv[2]
    filelist = sys.argv[3:]


    #
    # Load the region file
    #
    regions = []
    with open(region_fn, "r") as reg:
        lines = reg.readlines()
        for line in lines:
            if (line.startswith("box(")):
                items = line.split('box(')[1].split(')')[0].split(",")
                name = line.split('text={')[1].split('}')[0]
                print(items, name)

                if (items[0].find(":") > 0 or items[1].find(":") > 0):
                    coord = ephem.Equatorial(items[0], items[1])
                    ra = numpy.degrees(coord.ra)
                    dec = numpy.degrees(coord.dec)
                else:
                    ra = float(items[0])
                    dec = float(items[1])

                cos_dec = numpy.cos(numpy.radians(dec))
                width = get_size(items[2])/2. / cos_dec
                height = get_size(items[3])/2.

                # print(ra, dec, width, height, cos_dec)

                # now calculate the corner positions
                corners = numpy.array([
                    [ra-width, dec-height],
                    [ra+width, dec-height],
                    [ra-width, dec+height],
                    [ra+width, dec+height],
                ])
                # print(corners)

                regions.append((name, corners))

    for fn in filelist:

        print("Working on",fn)
        _,bn = os.path.split(fn)
        hdulist = pyfits.open(fn)

        # for now only use primary ext
        ext = hdulist[0]

        wcs = astWCS.WCS(ext.header, mode='pyfits')

        for name, corners in regions:
            data = ext.data

            xy = numpy.array(wcs.wcs2pix(corners[:,0], corners[:,1]))

            x1 = int(numpy.floor(numpy.max([0., numpy.min(xy[:,0])])))
            x2 = int(numpy.ceil(numpy.min([data.shape[1], numpy.max(xy[:,0])])))

            y1 = int(numpy.floor(numpy.max([0., numpy.min(xy[:,1])])))
            y2 = int(numpy.ceil(numpy.min([data.shape[0], numpy.max(xy[:,1])])))

            # print(xy, x1, x2, y1, y2)


            cutout = data[y1:y2+1, x1:x2+1]
            out_hdu = pyfits.PrimaryHDU(data=cutout, header=ext.header)
            out_hdu.header['CRPIX1'] -= x1
            out_hdu.header['CRPIX2'] -= y1

            out_fn = os.path.join(cutout_dir, bn[:-5]+".%s.fits"%(name))
            print("Writing",out_fn)
            out_hdu.writeto(out_fn, overwrite=True)
