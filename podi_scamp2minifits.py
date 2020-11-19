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
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import astropy.io.fits as pyfits
import numpy
import scipy
import ephem
from astLib import astWCS
import dev_ccmatch

gain_correct_frames = False
from podi_definitions import *

def scamp_header_to_minifits(filename, minifits_outputname, reference_fits, 
                             reference_extension="OTA33.SCI", distortion_level=4,
                             recenter_odi=True):

    if (not os.path.isfile(filename)):
        return None

    headfile = open(filename, "r")
    lines = headfile.readlines()

    values = {}
    values_list = []
    comments = {}
    comments_list = []

    extlist = []
    primhdu = pyfits.PrimaryHDU()
    extlist.append(primhdu)
    
    for line in lines:
        #print line

        key, value, comment = line[0:8].strip(), line[9:30].strip(), line[32:].strip()
        if (key in ("HISTORY", "COMMENT",
                    ) ):
            # Don't know what to do with those, so skip'em
            continue

        elif (key in ("FGROUPNO", "FLXSCALE", "MAGZEROP", 
                      "ASTINST",
                      "PHOTIRMS", "PHOTINST", "PHOTLINK",
                      ) ):
            # These are some scamp-specific headers, let's not copy them
            continue

        elif (key in ("CRVAL1", "CRVAL2",
                      "CRPIX1", "CRPIX2", 
                      "CD1_1", "CD1_2", "CD2_1", "CD2_2",
                      "EQUINOX",
                      ) or (key[0:2] == "PV" and key[3] == "_")):
            value = float(value)

        elif (key in ("RADECSYS", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2") ):
            # Strip the unnecessary quotes and spaces from the end
            value = value[1:-1].strip()

        elif (key == "END"):
            # This concludes one extension, add it to list and start new 
            # list for the next OTA

            values_list.append(values)
            comments_list.append(comments)
            values = {}
            comments = {}
            continue

        values[key] = value
        comments[key] = comment

    refhdu = pyfits.open(reference_fits)
    # Count how many image extensions exist in the reference frame
    n_imageext = 0
    ota_names = []
    for i in range(len(refhdu)):
        if (is_image_extension(refhdu[i])):
            n_imageext += 1
            ota_names.append(refhdu[i].name)
        # try:
        #     extname = refhdu[i].header['EXTNAME']
        #     if (extname[0:3] == "OTA" and extname[-3:] == "SCI"):
        #         n_imageext += 1
        # except:
        #     pass

    if (len(values_list) != n_imageext):
        stdout_write("Illegal scamp solution or wrong reference fits (%d vs %d)!\n" % (len(values_list), len(refhdu)-1))
        return -1

    # Now copy all scamp headers into the minifits
    # With this data at hand, work out the shift we need to apply to the scamp solution
    print("Creating the minfits headers")
    for ota in range(len(values_list)):

        # Create a new Image extension to hold the header
        hdu = pyfits.ImageHDU(name=ota_names[ota])

        # Copy all headers to the new HDU
        for key in values_list[ota]:
            hdu.header[key] = (values_list[ota][key], comments_list[ota][key])

        for key in ['NAXIS', 'NAXIS1', 'NAXIS2']:
            hdu.header[key] = refhdu[ota_names[ota]].header[key]

        # And add the HDU to the list that will form the minifits
        extlist.append(hdu)


    #
    # Work out the position of the origin -
    # take the center between the 33/34/43/44 OTAs as optical center
    #
    tmp_hdulist = pyfits.HDUList(extlist)

    if (recenter_odi):
        center_pos = {
            'OTA33.SCI': [4096, 4096],
            'OTA34.SCI': [4096,    0],
            'OTA43.SCI': [   0, 4096],
            'OTA44.SCI': [   0,    0],
        }
        centers = []
        for ota in center_pos:
            xy = center_pos[ota]
            ext = tmp_hdulist[ota]
            print(xy, ext.header)

            wcs = astWCS.WCS(ext.header, mode='pyfits')
            ra_dec = wcs.pix2wcs(xy[0], xy[1])
            print(ra_dec)
            centers.append(ra_dec)
        centers = numpy.array(centers)
        numpy.savetxt(sys.stdout, centers)
        optical_center = numpy.mean(centers, axis=0)
        print(optical_center)
        print(numpy.std(centers, axis=0))
        print(numpy.std(centers, axis=0)*3600.)
    else:
        wcs = astWCS.WCS(refhdu[reference_extension].header, mode='pyfits')
        optical_center = numpy.array(wcs.getCentreWCSCoords())

    #
    # Now modify and re-fit the WCS with the correct reference point
    #
    n_stars = 5000
    tmp_hdulist.info()
    output_hdulist = [pyfits.PrimaryHDU()]

    for ota in tmp_hdulist[1:]:
        print(ota.name)

        # create a number of random points scattered throughout the image
        xy = numpy.random.random((n_stars,2))
        xy *= [ota.header['NAXIS1'], ota.header['NAXIS2']]
        print(xy.shape)
        print(xy[:5])

        wcs = astWCS.WCS(ota.header, mode='pyfits')
        ra_dec = numpy.array(wcs.pix2wcs(xy[:,0], xy[:,1]))

        # Now create the catalog we need to re-compute the distortion parameters
        # Columns of catalog to be compatible with optimize_wcs_solution
        # 3/4: x/y
        # last 2 columns: reference ra/dec
        catalog = numpy.zeros((n_stars, 6))
        catalog[:,0:2] = ra_dec[:,:]
        catalog[:,2:4] = xy[:,:]
        catalog[:,4:6] = ra_dec[:,:]


        new_hdu = pyfits.ImageHDU(header=ota.header, name=ota.name)

        # define the new optical axis
        new_hdu.header['CRVAL1'] = optical_center[0]
        new_hdu.header['CRVAL2'] = optical_center[1]

        # reset some of the headers we do not need
        for key in ['PV1_0', 'PV1_2', 'PV2_0', 'PV2_2']:
            new_hdu.header[key] = 0.
        for key in ['PV1_1', 'PV2_1']:
            new_hdu.header[key] = 1.

        # set image size - otherwise WCS transformation won't work
        for key in ['NAXIS', 'NAXIS1', 'NAXIS2']:
            new_hdu.header[key] = ota.header[key]

        header_keywords  = ['CRPIX1', 'CRPIX2',
                            #'CRVAL1', 'CRVAL2',
                            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',]
        if (distortion_level >= 2):
            header_keywords.extend([
                'PV1_4', 'PV1_5', 'PV1_6', 'PV2_4', 'PV2_5', 'PV2_6',
            ])
        if (distortion_level >= 3):
            header_keywords.extend([
                'PV1_7', 'PV1_8', 'PV1_9', 'PV1_10',
                'PV2_7', 'PV2_8', 'PV2_9', 'PV2_10',
            ])
        if (distortion_level >= 4):
            header_keywords.extend([
                'PV1_12', 'PV1_13', 'PV1_14', 'PV1_15', 'PV1_16',
                'PV2_12', 'PV2_13', 'PV2_14', 'PV2_15', 'PV2_16',
            ])
        if (distortion_level >= 5):
            header_keywords.extend([
                'PV1_17', 'PV1_18', 'PV1_19', 'PV1_20', 'PV1_21', 'PV1_22',
                'PV2_17', 'PV2_18', 'PV2_19', 'PV2_20', 'PV2_21', 'PV2_22',
            ])

        #print header_keywords
        #print new_hdu.header

        print("re-fitting & optimizing WCS solution for OTA %s" % (ota.name))
        dev_ccmatch.optimize_wcs_solution(catalog, new_hdu.header, header_keywords)

        # Now re-compute the Ra/Dec positions using the new WCS system
        wcs_out = astWCS.WCS(new_hdu.header, mode='pyfits')
        radec_new = numpy.array(wcs_out.pix2wcs(xy[:,0]-1, xy[:,1]-1))
        catalog[:,0:2] = radec_new[:,:]
        numpy.savetxt("recomputed.%s" % (ota.name), catalog)

        # compute RMS value of old-new
        d_radec = catalog[:, 0:2] - catalog[:, 4:6]
        rms_radec = numpy.std(d_radec, axis=0) * 3600
        print("RMS of fit: %f/%f arcsec" % (rms_radec[0], rms_radec[1]))

        # Finally, delete the CRVAL header keywords and prepare the new ImageHDU
        for key in ['CRVAL1', 'CRVAL2']:
            new_hdu.header[key] = 0.

        # fix the RADESYS header
        new_hdu.header['RADESYS'] = 'ICRS'

        output_hdulist.append(new_hdu)

        # break

    # Create a proper HDUList from the list of new HDUs
    minifits_hdulist = pyfits.HDUList(output_hdulist)

    # return -1

    # Now go through both the minifits and the reference and assign the extname keywords
    # print "Correcting CRVALs for reference position"
    # ota33_found = False
    # for ext in range(1, len(refhdu)):
    #     if (not is_image_extension(refhdu[ext])):
    #         continue
    #     extname = refhdu[ext].header['EXTNAME']
    #     if (extname == reference_extension): ota33_found = True
    #     extlist[ext].name = extname

    # # Create a proper HDUList from the list of new HDUs
    # minifits_hdulist = pyfits.HDUList(extlist)


    # Next (and almost finally) we need to change the CRVAL keywords to account for the
    # offsets in pointing center and the center assumed by scamp

    # Adjust the solution for the known shift in CRPIX position
    # In pODI: True pointing of telescope is ~ at pixel 4200,4200 of OTA 3,3
    # From scamp: Typically somewhere around 2000,2000 in OTA 3,3
    # NB: 0.11 / 3600. is the pixel scale in degrees/pixel, ignoring rotation and distortion
    # if (ota33_found):
    #     dx_crpix = 4200 - minifits_hdulist[reference_extension].header["CRPIX1"]
    #     dy_crpix = 4200 - minifits_hdulist[reference_extension].header["CRPIX2"]
    #     d_ra  = dx_crpix * 0.11 / 3600.
    #     d_dec = dy_crpix * 0.11 / 3600.
    # else:
    #     d_ra, d_dec = 0, 0
    #
    # # fix potential problems with d_ra > 360 or d_ra < 0
    # d_ra = d_ra - math.floor(d_ra/360.)*360
    #
    # for ext in range(1, len(minifits_hdulist)):
    #     minifits_hdulist[ext].header['CRVAL1'] = d_ra
    #     minifits_hdulist[ext].header['CRVAL2'] = d_dec

    # Last step: write the minifits to file
    minifits_hdulist.writeto(minifits_outputname, overwrite=True)
    
    # That's it, all done!
    return 0


if __name__ == "__main__":

    scampfile = sys.argv[1]

    outputfile = sys.argv[2]

    reference = sys.argv[3]

    try:
        ref_ext = sys.argv[4]
    except:
        ref_ext = "OTA33.SCI"

    recenter = True
    try:
        if (sys.argv[5] == "simple"):
            recenter = False
    except:
        pass

    print(recenter)

    scamp_header_to_minifits(scampfile, outputfile, reference, ref_ext, recenter_odi=recenter)
    
