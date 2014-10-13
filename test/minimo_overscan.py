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


import sys
import os
import signal
import pyfits
import numpy
import scipy
import scipy.optimize
import ephem
import traceback
#import psutil
import datetime
import warnings

import Queue
import threading
import multiprocessing
import multiprocessing.reduction
import ctypes
import time
import logging
import itertools

from astLib import astWCS

from podi_definitions import *
from podi_commandline import *
import podi_imcombine
import podi_imarith
import podi_sitesetup as sitesetup
import podi_logging


def minimo_filter(ext):
    f = ext.header['FILTER']
    f = f.replace(" ", "_")
    f = f.replace("/", "_")
    return f

def merge_amps(hdulist):

    # Create new data buffer
    x1,x2,y1,y2 = break_region_string(hdulist[0].header['CCDSIZE'])
    buffer = numpy.zeros((y2-y1+1, x2-x1+1), dtype=numpy.float32)
    # print buffer.shape

    # for i in range(2):
    #     print hdulist[i].header['DATASEC'], hdulist[i].header['CCDSEC']
    #     insert_into_array(hdulist[i].data, '[1:1024, 1:4096]', #hdulist[i].header['DATASEC'],
    #                       buffer, hdulist[i].header['CCDSEC'])

    buffer[0:4096, 0:1024] = hdulist[0].data[:,:]
    buffer[0:4096, 1024:2048] = hdulist[1].data[:,:]

    imghdu = pyfits.ImageHDU(header=hdulist[0].header, data=buffer)
    for hdr in ['CCDSIZE',
                'CCDSUM',
                'CCDSEC',
                'AMPSEC',
                'DATASEC',
                'DETSEC',
                'BIASSEC',
                'TRIMSEC',
                'ATM1_1',
                'ATM2_2',
                'ATV1',
                'ATV2',
                'LTM1_1',
                'LTM2_2',
                'LTV1',
                'LTV2',
                'DTM1_1',
                'DTM2_2',
                'DTV1',
                'DTV2']:
        if (hdr in imghdu.header):
            del imghdu.header[hdr]

    # print imghdu.header

    # Try to fix the un-FITS-conform CD header 
    for hdr in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
        dum = float(imghdu.header[hdr])
        del imghdu.header[hdr]
        imghdu.header[hdr] = dum

    # cd11 = float(imghdu.header['CD1_1'])
    # print cd11
    # del imghdu.header['CD1_1']
    # imghdu.header['CD1_1'] = cd11

#    imghdu.header['CD1_1'] = float(imghdu.header['CD1_1'])

    return imghdu


def overscan_subtract(filename, outfile):

    hdulist = pyfits.open(inputfile)

    for ext in hdulist[1:]:

        biassec = extract_region(ext.data, ext.header['BIASSEC'])

        datasec = extract_region(ext.data, ext.header['DATASEC'])

        # Compute line-by-line bias-sec
        overscan = numpy.median(biassec, axis=1).reshape((-1,1))
        # print biassec.shape, overscan.shape, ext.header['DATASEC']

        data = datasec - overscan
        # Re-insert the corrected data
        ext.data = data
        # insert_into_array(data, ext.header['DATASEC'],
        #                   ext.data, ext.header['DATASEC'])

    # Now we have all data overscan-subtracted
    # Next we'll merge the 2 amplifiers of each chip

    out_hdulist = [hdulist[0]]
    # print "merged amps"
    out_hdulist.append(merge_amps(hdulist[1:3]))
    out_hdulist.append(merge_amps(hdulist[3:5]))

    out_hdulist = pyfits.HDUList(out_hdulist)

    out_hdulist[1].name = 'CCD01.SCI'
    out_hdulist[1].header['DETSEC'] = '[1:2048, 1:4096]'
    out_hdulist[1].header['OTA'] = 01

    out_hdulist[2].name = 'CCD02.SCI'
    out_hdulist[2].header['DETSEC'] = '[2101:4148, 1:4096]'
    out_hdulist[2].header['OTA'] = 02
    
    out_hdulist.writeto(outfile, clobber=True)



def make_calibrations(filelist):
    
    biaslist = []
    filterlist = []
    flatlist = {}
    
    for filename in filelist:
        hdulist = pyfits.open(filename)
        
        obstype = hdulist[0].header['OBSTYPE']
        filter = minimo_filter(hdulist[0])
        
        if (obstype == 'dome flat'):
            if (not filter in filterlist):
                filterlist.append(filter)
                flatlist[filter] = []
            flatlist[filter].append(filename)
        elif (obstype == 'zero'):
            biaslist.append(filename)

    print "\nBIAS: ".join(biaslist)
    print "\nFLAT: ".join(flatlist)

    # Now create all biases
    print "computing bias frame"
    biasfile = "bias.fits"
    if (not os.path.isfile(biasfile)):
        podi_imcombine.imcombine(biaslist, biasfile, 
                                 operation='sigmaclipmean')

    
#    return
    print "computing flat-fields"
    # Now create all flatfields
    for filter in filterlist:
        print "\nWorking on filter %s\n" % (filter)

        flat_filename = "flat_%s.fits" % (filter)
        if (os.path.isfile(flat_filename)):
            continue

        normflat_list = []
        for filename in flatlist[filter]:
            print "Prepping flat %s" % (filename),
            biassub = podi_imarith.imarith(filename, "-", biasfile, None)

            # Now compute a median intensity across the field
            fullfield = numpy.array([])
            for ext in biassub[1:]:
                fullfield = numpy.append(fullfield, ext.data)
            flatlevel = numpy.median(fullfield)
            # flatlevel = 1.
            print " median-level=%.2f" % (flatlevel)

            # Normalize all Chips
            # for ext in range(1, len(biassub)):
            #     biassub[ext].data /= flatlevel
            for ext in biassub[1:]:
                ext.data /= flatlevel
                ext.data[ext.data < 0.05] = numpy.NaN

            # Write the normalized filename
            norm_flatfile = "tmp/%s.flatnorm.fits" % (filename[:-5])
            biassub.writeto(norm_flatfile, clobber=True)
            normflat_list.append(norm_flatfile)

        flathdu = podi_imcombine.imcombine(input_filelist=normflat_list, 
                                           outputfile=None, 
                                           operation='sigmaclipmean', 
                                           return_hdu=True)
        flathdu.writeto(flat_filename, clobber=True)



def reduce(filelist, label = "redu", 
           wcs_distortion=None):

    for filename in filelist:
        print "working on %s" % (filename)

        hdulist = pyfits.open(filename)
        
        obstype = hdulist[0].header['OBSTYPE']
        filter = minimo_filter(hdulist[0])
        objectname = hdulist[0].header['OBJECT'].replace(" ", "_") if 'OBJECT' in hdulist[0].header else "no_object"
        exptime = hdulist[0].header['EXPTIME'] if 'EXPTIME' in hdulist[0].header else 0.0

        if (not obstype == 'object'):
            continue

        outfile = "%s.%s.%s-%s-%.2f.fits" % (
            filename[:-5], label, objectname, filter, exptime)
        if (os.path.isfile(outfile)):
            print " ---> exists. skipping", outfile
            continue

        biasfile = "bias.fits"
        flat_filename = "flat_%s.fits" % (filter)

        bias_hdu = None
        if (os.path.isfile(biasfile)):
            bias_hdu = pyfits.open(biasfile)
        
        flat_hdu = None
        if (os.path.isfile(flat_filename)):
            flat_hdu = pyfits.open(flat_filename)


        hdulist[0].header['BIASFILE'] = biasfile
        hdulist[0].header['FLATFILE'] = flat_filename

        # Now handle the CD matrix to account for rotation
        # work out the input rotation:
        wcs = pyfits.open(wcs_distortion)                
        wcs_header = wcs['CCD1'].header
        in_rot = numpy.arctan2(hdulist['CCD1'].header['CD1_1'], hdulist['CCD1'].header['CD1_2'])
        out_rot = numpy.arctan2(wcs_header['CD1_1'], wcs_header['CD1_2'])
        rot_angle = out_rot - in_rot
        rot_matrix = numpy.array( [[+math.cos(rot_angle), -math.sin(rot_angle)] ,
                                   [+math.sin(rot_angle), +math.cos(rot_angle)]] )

        for ext_id in range(1, len(hdulist)):

            ext = hdulist[ext_id]

            # Search for bias extension
            if (not bias_hdu == None):
                for bias_ext in bias_hdu[1:]:
                    if (ext.name == bias_ext.name):
                        ext.data -= bias_ext.data
                        ext.header['BIASFILE'] = biasfile
                        break
                    
            if (not flat_hdu == None):
                for flat_ext in flat_hdu[1:]:
                    if (ext.name == flat_ext.name):
                        ext.data /= flat_ext.data
                        ext.header['FLATFILE'] = flat_filename
                        break
                
            # # Add the OTA keyword to make the file compatible with ccmatch 
            # # and the per-detector optimization
            ext.header['OTA'] = ext_id

            if (not wcs_distortion == None and os.path.isfile(wcs_distortion)):
                #print "Adding header from WCS minifits (%s)" % (extname)
                wcs = pyfits.open(wcs_distortion)                
                wcs_header = wcs[ext.name].header
                

                # # Modify the WCS solution to properly represents the WCS solution of binned data
                # # print "Correcting the WCS solution for binning factor",binning,"..."
                # for hdr_name in ('CRPIX1', 'CRPIX2'):
                #     wcs_header[hdr_name] /= binning
                #     for hdr_name in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
                #         wcs_header[hdr_name] *= binning

                # reduction_files_used['wcs'] = options['wcs_distortion']

                cards = wcs_header.cards
                for (keyword, value, comment) in cards:
                    if (keyword in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']):
                        continue
                    elif (keyword == 'CRVAL1'):
                        ext.header["CRVAL1"] -= wcs_header['CRVAL1'] / math.cos(math.radians(ext.header['CRVAL2']))
                    elif (keyword == "CRVAL2"):
                        ext.header['CRVAL2'] -= wcs_header['CRVAL2']
                    else:
                        ext.header[keyword] = (value, comment)

                cd_wcs = numpy.array([ [wcs_header['CD1_1'], wcs_header['CD1_2']],
                                       [wcs_header['CD2_1'], wcs_header['CD2_2']] ])
                cd_out = cd_wcs #.dot(rot_matrix)
                ext.header['CD1_1'] = cd_out[0,0]
                ext.header['CD1_2'] = cd_out[0,1]
                ext.header['CD2_1'] = cd_out[1,0]
                ext.header['CD2_2'] = cd_out[1,1]

                # Make sure to write RAs that are positive
                if (ext.header["CRVAL1"] < 0):
                    ext.header['CRVAL1'] += 360.
            
                wcs.close()

        # Now we have bias and flat applied
        hdulist.writeto(outfile, clobber=True)
        print "results written to", outfile
    
def multi_sex4scamp(filelist):
    
    for filename in filelist:

        catfile = filename[:-5]+".cat"

        if (os.path.isfile(catfile)):
            continue

        sex_config_file = "sex4scamp.conf"
        parameters_file = "sex4scamp.param"
        sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
            sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
            filename)

        start_time = time.time()
        print "starting SourceExtractor for %s --> %s" % (filename, catfile)
        try:
            ret = subprocess.Popen(sexcmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
            (sex_stdout, sex_stderr) = ret.communicate()
            #os.system(sexcmd)
            # if (ret.returncode != 0):
            #     logger.warning("Sextractor might have a problem, check the log")
            #     logger.debug("Stdout=\n"+sex_stdout)
            #     logger.debug("Stderr=\n"+sex_stderr)
        except OSError as e:
            # podi_logging.log_exception()
            print >>sys.stderr, "Execution failed:", e
        end_time = time.time()
        # logger.debug("SourceExtractor returned after %.3f seconds" % (end_time - start_time))
        print "SourceExtractor returned after %.3f seconds" % (end_time - start_time)
        pass



def minimo_fixwcs(infile, outfile):
    import dev_ccmatch

    hdulist = pyfits.open(infile)

    #
    # Run Sextractor
    #
    catfile = infile[:-5]+".wcscat"
    sex_config_file = "%s/.config/wcsfix.sex" % (sitesetup.exec_dir)
    parameters_file = "%s/.config/wcsfix.sexparam" % (sitesetup.exec_dir)
    sexcmd = "%s -c %s -PARAMETERS_NAME %s -CATALOG_NAME %s %s" % (
        sitesetup.sextractor, sex_config_file, parameters_file, catfile, 
        infile)

    if (not os.path.isfile(catfile)):
        print "Starting SourceExtractor"
        start_time = time.time()
        try:
            ret = subprocess.Popen(sexcmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
            (sex_stdout, sex_stderr) = ret.communicate()
            #os.system(sexcmd)
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e
        end_time = time.time()
        print "sextractor done after", end_time-start_time,"seconds"
    else:
        print "Catalog exists, reusing it", catfile

    # load catalog
    src_catalog = numpy.loadtxt(catfile)
    print src_catalog.shape

    #
    # run ccmatch
    #
    print "Running ccmatch"
    ccret = dev_ccmatch.ccmatch(source_catalog=src_catalog, 
                                reference_catalog=None, 
                                input_hdu=hdulist,
                                mode="rotation",
                                max_pointing_error=[25],
                                max_rotator_error=[-5,5],
                                angle_steps=30,
                                fov=0.3)

    hdu_out = ccret['hdulist']
    odi_2mass_matched = ccret['matched_src+2mass']

    # Create the WCS scatter plot
    try:
        import podi_diagnosticplots
        title_info = hdu_out[0].header.copy()

        plotfilename = "wcs1"
        podi_diagnosticplots.wcsdiag_scatter(matched_radec_odi=odi_2mass_matched[:,0:2], 
                                         matched_radec_2mass=odi_2mass_matched[:,-2:],
                                         matched_ota=odi_2mass_matched[:,SXcolumn['ota']],
                                         filename=plotfilename, 
                                         options=None,
                                         ota_wcs_stats=ota_wcs_stats,
                                         also_plot_singleOTAs=False,
                                         title_info=title_info)
    except:
        pass

        # # Create the WCS shift plot
        # plotfilename = create_qa_filename(outputfile, "wcs2", options)
        # podi_diagnosticplots.wcsdiag_shift(matched_radec_odi=odi_2mass_matched[:,0:2],
        #                                    matched_radec_2mass=odi_2mass_matched[:,-2:],
        #                                    matched_ota=odi_2mass_matched[:,SXcolumn['ota']],
        #                                    filename=plotfilename, #outputfile[:-5]+".wcs2", 
        #                                    options=options,
        #                                    ota_wcs_stats=ota_wcs_stats,
        #                                    ota_outlines=ota_outlines,
        #                                    also_plot_singleOTAs=options['otalevelplots'],
        #                                    title_info=title_info)

    # Write output file
    hdu_out.writeto(outfile, clobber=True)
    print "output written to",outfile


if __name__ == "__main__":

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)

    if (cmdline_arg_isset("-overscan")):
        for inputfile in get_clean_cmdline()[1:]:
            print inputfile
            overscan_subtract(inputfile, 
                              inputfile[:-5]+".merged.fits"
                          )
    elif (cmdline_arg_isset("-calibs")):
        make_calibrations(get_clean_cmdline()[1:])

    elif (cmdline_arg_isset("-reduce")):
        wcs_dist = cmdline_arg_set_or_default("-wcs", None)
        label = cmdline_arg_set_or_default("-label", 'redu')
        reduce(get_clean_cmdline()[1:], label=label, wcs_distortion=wcs_dist)

    elif (cmdline_arg_isset("-sex4scamp")):
        multi_sex4scamp(get_clean_cmdline()[1:])

    elif (cmdline_arg_isset("-fixwcs")):
        label = cmdline_arg_set_or_default("-label", 'wcsfix')
        for filename in get_clean_cmdline()[1:]:
            print
            print "Starting work on", filename
            print
            outfile = "%s.%s.fits" % (filename[:-5], label)
            minimo_fixwcs(filename, outfile)

    print "done!"

    podi_logging.shutdown_logging(options)
