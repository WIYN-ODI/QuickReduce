#!/usr/bin/env python3

import os
import sys
import numpy
import pyfits

try:
    import cPickle as pickle
except:
    import pickle

#
# Airmass terms from Table 5 in Landolt 2007,
# http://adsabs.harvard.edu/abs/2007ASPC..364...27L
#

known_filters = {
    "Us_solid": (24.0, 0.586, "purple"),
    "odi_g":    (26.20, 0.264, "blue"),
    "odi_r":    (26.13, 0.122, "green"),
    "odi_i":    (26.0, 0.077, "orange"),
    "odi_z":    (24.8, 0.050, "red"),
}


def seconds2mjd(sec):
    return sec/86400.


def read_data_from_files(filelist):

    
    obstype, exptime, filtername, photzp, photzpe, mjd, dateobs, airmass = [], [], [], [], [], [], [], []

    
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


    print filelist
    for filename in filelist:

        if (not os.path.isfile(filename)):
            continue

        if (filename in direntry):
            # We know about this file from the pickle
            file_dir = direntry[filename]
            #(_obstype, _exptime, expmea_filtername, _photzp, _photzpe, _mjdobs, _dateobs) = direntry[filename]

        else:
            print "Getting info for file", filename
            try:
                hdulist = pyfits.open(filename)
                hdr = hdulist[0].header
            except IOError:
                pass
                continue
            except:
                raise

            #print "reading headers"

            file_dir = {
                "OBSTYPE" : "???", 
                "EXPTIME" : -1, 
                "EXPMEAS" : -1,
                "FILTER"  : "???", 
                "PHOTZP"  : -99, 
                "PHOTZPER" : -99, 
                "MJD-OBS" : -99, 
                "DATE-OBS": "???",
                "AIRMASS" : 1.0,
            }

            for header in file_dir:
                try:
                    file_dir[header] = hdr[header]
                except KeyError:
                    pass
                except:
                    raise

            direntry[filename] = file_dir

        #print file_dir #['OBSTYPE']
        obstype.append(file_dir['OBSTYPE'])
        exptime.append(file_dir['EXPTIME'])
        filtername.append(file_dir['FILTER'])
        photzp.append(file_dir['PHOTZP'])
        photzpe.append(file_dir['PHOTZPER'])
        mjd.append(file_dir['MJD-OBS'])
        dateobs.append(file_dir['DATE-OBS'])
        airmass.append(file_dir['AIRMASS'])

    # Now pickle the dir_entry for later use
    picklejar = "index.pickle"
    print "Pickling data..."
    with open(picklejar, "wb") as pf:
        pickle.dump(direntry, pf)
        pf.close()

    #print mjd

    return direntry, (obstype, exptime, filtername, photzp, photzpe, mjd, dateobs, airmass)
