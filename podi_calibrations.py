#!/usr/bin/env python

import pyfits
import os
import logging

from podi_commandline import *
import podi_focalplanelayout
import podi_sitesetup as sitesetup


def check_filename_directory(given, default_filename):
    """
    Some of the options support either a directory or a filename. This function
    checks if the input is a directory or a filename. In the first case, add
    the specified default filename and return the filename. In the latter case,
    simply return the filename.

    Parameters
    ----------
    given : string

        The value specified by the user, either a directory or filename

    default_filename : string

        In the case the user specified only a directory, append the name of this
        filename

    Returns
    -------
    filename

    Example
    -------
    This function is called e.g. during bias-subtraction. Using the -bias
    option, the user can specify the bias file to be used during the reduction.
    However, the default way of specifying the bias file is to simply give the
    directory with the calibration files, and collectcells adds the 'bias.fits'
    during execution.

    """

    if (given == None):
        return None

    if (not type(given) == list):
        _g = [given]
    else:
        _g = given

    for g in _g:
        if (os.path.isfile(g)):
            return g
        else:
            fn = "%s/%s" % (g, default_filename)
            if (os.path.isfile(fn)):
                return fn

    return ""


class ODICalibrations(object):

    def __init__(self, cmdline_options=None, hdulist=None):

        # set command-line options to a safe value
        if (cmdline_options is None):
            cmdline_options = set_default_options()
        self.options = cmdline_options

        self.filtername = None
        self.ota = None

        self.fpl = None
        self.detector_glow = "yes"
        self.binning = 1

        if (hdulist is not None):
            self.filtername = hdulist[0].header['FILTER']
            self.fpl = podi_focalplanelayout.FocalPlaneLayout(inp=hdulist)

    #
    # general class utility functions
    #
    def verify_fn(self, fn):
        if (os.path.isfile(fn)):
            return fn
        return None

    #
    # BIAS #######################
    #
    def bias(self):
        return self.verify_fn(
            check_filename_directory(
                self.options['bias_dir'],
                "bias_bin%d.fits" % (self.binning)
            )
        )
    def apply_bias(self):
        return (self.options['bias_dir'] is not None)

    #
    # DARK #######################
    #
    def apply_dark(self):
        return (self.options['dark_dir'] is not None)
    def dark(self):
        return self.verify_fn(
            check_filename_directory(
                self.options['dark_dir'],
                "dark_%s_bin%d.fits" % (self.detector_glow, self.binning)
            )
        )


    #
    # FLAT #######################
    #
    def apply_flat(self):
        return (self.options['flat_dir'] is not None)
    def flat(self, flat_order=None):
        if (flat_order is None):
            flat_order = ['dflat']
        for ft in flat_order:
            flatfield_filename = check_filename_directory(
                self.options['flat_dir'],
                "%s_%s_bin%d.fits" % (ft, self.filtername, self.binning)
            )
            if (os.path.isfile(flatfield_filename)):
                return self.verify_fn(flatfield_filename)
        return None


    #
    # FRINGE #######################
    #
    def apply_fringe(self):
        return False

    def fringe(self, mjd=0, filtername=None, binning=None, detector=None):

        if (filtername is None):
            filtername = self.filtername
        if (binning is None):
            binning = self.binning
        if (detector is None):
            detector = fpl.layout.lower()
        history_fn = "%s/mastercals/fringe/%s/%s/fringe_bin%d.history" % (
            sitesetup.exec_dir, detector, filtername, binning
        )
        return None

    #
    # PUPILGHOST #######################
    #
    def apply_pupilghost(self):
        return False
    def pupilghost(self):
        return None


class CalibrationHistory(object):

    def __init__(self, filename):
        self.filename = filename

        self.mjd = []
        self.filenames = []
        self.url = []

        self.read()

        pass

    def read(self):
        with open(self.filename) as f:
            lines = f.readlines()

            for line in lines:
                if (line.startswith("#") or len(line) <= 0):
                    continue
                items = line.split()
                if (len(items) < 3):
                    continue

                self.mjd.append(float(items[0]))
                self.filenames.append(items[1])
                self.url.append(items[2])

    def find(self, mjd):
        pass

    def files_urls(self):
        return zip(self.filenames, self.url)

def download(url, local_fn):
    wget_cmd = "wget -O %s -nv %s" % (local_fn, url)
    print wget_cmd
    os.system(wget_cmd)

def update_local_mastercals(local_cache_dir=None):

    if (local_cache_dir is None):
        local_cache_dir = sitesetup.mastercal_cache
    print("Saving all WIYN-delivered master-cals to %s" % (local_cache_dir))
    if (not os.path.isdir(local_cache_dir)):
        os.mkdir(local_cache_dir)

    product_order = [
        ('fringe', ['odi_z', 'odi_i']),
        ('pupilghost', None),
    ]
    detector_layouts = ['podi', 'odi_5x6']
    binning = [1,2]

    mastercal_base = "%s/mastercals" % (sitesetup.exec_dir)
    for product, filterlist in product_order:
        for layout in detector_layouts:
            if (filterlist is None):
                dirlist = ["%s/%s/%s" % (mastercal_base, product, layout)]
            else:
                dirlist = ["%s/%s/%s/%s" % (mastercal_base, product, layout, filtername) for filtername in filterlist]

            for dirname in dirlist:

                for bin in binning:
                    print("Checking %s/%s in %s" % (product, layout, dirname))

                    history_fn = "%s/%s_bin%d.history" % (dirname, product, bin)
                    if (not os.path.isfile(history_fn)):
                        print("    Unable to find calibration history file (%s)" % (history_fn))
                        continue
                    else:
                        print("  Reading calibration history from %s" % (history_fn))

                    history = CalibrationHistory(history_fn)
                    for fn, url in history.files_urls():
                        local_fn = "%s/%s" % (local_cache_dir, fn)
                        if (not os.path.isfile(local_fn)):
                            download(url, local_fn)
                        else:
                            print("    --> file (%s) already available locally" % (fn))

if __name__ == "__main__":

    # use this as a test-case and stand-alone tool to download calibrations from the web
    if (sys.argv[1] == "update"):
        update_local_mastercals()
    pass