#!/usr/bin/env python

import pyfits
import os
import logging

from podi_commandline import *
import podi_focalplanelayout

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

        self.mastercal_dir = "%s/mastercals" % (sitesetup.exec_dir)

    #
    # general class utility functions
    #
    def verify_fn(self, fn):
        if (fn is None):
            return None
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
        x = not self.options['fringe_dir'] == False
        # print "Apply fringe?", x
        return (not self.options['fringe_dir'] == False)


    def fringe(self, mjd=0, filtername=None, binning=None, detector=None):

        if (filtername is None):
            filtername = self.filtername
        if (binning is None):
            binning = self.binning
        if (detector is None):
            detector = self.fpl.layout.lower()

        fringe_param = self.options['fringe_dir']
        if (fringe_param is not None and os.path.isfile(fringe_param)):
            full_filename = fringe_param
        else:
            if (fringe_param is not None and os.path.isdir(fringe_param)):
                fringe_dir = fringe_param
            else:
                fringe_dir = "%s/fringe/%s/%s" % (self.mastercal_dir, detector, filtername)

            history_fn = "%s/fringe_bin%d.history" % (fringe_dir, binning)
            if (not os.path.isfile(history_fn)):
                return None
            hist = CalibrationHistory(history_fn)
            fn = hist.find(mjd)
            if (fn is None):
                return None

            full_filename = "%s/%s" % (sitesetup.mastercal_cache, fn)

        # print "Fringe template", full_filename
        return self.verify_fn(full_filename)


    def fringevector(self, ota, filtername=None, binning=None):

        if (filtername is None):
            filtername = self.filtername
        if (binning is None):
            binning = self.binning

        # print "fringe vector dir:", self.options['fringe_vectors']
        fringe_vector_dir = self.options['fringe_vectors']
        if (fringe_vector_dir is None or not os.path.isdir(fringe_vector_dir)):
            fringe_vector_dir = "%s/fringevectors/%s/%s" % (self.mastercal_dir, self.fpl.layout.lower(), filtername)

        # print fringe_vector_dir
        if (not os.path.isdir(fringe_vector_dir)):
            return None

        full_filename = "%s/fringevectors__%s__OTA%02d.reg" % (fringe_vector_dir, self.filtername, ota)
        # print "fringe vector:", full_filename
        return self.verify_fn(full_filename)


    #
    # PUPILGHOST #######################
    #
    def apply_pupilghost(self):
        return False
    def pupilghost(self):
        return None


    #
    # BAD PIXEL MASK #######################
    #
    def apply_bpm(self):
        #print self.options['bpm_dir']
        return (not self.options['bpm_dir'] == False)
    def bpm(self, ota, fpl=None, layout=None):

        if (fpl is None):
            fpl = self.fpl
        if (layout is None):
            layout = fpl.layout

        fn = "bpm_%s.reg" % (ota)

        bpm_dir = self.options['bpm_dir']
        if (bpm_dir is None):
            return None

        elif (os.path.isfile(bpm_dir)):
            return bpm_dir

        elif (os.path.isdir(bpm_dir)):
            # print("BPM %s is directory" % (bpm_dir))
            ffn = "%s/%s" % (bpm_dir, fn)
            return ffn if os.path.isfile(ffn) else None

        else:
            ffn = "%s/bpm/%s/%s" % (
                self.mastercal_dir, layout.lower(), fn)
            # print(ffn)
            return ffn if os.path.isfile(ffn) else None

        return None


    #
    # NON-LINEARITY CORRECTION #######################
    #
    def apply_nonlinearity(self):
        return (not self.options['nonlinearity'] == False) #is not None)

    def nonlinearity(self, mjd):

        """

        Select the appropriate non-linearity coefficient file based on a history
        file and a given MJD timestamp.

        """
        if (self.options['nonlinearity'] is not None and
                os.path.isfile(self.options['nonlinearity'])):
            full_filename = self.options['nonlinearity']
        else:
            if (self.options['nonlinearity'] is not None and
                    os.path.isdir(self.options['nonlinearity'])):
                nl_basedir = self.options['nonlinearity']
            else:
                # Construct the name of the history file
                nl_basedir = "%s/nonlinearity/" % (self.mastercal_dir)
            history_filename = "%s/nonlinearity.history" % (nl_basedir)

            # Read the file and determine which coefficient file is the best one.
            history_file = open(history_filename, "r")
            hist = CalibrationHistory(history_filename)
            fn = hist.find(mjd)
            # hist.info()

            if (fn is None):
                return None

            # Now assemble the entire filename
            full_filename = "%s/%s" % (nl_basedir, fn)

        return self.verify_fn(full_filename)


    def apply_relative_gain(self):
        return (self.options['gain_method'] == "relative")
    def relative_gain(self, mjd):
        return self.nonlinearity(mjd)



    #
    # WCS & DISTORTION MODEL #######################
    #
    def apply_wcs(self):
        return (not self.options['wcs_distortion'] == False)
    def wcs(self, mjd=0):
        wcs_param = self.options['wcs_distortion']
        if (wcs_param is not None and os.path.isfile(wcs_param)):
            full_fn = wcs_param
        else:
            if (wcs_param is not None and os.path.isdir(wcs_param)):
                history_dir = wcs_param
            else:
                history_dir = "%s/wcs/%s" % (self.mastercal_dir, self.fpl.layout.lower())

            # read history file and pick a suitable file
            history_fn = "%s/wcs.history" % (history_dir)
            # print "Reading WCS history file", history_fn
            hist = CalibrationHistory(history_fn)
            # print hist.mjd, hist.filenames
            fn = hist.find(mjd)
            # print "FOUND WCS", fn

            full_fn = "%s/%s" % (history_dir, fn)

        # print full_fn
        return self.verify_fn(full_fn)





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
            # print lines
            for line in lines:
                if (line.startswith("#") or len(line) <= 0):
                    continue
                items = line.split()
                if (len(items) < 2):
                    continue

                self.mjd.append(float(items[0]))
                self.filenames.append(items[1])

                _url = items[2] if len(items) > 2 else None
                self.url.append(_url)

    def find(self, mjd):

        fn = None
        for i, i_mjd in enumerate(self.mjd):
            # print "FIND: ", mjd, i_mjd, self.filenames[i]
            if (mjd >= i_mjd):
                fn = self.filenames[i]
                break
        if (fn is not None):
            return fn
        return None

    def files_urls(self):
        return zip(self.filenames, self.url)


    def info(self):
        print("\nHISTORY - %s" % (self.filename))
        for i in range(len(self.mjd)):
            print("%9.3f: %s (%s)" % (self.mjd[i], self.filenames[i], self.url[i]) )
        print("-"*10+"\n")



def download(url, local_fn):

    logger = logging.getLogger("Download")
    wget_cmd = "wget --quiet --progress=bar --show-progress -O %s %s" % (local_fn, url)
    logger.debug(wget_cmd)
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


