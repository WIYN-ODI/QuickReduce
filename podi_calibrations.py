#!/usr/bin/env python3

import astropy.io.fits as pyfits
import os
import logging

from podi_commandline import *
import podi_focalplanelayout
from podi_definitions import get_binning, get_filter_level

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

    if (given is None):
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
        self.logger = logging.getLogger("Calibrations")

        if (hdulist is not None):
            self.filtername = hdulist[0].header['FILTER']
            self.fpl = podi_focalplanelayout.FocalPlaneLayout(inp=hdulist)
            self.binning = get_binning(hdulist[0].header)
            self.filter_level = get_filter_level(hdulist[0].header)

        self.mastercal_dir = "%s/mastercals" % (sitesetup.exec_dir)

    #
    # general class utility functions
    #
    def verify_fn(self, fn):
        self.logger.debug("Checking if file exists: %s" % (fn))
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
                self.logger.debug("Could not find fringe history file (%s)" % (history_fn))
                return None

            self.logger.debug("Checking fringe history file in %s ..." % (history_fn))
            hist = CalibrationHistory(history_fn, self.logger)
            fn = hist.find(mjd)
            if (fn is None):
                return None

            # check if the file is installed locally
            check_dirs = [fringe_dir, sitesetup.mastercal_cache]
            for cd in check_dirs:
                full_filename = os.path.join(cd, fn)
                if (os.path.isfile(full_filename)):
                    self.logger.debug("Found fringe frame: %s" % (full_filename))
                    break

                # full_filename = "%s/%s" % (sitesetup.mastercal_cache, fn)
                # if (os.path.isfile(full_filename)):
                #     self.logger.debug("Found fringe frame: %s" % (full_filename))
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
        # print "considering PG correction",(not self.options['pupilghost_dir'] == False)
        return (not self.options['pupilghost_dir'] == False)

    def pupilghost(self, mjd=0, filter_level=None, binning=None):

        if (filter_level is None):
            filter_level = self.filter_level
        if (binning is None):
            binning = self.binning

        # print "Looking for pupilghost filename"

        pg_dir = self.options['pupilghost_dir']
        if (pg_dir is not None and os.path.isfile(pg_dir)):
            full_filename = pg_dir
        else:
            if (pg_dir is None or not os.path.isdir(pg_dir)):
                pg_dir = "%s/pupilghost/%s/level%d/" % (self.mastercal_dir, self.fpl.layout.lower(), filter_level)

            history_fn = "%s/pupilghost_bin%d.history" % (pg_dir, binning)
            hist = CalibrationHistory(history_fn)
            fn = hist.find(mjd)

            if (fn is None):
                return None

            full_filename = "%s/%s" % (sitesetup.mastercal_cache, fn)

        return self.verify_fn(full_filename)


    #
    # BAD PIXEL MASK #######################
    #
    def apply_bpm(self):
        # print self.options['bpm_dir']
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
            # print mjd, full_filename

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
        if (wcs_param == "plain"):
            return "plain"
        elif (wcs_param is not None and os.path.isfile(wcs_param)):
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


    #
    # CROSSTALK #######################
    #
    def apply_crosstalk(self):
        # print "XTALK:", self.options['crosstalk']
        return (not self.options['crosstalk'] == False)

    def crosstalk(self, mjd=0, ota=None):

        xtalk = self.options['crosstalk']
        if (xtalk is not None and os.path.isfile(xtalk)):
            full_filename = xtalk
        else:

            if (xtalk is not None and os.path.isdir(xtalk)):
                xtalk_dir = xtalk
            else:
                xtalk_dir = "%s/crosstalk/%s" % (self.mastercal_dir, self.fpl.layout.lower())

            # print "XTALK-DIR: ", xtalk_dir
            history_fn = "%s/crosstalk.history" % (xtalk_dir)
            hist = CalibrationHistory(history_fn)
            fn = hist.find(mjd)

            if (fn is None):
                return None

            full_filename = "%s/%s" % (xtalk_dir, fn)
            # print full_filename

        return self.verify_fn(full_filename)


    #
    # PHOTOMETRIC FLAT #######################
    #
    def apply_photflat(self):
        return (self.options['photflat'] is not None)

    def photflat(self):
        photflat_filename = check_filename_directory(
                self.options['photflat'],
                "photflat_%s_bin%d.fits" % (self.filtername, self.binning)
        )
        if (os.path.isfile(photflat_filename)):
            return self.verify_fn(photflat_filename)
        return None




    def _info(self, step, apply_fct, filename_fct):

        try:
            apply = apply_fct()
        except:
            apply = "ERROR"
            pass

        try:
            filename = filename_fct()
        except:
            filename = "ERROR"
            pass

        if (filename is None):
            filename = "--undefined--"
        self.logger.info("%s: %5s (%s)" % (step, apply, filename))

    def info(self):
        self._info("BIAS", self.apply_bias, self.bias)
        self._info("DARK", self.apply_dark, self.dark)
        self._info("FLAT", self.apply_flat, self.flat)
        self._info("FRINGE", self.apply_fringe, self.fringe)
        self._info("PUPILGHOST", self.apply_pupilghost, self.pupilghost)
        self._info("BPM", self.apply_bpm(), self.bpm)
        self._info("NONLINEARITY", self.apply_nonlinearity, self.nonlinearity)
        self._info("WCS", self.apply_wcs, self.wcs)
        #self._info("CROSSTALK", self.apply_crosstalk, self.crosstalk)

class CalibrationHistory(object):

    def __init__(self, filename, logger=None):
        self.filename = filename

        self.mjd = []
        self.filenames = []
        self.url = []

        if (logger is None):
            logger = logging.getLogger("CalibHistory")
        self.logger = logger

        self.read()

        pass

    def read(self):
        if (not os.path.isfile(self.filename)):
            self.logger.debug("Unable to read calibration history in %s" % (self.filename))
            return

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
        self.logger.debug("Read a total of %d entries from %s" % (len(self.filenames), self.filename))

    def find(self, mjd):
        self.logger.debug("Searching for calibration product for MJD = %.5f" % (mjd))
        fn = None
        for i, i_mjd in enumerate(self.mjd):
            # print "FIND: ", mjd, i_mjd, self.filenames[i]
            if (mjd >= i_mjd):
                fn = self.filenames[i]
            else:
                break
        self.logger.debug("Found: %s" % (fn))
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

    progress = ["--quiet --progress=bar --show-progress",
                "--quiet --progress=bar",
                ""]
    for p in progress:
        wget_cmd = "wget %s -O %s %s" % (p, local_fn, url)
        logger.debug(wget_cmd)
        ret = os.system(wget_cmd)
        if (ret == 0):
            break

    if (ret != 0):
        logger.error("Unable to run wget")


def update_local_mastercals(local_cache_dir=None):

    if (local_cache_dir is None):
        local_cache_dir = sitesetup.mastercal_cache
    print("Saving all WIYN-delivered master-cals to %s" % (local_cache_dir))
    if (not os.path.isdir(local_cache_dir)):
        os.mkdir(local_cache_dir)

    product_order = [
        ('fringe', ['odi_z', 'odi_i']),
        ('pupilghost', ['level1', 'level2', 'level3']),
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
                        # print("    Unable to find calibration history file (%s)" % (history_fn))
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


