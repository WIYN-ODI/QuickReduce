

from collections import OrderedDict
from podi_definitions import add_fits_header_title

import logging

class ReductionLog:

    _attempt = 'attempted'
    _fail = 'failed'
    _success = 'successful'
    _not_selected = 'not selected'
    _skip_after_fail = 'skipped after fail'
    _partial_fail = 'partially failed'
    _missing_data = 'missing data'
    _partial_success = 'partial success'
    _not_required = 'not required'
    _no_data = 'no data'

    def __init__(self):

        # self._attempt = 'attempted'
        # self._fail = 'failed'
        # self._success = 'successful'
        # self._not_selected = 'not selected'
        # self._skip_after_fail = 'skipped after fail'

        self.steps = OrderedDict()
        self.steps['overscan']    = (None, '_OVERSCN', 'overscan subtraction')
        self.steps['crosstalk']   = (None, '_XTALK'  , 'crosstalk correction')
        self.steps['gain']        = (None, '_GAIN'   , 'explicit gain correction')
        self.steps['trimcell']    = (None, '_TRIMCEL', 'trim edges of all cells')
        self.steps['nonlinearity']= (None, '_NONLIN' , 'non-linearity correction')
        self.steps['bias']        = (None, '_BIAS'   , 'bias subtraction')
        self.steps['dark']        = (None, '_DARK'   , 'dark correction')
        self.steps['flat']        = (None, '_FLAT'   , 'flatfielding')
        self.steps['illumcorr']   = (None, '_ILLUMCR', 'illumination correction')
        self.steps['badpixels']   = (None, '_BADPXLS', 'bad pixel masks applied')
        self.steps['exptime_norm']= (None, '_EXTNORM', 'normalize to counts per sec')
        self.steps['persistency'] = (None, '_PERSIST', 'persistency masking')
        self.steps['saturation']  = (None, '_SATURAT', 'saturation trail masking') 
        self.steps['wcs_dist']    = (None, '_WCSDIST', 'import WCS distorion')
        self.steps['crj']         = (None, '_COSMICS', 'cosmic ray rejection')

        self.steps['pupilghost']  = (None, '_PGHOST' , 'pupilghost removal')
        self.steps['fringe']      = (None, '_FRINGE' , 'fringe removal')
        self.steps['wcscal']      = (None, '_WCSCAL' , 'astrometric calibration')
        self.steps['nonsidereal'] = (None, '_NONSDRL', 'non-sidereal correction')
        self.steps['photcal']     = (None, '_PHOTCAL', 'photometric ZP calibration')
        self.steps['softwarebin'] = (None, '_SOFTBIN', 'software binning')

        
    def set(self, task, outcome):
        if (task in self.steps):
            _, hdrkey, comment = self.steps[task]
            self.steps[task] = (outcome, hdrkey, comment)
            logging.getLogger("ReductionLog").debug("Setting reduction outcome: %s --> %s" % (task, outcome))
            return True
        logging.getLogger("ReductionLog").warning("Trying to set unknown reduction parameter: %s" % (task))
        return False

    def attempt(self, task):
        return self.set(task, self._attempt)
    def success(self, task):
        return self.set(task, self._success)
    def fail(self, task):
        return self.set(task, self._fail)
    def not_selected(self, task):
        return self.set(task, self._not_selected)
    def skip_after_fail(self, task):
        return self.set(task, self._skip_after_fail)
    def partial_fail(self, task):
        return self.set(task, self._partial_fail)
    def missing_data(self, task):
        return self.set(task, self._missing_data)
    def not_required(self, task):
        return self.set(task, self._not_required)
    def no_data(self, task):
        return self.set(task, self._no_data)

    def dump(self):
        for step in self.steps:
            try:
                outcome, hdrkey, comment = self.steps[step]
                logging.getLogger("ReductionLog").info("%8s = %s / %s" % (hdrkey, outcome, comment))
            except:
                pass

    def write_to_header(self, hdr):
        logging.getLogger("ReductionLog").debug("Adding reduction log to FITS header")
        firstkey = None
        for step in self.steps:
            try:
                outcome, hdrkey, comment = self.steps[step]
                firstkey = hdrkey if firstkey == None else firstkey
                hdr.append((hdrkey, outcome, comment))
            except:
                pass
        # Add some fits header title
        add_fits_header_title(hdr, "Reduction log", firstkey)
        logging.getLogger("ReductionLog").debug("All reduction steps added")
        return

    def get(self, step):
        if (step in self.steps):
            outcome, hdrkey, comment = self.steps[step]
            return outcome
        return None


    def combine(self, redlog):

        logging.getLogger("ReductionLog").debug("Combining two reduction logs")
        for step in self.steps:
            other = redlog.get(step)
            this = self.get(step)

            if (this == other):
                # They are the same, so don't change anything
                pass

            elif (this == None):
                # If the current state is undefined, use the other one
                self.set(step, other)

            elif (this in [self._success,
                           self._partial_fail] and 
                  other in [self._fail,
                            self._partial_fail,
                            self._skip_after_fail]):
                # Set to partial fail if mix of success and fails
                self.set(step, self._partial_fail)

            elif (this in [self._not_selected, 
                           self._partial_success] and 
                  other in [self._success]):

                # if some are done and some or not, set to partial 
                # (NOT the same as partial fail)
                self.set(step, self._partial_success)
                
            elif (this in [self._not_required,
                           self._success] and
                  other in [self._success,
                            self._not_required]):

                self.set(step, self._success)
                
            else:
                pass
                logging.getLogger("ReductionLog").warning("Unable to handle outcome combination: %s & %s" % (
                    this, other))

        logging.getLogger("ReductionLog").debug("done combining reduction logs")


