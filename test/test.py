#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy

from podi_definitions import *

FWHM = 2 *math.sqrt(2*math.log(2))

    
if __name__ == "__main__":

    

    inputfile = sys.argv[1]
    #findstar_allota(inputfile)

    import cProfile
    cProfile.run('findstar_allota(inputfile)')
    
    
        
