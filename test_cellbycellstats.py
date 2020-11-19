#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import podi_diagnosticplots

if __name__ == "__main__":

    fn = sys.argv[1]
    hdulist = pyfits.open(fn)

    
    podi_diagnosticplots.plot_cellbycell_stats(
        hdulist=hdulist,
        title="test plot",
        vmin=0.95, vmax=1.05,
        plotfile="difftest.pdf",
    )