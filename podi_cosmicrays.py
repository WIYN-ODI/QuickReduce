#!/usr/bin/env python

"""

Wrapper around the LACosmics python implementation by Malte Tewes. Also provides
a stand-alone funtionality, based on Malte's demo program.

"""

import ext_cosmics as cosmics
import pyfits
import numpy
from podi_definitions import *
import itertools



def remove_cosmics(hdu):

    gain = hdu.header['GAIN']
    readnoise = 8 #hdu.header['RON']

    corrected = hdu.data.copy()
    mask = numpy.zeros_like(hdu.data)

    # Break up each entire OTA into its cells to see if this is faster
    for cx, cy in itertools.product(range(8),repeat=2):
        print "working on cell",cx, cy

        x1, x2,y1, y2 = cell2ota__get_target_region(cx, cy, binning=1)
        cell_data = hdu.data[x1:x2, y1:y2]

        c = cosmics.cosmicsimage(cell_data, 
                                 gain=gain, readnoise=readnoise, 
                                 sigclip = 5.0, sigfrac = 0.3, objlim = 5.0,
                                 verbose=False)
        c.run(maxiter=4)

        # # Re-insert the data into the full frame
        corrected[x1:x2, y1:y2] = c.cleanarray
        mask[x1:x2, y1:y2][c.mask] = 1

    return corrected, mask

    # # Prepare all structures
    # c = ext_lacosmics.cosmicsimage(hdu.data, 
    #                                gain=gain, readnoise=readnoise, 
    #                                sigclip = 5.0, sigfrac = 0.3, objlim = 5.0)

    # # Run the cosmic removal
    # c.run(maxiter=4)

    # return c.cleanarray, c.mask



if __name__ == "__main__":

    if (cmdline_arg_isset("-standalone")):
        inputfile = get_clean_cmdline()[1]
        hdulist = pyfits.open(inputfile)

        import yappi
        yappi.start()

        for i in range(2): #len(hdulist)):
            if (not is_image_extension(hdulist[i])):
                continue

            if ('OTA' in hdulist[i].header):
                print "Removing cosmics from OTA",hdulist[i].header['OTA']
            else:
                print "removing cosmics"

            fixed, mask = remove_cosmics(hdulist[i])
            hdulist[i].data = fixed
            
            ota = hdulist[i].header['OTA']
            pyfits.HDUList([pyfits.PrimaryHDU(data=fixed, header=hdulist[i].header)]).writeto("cleaned_OTA%d.fits" % (ota), clobber=True)
            pyfits.HDUList([pyfits.PrimaryHDU(data=mask, header=hdulist[i].header)]).writeto("mask_OTA%d.fits" % (ota), clobber=True)
            # ignore the mask for now

        outputfile = get_clean_cmdline()[2]
        hdulist.writeto(outputfile, clobber=True)

        yappi.get_thread_stats().print_all()
        yappi.get_func_stats().debug_print() #print_all()

    else:
        print __doc__
        print cosmics.__doc__
