#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import sys

"""
crpix1[cell=0,0] (otaX) = -4222.7484375 * otaX + 15973.0199279
crpix2[cell=0,0] (otaY) = -4315.55125 * otaY + 12702.5552885

crpix1[x,y] = crpix[0,0] + -508.0 *cellX
crpix2[y,y] = crpix[0,0] -505.0 * cellY
"""

ltv = """OTA 61, EXTNAME xy00:   LTV1 = 0        LTV2 = -3535
OTA 61, EXTNAME xy10:   LTV1 = -508     LTV2 = -3535
OTA 61, EXTNAME xy20:   LTV1 = -1016    LTV2 = -3535
OTA 61, EXTNAME xy30:   LTV1 = -1524    LTV2 = -3535
OTA 61, EXTNAME xy40:   LTV1 = -2032    LTV2 = -3535
OTA 61, EXTNAME xy50:   LTV1 = -2540    LTV2 = -3535
OTA 61, EXTNAME xy60:   LTV1 = -3048    LTV2 = -3535
OTA 61, EXTNAME xy70:   LTV1 = -3556    LTV2 = -3535
OTA 61, EXTNAME xy01:   LTV1 = 0        LTV2 = -3030
OTA 61, EXTNAME xy11:   LTV1 = -508     LTV2 = -3030
OTA 61, EXTNAME xy21:   LTV1 = -1016    LTV2 = -3030
OTA 61, EXTNAME xy31:   LTV1 = -1524    LTV2 = -3030
OTA 61, EXTNAME xy41:   LTV1 = -2032    LTV2 = -3030
OTA 61, EXTNAME xy51:   LTV1 = -2540    LTV2 = -3030
OTA 61, EXTNAME xy61:   LTV1 = -3048    LTV2 = -3030
OTA 61, EXTNAME xy71:   LTV1 = -3556    LTV2 = -3030
OTA 61, EXTNAME xy02:   LTV1 = 0        LTV2 = -2525
OTA 61, EXTNAME xy12:   LTV1 = -508     LTV2 = -2525
OTA 61, EXTNAME xy22:   LTV1 = -1016    LTV2 = -2525
OTA 61, EXTNAME xy32:   LTV1 = -1524    LTV2 = -2525
OTA 61, EXTNAME xy42:   LTV1 = -2032    LTV2 = -2525
OTA 61, EXTNAME xy52:   LTV1 = -2540    LTV2 = -2525
OTA 61, EXTNAME xy62:   LTV1 = -3048    LTV2 = -2525
OTA 61, EXTNAME xy72:   LTV1 = -3556    LTV2 = -2525
OTA 61, EXTNAME xy03:   LTV1 = 0        LTV2 = -2020
OTA 61, EXTNAME xy13:   LTV1 = -508     LTV2 = -2020
OTA 61, EXTNAME xy23:   LTV1 = -1016    LTV2 = -2020
OTA 61, EXTNAME xy33:   LTV1 = -1524    LTV2 = -2020
OTA 61, EXTNAME xy43:   LTV1 = -2032    LTV2 = -2020
OTA 61, EXTNAME xy53:   LTV1 = -2540    LTV2 = -2020
OTA 61, EXTNAME xy63:   LTV1 = -3048    LTV2 = -2020
OTA 61, EXTNAME xy73:   LTV1 = -3556    LTV2 = -2020
OTA 61, EXTNAME xy04:   LTV1 = 0        LTV2 = -1515
OTA 61, EXTNAME xy14:   LTV1 = -508     LTV2 = -1515
OTA 61, EXTNAME xy24:   LTV1 = -1016    LTV2 = -1515
OTA 61, EXTNAME xy34:   LTV1 = -1524    LTV2 = -1515
OTA 61, EXTNAME xy44:   LTV1 = -2032    LTV2 = -1515
OTA 61, EXTNAME xy54:   LTV1 = -2540    LTV2 = -1515
OTA 61, EXTNAME xy64:   LTV1 = -3048    LTV2 = -1515
OTA 61, EXTNAME xy74:   LTV1 = -3556    LTV2 = -1515
OTA 61, EXTNAME xy05:   LTV1 = 0        LTV2 = -1010
OTA 61, EXTNAME xy15:   LTV1 = -508     LTV2 = -1010
OTA 61, EXTNAME xy25:   LTV1 = -1016    LTV2 = -1010
OTA 61, EXTNAME xy35:   LTV1 = -1524    LTV2 = -1010
OTA 61, EXTNAME xy45:   LTV1 = -2032    LTV2 = -1010
OTA 61, EXTNAME xy55:   LTV1 = -2540    LTV2 = -1010
OTA 61, EXTNAME xy65:   LTV1 = -3048    LTV2 = -1010
OTA 61, EXTNAME xy75:   LTV1 = -3556    LTV2 = -1010
OTA 61, EXTNAME xy06:   LTV1 = 0        LTV2 = -505
OTA 61, EXTNAME xy16:   LTV1 = -508     LTV2 = -505
OTA 61, EXTNAME xy26:   LTV1 = -1016    LTV2 = -505
OTA 61, EXTNAME xy36:   LTV1 = -1524    LTV2 = -505
OTA 61, EXTNAME xy46:   LTV1 = -2032    LTV2 = -505
OTA 61, EXTNAME xy56:   LTV1 = -2540    LTV2 = -505
OTA 61, EXTNAME xy66:   LTV1 = -3048    LTV2 = -505
OTA 61, EXTNAME xy76:   LTV1 = -3556    LTV2 = -505
OTA 61, EXTNAME xy07:   LTV1 = 0        LTV2 = 0
OTA 61, EXTNAME xy17:   LTV1 = -508     LTV2 = 0
OTA 61, EXTNAME xy27:   LTV1 = -1016    LTV2 = 0
OTA 61, EXTNAME xy37:   LTV1 = -1524    LTV2 = 0
OTA 61, EXTNAME xy47:   LTV1 = -2032    LTV2 = 0
OTA 61, EXTNAME xy57:   LTV1 = -2540    LTV2 = 0
OTA 61, EXTNAME xy67:   LTV1 = -3048    LTV2 = 0
OTA 61, EXTNAME xy77:   LTV1 = -3556    LTV2 = 0"""

def plot_crval () :
    
    
    data = np.loadtxt("/tmp/test", unpack="True")
       
    otaX = data[0] // 10 
    otaY = data[0] % 10
    
    cellX = data[1] // 10
    cellY = data[1] % 10
    
    fig = plt.figure()
    
    y = data[3][(otaX == 1) ]
    x = cellY[(otaX == 1) ]
    
    
    plt.plot(x, y, ".")
    m, b = np.polyfit(x, y, 1) 
    print m, b
    plt.plot (x, m * x + b, '--k')
    plt.savefig ("crval.png")
    return


def printWCS ():

    otalist = [pyfits.PrimaryHDU()]

    cx, cy = 5,5
    for yy in range(0, 8):
        for xx in range(0, 8):
            imghdu = pyfits.ImageHDU()
            
            imghdu.header['WCSASTRM'] = 'baseline wcs D. Harbeck Jan 2015'
            imghdu.header['CTYPE1'] = 'RA---TAN'
            imghdu.header['CTYPE2'] = 'DEC--TAN'
            imghdu.header['CRVAL1'] = 0.0
            imghdu.header['CRVAL2'] = 0.0
            imghdu.header['CUNIT1']  = ('deg     ', 'Axis unit')
            imghdu.header['CUNIT2']  = ('deg     ', 'Axis unit')
            imghdu.header['CD1_1'] = 3.11588947907302E-5
            imghdu.header['CD2_2'] = 3.11503811254695E-5
            imghdu.header['CD1_2'] = -4.9413159234116E-8
            imghdu.header['CD2_1'] = -3.9298844602068E-8

            imghdu.header['CRPIX1'] = -4222.7484375 * xx + 15973.0199279 - 508.0 * cx
            imghdu.header['CRPIX2'] = -4315.55125 * yy + 12702.5552885 - 505.0 * cy

            imghdu.name = "OTA%d%d.SCI" % (xx,yy)
            
            otalist.append(imghdu)

    return otalist

            # print "OTA %i%i: WCSASTRM = 'baseline wcs D. Harbeck Jan 2015'" % (xx, yy)
            # print "OTA %i%i: CTYPE2 = 'DEC--ZPX'" % (xx, yy)
            # print "OTA %i%i: CD1_1 = 3.11588947907302E-5" % (xx, yy)
            # print "OTA %i%i: CD2_2 = 3.11503811254695E-5" % (xx, yy)
            # print "OTA %i%i: CD1_2 = -4.9413159234116E-8" % (xx, yy)
            # print "OTA %i%i: CD2_1 = -3.9298844602068E-8" % (xx, yy)
            
            # for cx in range(0, 7):
            #     for cy in range (0, 7):
            #         cr1 = -4222.7484375 * xx + 15973.0199279 - 508.0 * cx
            #         cr2 = -4315.55125 * yy + 12702.5552885 - 505.0 * cy
            #         print "OTA %i%i, EXTNAME xy%i%i: CRPIX1 = %.2f    CRPIX2 = %.2f" % (xx, yy, cx, cy, cr1, cr2)
           
            # otaid = "OTA %i%i" % (xx, yy)
            # newltvs = ltv.replace("OTA 61", otaid)
            # print newltvs

if __name__ == "__main__":
    
    otalist = printWCS()
    hdulist = pyfits.HDUList(otalist)
    hdulist.writeto(sys.argv[1])
