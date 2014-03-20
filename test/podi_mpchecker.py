#!/usr/bin/env python

import os, sys, urllib2, urllib
import astropy.io.fits
import datetime
import ephem
import math
import numpy




def get_mpc_catalog(fitsfile):

    data = urllib.urlencode({
        'year': 2014,
        'month': 3,
        'day': 17.61,
        'which': 'pos',
        'ra': '13 22 17',
        'decl': '-15 33 16',
        'radius': 25,
        'TextArea': '',
        'limit': 20.0,
        'oc': 695,
        'sort': 'd',
        'mot': 'h',
        'tmot': 's',
        'pdes': 'u',
        'needed': 'f',
        'type': 'p',
        'ps': 'n',
        'submit': " Produce list ",
    })

    hdulist = astropy.io.fits.open(fitsfile)
    # print hdulist[0].header

    date_format = "%Y-%m-%dT%H:%M:%S.%f"
    date_obs = datetime.datetime.strptime(hdulist[0].header['DATE-OBS'], date_format)
    # print date_obs

    #print date_obs.year
    #print date_obs.month
    #print date_obs.day + date_obs.hour/24. + date_obs.minute/1440.

    ut_day = date_obs.day + date_obs.hour/24. + date_obs.minute/1440.

    ra = hdulist[0].header['CRVAL1']
    dec = hdulist[0].header['CRVAL2']
    j2k = ephem.Equatorial(math.radians(ra), math.radians(dec), epoch=ephem.J2000)
    print j2k.ra, j2k.dec
    print str(j2k.ra).replace(":", " ")
    print str(j2k.dec).replace(":", " ")

    data2 = [
        {'year': date_obs.year},
        {'month': date_obs.month},
        {'day': ut_day},
        {'which': 'pos'},
        {'ra': str(j2k.ra).replace(":", " ")},
        {'decl': str(j2k.dec).replace(":", " ")},
        {'TextArea': ''},
        {'radius': 40},
        {'limit': 23.0},
        {'oc': 695},
        {'sort': 'd'},
        {'mot': 'h'},
        {'tmot': 's'},
        {'pdes': 'u'},
        {'needed': 'f'},
        {'ps': 'n'},
        {'type': 'p'},
    ]

    url2 = ''
    for i in data2:
        url2 += urllib.urlencode(i)+"&"
        # print url2[:-1]


    web_server = "http://scully.cfa.harvard.edu"
    web_address = "/cgi-bin/mpcheck.cgi"

    headers = {"Content-type": "application/x-www-form-urlencoded", 
               "Accept": "text/plain",
               "Referer": "http://www.minorplanetcenter.net/cgi-bin/checkmp.cgi",
    }

    #u = urllib2.urlopen("http://scully/cfa/harvard.edu/cgi-bin/", data)
    #conn = u.request('POST', web_address, data, headers)

    # print data

    if (not os.path.isfile('dump')):
        req = urllib2.Request(web_server+web_address, url2, headers)
        response = urllib2.urlopen(req)
        the_page = response.read()

        dump = open("dump", "w")
        dump.write(the_page)
        dump.close()
    else:
        dump = open('dump', 'r')
        the_page = dump.read()
        dump.close()

    #print the_page

    # Search for the <pre> tag
    pre_start = the_page.find("<pre>")
    pre_end = the_page.find("</pre>")
    pre_block = the_page[pre_start:pre_end]

    lines = pre_block.split('\n')
    # print lines

    # The first 4 lines are the pre-tag and some header
    # The last line is empty
    datablock = lines[4:-1]

    print "\nResults found:\n---"
    print "\n".join(datablock)
    print "---"


    example = """
                                                                                                    1         1
          1         2         3         4         5         6         7         8         9         0         1
012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
(308501) 2005 TN158      13 27 50.6 -14 20 17  21.5   4.2E   0.0S     4+     1-    7o  None needed at this time.
(296241) 2009 CN60       13 27 47.3 -14 23 41  21.9   3.4E   3.4S     4+     5-    8o  None needed at this time.
(298661) 2004 CW32       13 27 44.9 -14 26 21  22.1   2.9E   6.1S     3+     8-    7o  None needed at this time.

(96709) 1999 JK127      13 21 58.0 -15 28 46  19.5   4.6W   4.5N    22-     9+   10o  None needed at this time.
(99996) 1981 EJ45       13 22 48.1 -15 27 06  19.9   7.5E   6.2N    20-     5+   16o  None needed at this time.
(92919) 2000 RR17       13 21 31.5 -15 38 05  18.8  10.9W   4.8S    23-    11+   11o  None needed at this time.
(143471)2003 CT5        13 23 35.6 -15 34 53  19.4  18.9E   1.6S    24-     0+   11o  None needed at this time.
(61081) 2000 LZ19       13 20 58.2 -15 29 14  18.5  19.0W   4.0N    16-    10+   12o  None needed at this time.
    """




    name, ra, dec, d_ra, d_dec, comment = [], [], [], [], [], []
    for obs in datablock:
        id = obs[:8].strip()

        _name = obs[9:25].strip() 
        if (len(id) <= 1 or id.find("(") < 0):
            _name = obs[:25].strip()

        _ra = obs[25:35].strip()
        _dec = obs[35:45].strip()
        rate_ra = obs[65:71].strip()
        rate_dec = obs[72:78].strip()

        _d_ra = float(rate_ra) if obs[71] == "+" else -1*float(rate_ra)
        _d_dec = float(rate_dec) if obs[78] == "+" else -1*float(rate_dec)

        # print "_%s_ _%s_ _%s_ _%s_ _%s_ _%s_ %.1f %.1f"  %(id, _name, _ra, _dec, rate_ra, rate_dec, _d_ra, _d_dec)

        _comment = obs[87:].strip()

        name.append(_name)
        ra.append(_ra)
        dec.append(_dec)
        d_ra.append(_d_ra)
        d_dec.append(_d_dec)
        comment.append(_comment)


    print name

    results = {
        'RA': numpy.array(ra),
        'DEC': numpy.array(dec),
        'Name': name,
        'dracosdec': numpy.array(d_ra),
        'ddec': numpy.array(d_dec),
        'comment': comment,
        }

    return results


if __name__ == "__main__":

    
    fitsfile = sys.argv[1]
    if (not os.path.isfile(fitsfile)):
        pass
        sys.exit(0)
        


        
    hdulist = astropy.io.fits.open(fitsfile)

    ra = hdulist[0].header['CRVAL1']
    dec = hdulist[0].header['CRVAL2']
    print ra, dec, ra/15.
    j2k = ephem.Equatorial(math.radians(ra), math.radians(dec), epoch=ephem.J2000)

    print j2k.ra, j2k.dec

    print str(j2k.ra)
    print str(j2k.ra).replace(":", " ")
    print str(j2k.dec).replace(":", " ")

    # print hdulist[0].header

    name, ra, dec, d_ra, d_dec = get_mpc_catalog(fitsfile)
