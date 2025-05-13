#!/usr/bin/env python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#


import struct, datetime, pprint, sys, numpy, os
from podi_definitions import *

def my_float(s):
    try:
        return float(s)
    except:
        return 999.99 #numpy.nan
    return 0

# functions for converting input fields to usable data
cnv_text = lambda s: s.rstrip()
cnv_int = lambda s: int(s)
cnv_float = lambda s: my_float(s)
cnv_date_dmy = lambda s: datetime.datetime.strptime(s, "%d%m%Y") # ddmmyyyy
# etc

# field specs (field name, start pos (1-relative), len, converter func)
fieldspecs = [
    ('usnob1_id',   1,  12, cnv_text),
    ('ra',         15,  10, cnv_float),
    ('dec',        26,  10, cnv_float),
    ('bmag1',      77,   5, cnv_float),
    ('rmag1',      84,   5, cnv_float),
    ('bmag2',      91,   5, cnv_float),
    ('rmag2',      98,   5, cnv_float),
    ]

fieldspecs.sort(key=lambda x: x[1]) # just in case

# build the format for struct.unpack
unpack_len = 0
unpack_fmt = ""
for fieldspec in fieldspecs:
    start = fieldspec[1] - 1
    end = start + fieldspec[2]
    if start > unpack_len:
        unpack_fmt += str(start - unpack_len) + "x"
    unpack_fmt += str(end - start) + "s"
    unpack_len = end
field_indices = range(len(fieldspecs))
#print unpack_len, unpack_fmt
unpacker = struct.Struct(unpack_fmt).unpack_from

#print "CALCSIZE=",struct.calcsize(unpack_fmt)



def query_usno(ra, dec, radius, maxcount, catfile, download=True):

    if (ra < 0): ra += 360.0
        
    # Assemble the command string
    findusno = "findusnob1 %.8f %.8f -r %d -e b -m %d -sm > %s" % (ra, dec, radius, maxcount, catfile)
    print(findusno)
    stdout_write("Downloading USNO catalog from CDS/Vizier...\n")
    if (download):
        os.system(findusno)

    # Read the resulting catalog file
    file = open(catfile, "rb")
    raw_data = file.read()

    results = numpy.zeros(shape=(maxcount, 5))

    # Now interpret the resulting catalog line by line
    current_star = 0
    for line in file.readlines():

        # Ignore comment lines
        if (line[0] == "#"):
            continue
    

        raw_fields = unpacker(line[:struct.calcsize(unpack_fmt)])

        # Activate this for debugging
        #for x in field_indices:
        #    print fieldspecs[x][0], fieldspecs[x][3](raw_fields[x]),
        #print 

        # Put the data into the resulting array
        results[current_star,0] = fieldspecs[1][3](raw_fields[1])
        results[current_star,1] = fieldspecs[2][3](raw_fields[2])

        bmag = numpy.array([fieldspecs[3][3](raw_fields[3]), 
                            fieldspecs[5][3](raw_fields[5])])
        rmag = numpy.array([fieldspecs[4][3](raw_fields[4]), 
                            fieldspecs[6][3](raw_fields[6])])

        results[current_star,2] = numpy.min(bmag)
        results[current_star,3] = numpy.min(rmag)
        results[current_star,4] = numpy.min([bmag,rmag])
        current_star += 1

    return results[:current_star,:]



if __name__ == "__main__":

    results = query_usno(10.74225833333333, 41.37409777777778, 45, 1000, "usno.cat")
    print("Read",results.shape[0],"stars from USNO-B1")

    #for i in range(results.shape[0]):
    #    print results[i,0], results[i,1], results[i,2], results[i,3], results[i,4]

    #print results

