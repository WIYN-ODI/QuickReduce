#!/usr/bin/env python

import sys
import os
import numpy

from podi_definitions import get_collected_image_dimensions

size_x, size_y = get_collected_image_dimensions(1)

if __name__ == "__main__":
    
    output_file = sys.argv[-1]
    if (os.path.isfile(output_file)):
        print "output file",output_file,"exists"
        sys.exit(0)

    outfile = open(output_file, "w")
    print >> outfile, """\
# Region file format: DS9 version 4.1
global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical
"""

    for filename in sys.argv[1:-1]:
        #print filename

        ota = filename[-6:-4]
        otax = int(ota[0])
        otay = int(ota[1])
        print filename, ota, otax, otay 
        offset = numpy.array([otax * size_x, otay*size_y, 0, 0])
        print offset

        file = open(filename)
        for line in file:
            if (line[0:3] == "box"):
                coords = line[4:-2]
                coord_list = coords.split(",")
                      
                c = numpy.array([int(float(s)) for s in coord_list[:4]])
                print c
                # add offset
                c += offset

                print >>outfile, "-box(%d,%d,%d,%d)" % (c[0], c[1], c[2], c[3])
#                numpy.savetxt(outfile, c, "box(%d,%d,%d,%d)")

                #print coord_list

        file.close()

    outfile.close()
