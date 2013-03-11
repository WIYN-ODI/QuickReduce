#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os



if __name__ == "__main__":

    for dirname in sys.argv[1:]:
        if (not os.path.isdir(dirname)):
            continue

        dirlist = os.listdir(dirname)

        # Check for the .33.fits file
        for filename in dirlist:
            if (filename.find("33.fits") > 0):
                pos = filename.find(".33.fits")
                basename = filename[:pos]
                #print dirname,"-->",basename

                if (os.path.isdir(basename)):
                    print "Directory",basename,"already exists, skipping!"
                else:
                    print "Renaming directory %s --> %s" % (dirname, basename)
                    mv_cmd = "mv %s %s" % (dirname, basename)
                    os.system(mv_cmd)

                break
    print "All done!\n"



