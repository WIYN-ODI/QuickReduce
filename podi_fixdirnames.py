#!/usr/local/bin/python
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



