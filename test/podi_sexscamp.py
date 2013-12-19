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

sex_params = {
    "CATALOG_TYPE":    "FITS_LDAC",
    "DETECT_THRESH":   5,
    "ANALYSIS_THRESH": 5,
    "FILTER":          "N",
}
# PARAMETERS_NAME
# CATALOG_NAME

scamp_params = {
    "ASTREF_CATALOG": "2MASS",
    "MATCH":          "N",
    "CROSSID_RADIUS": 3.0,
    "SOLVE_PHOTOM":   "N",
    "CHECKPLOT_DEV":  "NULL",
}


if __name__ == "__main__":

    dir, exe = os.path.split(sys.argv[0])
    param_file = "%s/.config/sex4scamp.param" % (dir)

    for filename in sys.argv[1:]:
        
        logfile = filename[:-5]+".log"

        sex_opt = ""
        for p in sex_params:
            sex_opt += " -%s %s" % (p, sex_params[p])

        catfile = filename[:-5]+".cat"
        sex_cmd = "sex -PARAMETERS_NAME %s -CATALOG_NAME %s %s %s >& %s" % (param_file, catfile, sex_opt, filename, logfile)

        #print sex_cmd
        print "Running SExtractor on frame",filename
        os.system(sex_cmd)


        scamp_opt = ""
        for p in scamp_params:
            scamp_opt += " -%s %s" % (p, scamp_params[p])

        print "Running Scamp on frame",filename
        scamp_cmd = "scamp %s %s >> %s 2>&1" % (scamp_opt, catfile, logfile)
        #print scamp_cmd
        os.system(scamp_cmd)

        print "done!\n"



