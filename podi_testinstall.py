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


"""

podi_testinstall is a small tool that checks if all package dependencies at met.

See the podi-website at http://members.galev.org/rkotulla/research/podi-pipeline
for a full list of currently required packages.

"""



def check_package(name):
    """
    Try to import a package and print feedback message for the user whether or
    not the package has been found.
    """

    try:
        import_cmd = "import %s" % (name)
        exec(import_cmd)
    except ImportError:
        print "\nProblem importing %s" % (name)
    except:
        print "\nSome error occured while trying to import %s" % (name)
    else:
        print "Found working version of %s!" % (name)
        return True

    return False


if __name__ == "__main__":
    print
    print "Testing if all packages are installed"
    print

    print "\nchecking standard packages ..."
    check_package('os')
    check_package('sys')
    check_package('math')
    check_package('time')
    check_package('types')
    check_package('ctypes')
    check_package('itertools')

    print "\nchecking multi-processor packages ..."
    check_package('multiprocessing')
    check_package('Queue')
    check_package('threading')
    check_package('subprocess')


    print "\nchecking numerical processing packages ..."
    check_package('numpy')
    check_package('scipy')
    check_package('scipy.stats')
    check_package('scipy.optimize')
    check_package('scipy.interpolate')
    check_package('scipy.ndimage')
    check_package('bottleneck')

    print "\nchecking plotting packages ..."
    check_package('matplotlib')
    check_package('Image')
    check_package('ImageDraw')


    print "\nchecking astronomy-related packages ..."
    check_package('pyfits')
    check_package('ephem')
    check_package('astLib')
    check_package('pywcs')

    if (not check_package('podi_sitesetup')):
        print """\
        Module podi_sitesetup is a global configuration file for
        this podi pipeline. Copy the existing file
        podi_sitesetup.py.example to podi_sitesetup.py, open it
        in a text-editor and make sure the global settings for the 
        WCS and photometric reference catalogs are set correctly.
        Then re-run this program.

    """
