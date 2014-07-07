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

def ask_for_option(var_name, question, backup, needs_quotes, config_array):

    suggestion = backup
    if (hasattr(podi_sitesetup, var_name)):
        suggestion = eval('podi_sitesetup.%s' % (var_name))
        changed = False
    else:
        changed = True

    # Now ask user
    print "----"
    print question
    print "Default:", suggestion
    try:
        answer = raw_input("Answer:  ")
    except KeyboardInterrupt:
        print "\nTerminating\n"
        pass
        sys.exit(0)

    if (len(answer.strip()) <= 0):
        answer = suggestion
    
    # New configuration string
    if (needs_quotes):
        config_string = "%s = \"%s\"" % (var_name, answer)
    else:
        config_string = "%s = %s" % (var_name, answer)

    config_array.append(config_string)

    return (changed or (not answer == suggestion))

def update_sitesetup():

    import sys, os, podi_sitesetup
    import multiprocessing
    import datetime

    
    # Now go over all options, ask user for input
    # Use the current configuration as default value
    
    config_array = ["# Your setup parameters are next:"]

    changes = False
    changes = changes | ask_for_option('max_cpu_count', 
                   "Maximum number of CPU cores available for use in reduction", 
                   multiprocessing.cpu_count(), False, config_array)
    changes = changes | ask_for_option('wcs_ref_dir', 
                   "Directory holding the local 2MASS catalog in FITS format", 
                   "/some/dir", True, config_array)
    changes = changes | ask_for_option('wcs_ref_type', 
                   "Type of astrometric reference catalog (leave unchanged)", 
                   "2mass_nir", True, config_array)

    changes = changes | ask_for_option('sdss_ref_type', 
                   "Source of SDSS catalog (choose from local, web, stripe82)", 
                   "local", True, config_array)
    changes = changes | ask_for_option('sdss_ref_dir', 
                   "If local SDSS catalog, specify source directory", 
                   "/some/dir", True, config_array)

    changes = changes | ask_for_option('ucac4_ref_dir', 
                   "Directory of the UCAC4 directory - specify none if it does not exist", 
                   "none", True, config_array)

    changes = changes | ask_for_option('ippref_ref_dir', 
                   "Directory of the IPPRef directory - specify none if it does not exist", 
                   "none", True, config_array)

    changes = changes | ask_for_option('scratch_dir', 
                   "Scratch-directory for temporary files (the faster the better)", 
                   "/tmp", True, config_array)

    changes = changes | ask_for_option('sextractor', 
                   "Path to SourceExtractor executable (Hint: run 'which sex' in another terminal)", 
                   "/usr/local/bin/sex", True, config_array)
    changes = changes | ask_for_option('sex_redirect', 
                   "Command to re-direct stdout and stderr to /dev/null", 
                   " >& /dev/null", True, config_array)
    changes = changes | ask_for_option('sex_delete_tmps', 
                   "Delete SourceExtractor temporary files after use (choose from True (recommended) or False)", 
                   "True", False, config_array)

    changes = changes | ask_for_option('diagplot__zeropoint_ZPrange', 
                   "Vertical range around median in photometric calibation plots (format: [0.3,0.3]", 
                   "[0.5,0.5]", False, config_array)
    changes = changes | ask_for_option('diagplot__zeropoint_magrange', 
                   "SDSS Magnitude range in photometric calibation plots (format: [11, 22]", 
                   "[11, 21]", False, config_array)
    changes = changes | ask_for_option('diagplot__zeropointmap_range', 
                   "Spread around median ZP in photometric zeropoint map (format: [0.3,0.3]", 
                   "[0.2,0.2]", False, config_array)

    changes = changes | ask_for_option('debug_log_filename', 
                   "Filename for debug logs (e.g. /tmp/debug.log) ", 
                   "debug.log", True, config_array)
    changes = changes | ask_for_option('debug_log_append', 
                   "Keep adding to debug-log (choose True) or create a new file for each run (False)", 
                   "True", False, config_array)

    changes = changes | ask_for_option('crj_sigclip', 
                   "Cosmic Ray rejection: Sigma-clipping threshold", 
                   5.0, False, config_array)
    changes = changes | ask_for_option('crj_sigfrac', 
                   "Cosmic Ray rejection: Fraction of sigma-clipping for nearby neighbors", 
                   0.3, False, config_array)
    changes = changes | ask_for_option('crj_objlim', 
                   "Cosmic Ray rejection: Minimum CR significance above sources", 
                   5.0, False, config_array)
    changes = changes | ask_for_option('crj_saturation', 
                   "Cosmic Ray rejection: Saturation limit (negative disables this feature)", 
                   55000, False, config_array)
    changes = changes | ask_for_option('crj_niter', 
                   "Cosmic Ray rejection: Number of CR rejection iterations (typically in the range 1-4)", 
                   4, False, config_array)

    changes = changes | ask_for_option('swarp_exec', 
                   "Path to Swarp executable (Hint: run 'which swarp' in another terminal)", 
                   "swarp", True, config_array)
    changes = changes | ask_for_option('swarp_singledir', 
                   "Path for swarp intermediate files (you'll need several GB here!)", 
                   ".", True, config_array)

    changes = changes | ask_for_option('exec_dir', 
                   "Path where all the podi-scripts are installed", 
                   os.path.abspath("."), True, config_array)
    
    changes = changes | ask_for_option('max_pointing_error', 
                   "Maximum allowed pointing error for astrometric calibration", 
                   8.0, False, config_array)
    changes = changes | ask_for_option('max_rotator_error', 
                   "Maximum allowed range for rotator error for astrometric calibration", 
                   [-3,3.5], False, config_array)
    
    config_array.append("")

    # print config_array

    blank_setup = open("podi_sitesetup.py.blank", "r")
    lines = blank_setup.readlines()
    # Now find where to insert variables
    insert_at = -1
    for i in range(len(lines)):
        # print "____%s___" % (lines[i].strip())
        if (lines[i].startswith("###AUTO-CONFIG-INSERT-HERE")):
            # print "Found insert point"
            insert_at = i
            break

    if (changes):
        print "\n\nThere were some changes to the configuration!"
        while (True):
            answer = raw_input("Are you sure you want to update the configuration (y/n)? ")
            if (len(answer) > 0):
                break
            else:
                print "I really need an answer!",
        if (answer == "y" or answer == "Y"):
            backup_file = "podi_sitesetup.py.backup_from_%s" % (datetime.datetime.now().strftime("%Y%m%dT%H%M%S"))
            os.system("cp podi_sitesetup.py %s" % (backup_file))

            new_config = open("podi_sitesetup.py", "w")
            new_config.write("".join(lines[:insert_at]))
            new_config.write(os.linesep.join(config_array))
            new_config.write("".join(lines[insert_at+1:]))
            new_config.close()
            print "Changes saved!"
        else:
            print "Keeping configuration unchanged"
    else:
        print "\nNo changes found, keeping current configuration"


def check_component(pkg, fct):
    found = hasattr(pkg, fct)
    print "   * %-20s: %s" % (fct, found)
    return found

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

    print "\nChecking cython-optimized package for pODI"
    if (not check_package('podi_cython')):
        print """\
 There was a problem import podi_cython. This module contains optimized
 code that needs to be compiled first. To do so, simply run:

 python setup.py build_ext --inplace

"""
    else:
        # print "Checking podi_cython components:"
        import podi_cython
        all_found = True
        all_found = all_found and check_component(podi_cython, "sigma_clip_mean")
        all_found = all_found and check_component(podi_cython, "sigma_clip_median")
        all_found = all_found and check_component(podi_cython, "lacosmics")
        if (all_found):
            print "All routines found"
        else:
            print "Some podi-cython routines could not be found!"
            print "Please re-compile via python setup.py build_ext --python"
    
    print "\nCheck done!\n"

    answer = raw_input("Do you want to run the sitesetup assistant (y/N)?")
    if (answer.lower() == "y"):
        print "\n"*4,"     Starting auto-configuration!","\n"*4
        import sys, os, podi_sitesetup
        update_sitesetup()



