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

from podi_definitions import *
from podi_commandline import *
number_cpus = "auto"

catalog_directory = {}

#
#
# DO NOT MODIFY ANYTHING ABOVE THIS LINE
#
#



# Your setup parameters are next:
max_cpu_count = 7
wcs_ref_dir = ""
wcs_ref_type = "2mass_nir"
sdss_ref_type = "local"
sdss_ref_dir = ""
ucac4_ref_dir = ""
ippref_ref_dir = "none"
scratch_dir = "/scratch/"
sextractor = "/usr/bin/sex"
sex_redirect = " >& /dev/null"
sex_delete_tmps = True
diagplot__zeropoint_ZPrange = [0.5, 0.5]
diagplot__zeropoint_magrange = [11, 21]
diagplot__zeropointmap_range = [0.2, 0.2]
debug_log_filename = "/output/debug.log"
debug_log_append = True
crj_sigclip = 5.0
crj_sigfrac = 0.3
crj_objlim = 5.0
crj_saturation = 55000
crj_niter = 4
swarp_exec = "/usr/bin/swarp"
swarp_singledir = "/scratch/"
exec_dir = "/app/"
fixwcs_mode = "otashear"
max_pointing_error = [2, 5, 8, 12, 20]
max_rotator_error = [-3, 3.5]
min_wcs_quality = 3.0
wcs_match_multiplier = 0.6
saturation_limit = 65535
persistency_duration = 600
log_shell_output = False
staging_dir = "/scratch/"
sextractor_cache_dir = "/scratch/"
flat_order = ['tflat', 'dflat', 'flat']
per_ota_timeout = 300
mastercal_cache = "/mastercals/"
photcalib_saturation = 54000



# Catalog configuration
catalog_directory['ucac4'] = ("/catalogs/ucac4/", 5)
catalog_directory['sdss_dr8'] = ("/catalogs/sdss_dr8/", 5)
catalog_directory['ippref'] = ("/catalogs/ippref/", 5)
catalog_directory['gaia'] = ("/catalogs/gaia", 5)
catalog_directory['sdss_dr13'] = ("/catalogs/sdss_dr13/", 7)
catalog_directory['2mass'] = ("/catalogs/2mass/", 7)
catalog_directory['panstarrs'] = ("/catalogs/panstarrs/", 7)


# Add more options here

catalog_mags = {}
catalog_mags['sdss_dr8'] = ['ra', 'dec', 'u', 'err_u', 'g', 'err_g', 'r', 'err_r', 'i', 'err_i', 'z', 'err_z']
catalog_mags['ippref'] = ['ra', 'dec', ]
catalog_mags['2mass'] = ['ra', 'dec', 'err_ra', 'err_dec', ]
catalog_mags['sdss_dr13'] = ['ra', 'dec', 'u', 'err_u', 'g', 'err_g', 'r', 'err_r', 'i', 'err_i', 'z', 'err_z']
catalog_mags['ucac4'] = ['ra', 'dec', 'err_ra', 'err_dec', '', '', '', '', '', '',
                         'g', 'err_g', 'r', 'err_r', 'i', 'err_i',]
catalog_mags['panstarrs'] = ['ra', 'dec', 'err_ra', 'err_dec',
                             'g', 'err_g', 'r', 'err_r', 'i', 'err_i',
                             'z', 'err_z', 'y', 'err_y',]


# Add more options here




#
#
# DO NOT MODIFY ANYTHING BELOW THIS LINE
#
#

# This the order in which the catalogs are queried for the photometric 
# calibration. This is the currently best order, so mess with it on you 
# own risk. 
photcalib_order = ['panstarrs', 'sdss_dr13', 'ippref', 'ucac4']
photcalib_error_cutoff = {
    'SDSS': 0.02,
    'IPPRef': 0.05, 
    'UCAC4': 0.05,
}
wcscalib_order = ['gaia', 'sdss_dr13', '2mass']


if (cmdline_arg_isset("-ncpu")):
    number_cpus = int(cmdline_arg_set_or_default("-ncpu", number_cpus))
    # print "Using user-defined CPU count of",number_cpus
else:
    try:
        import multiprocessing
        number_cpus = multiprocessing.cpu_count()
        # print "Yippie, found %d CPUs to use in parallel!" % (number_cpus)
        if (number_cpus > max_cpu_count and max_cpu_count > 1):
            number_cpus = max_cpu_count
            # print "... but using only %d of them!" % (number_cpus)
    except:
        # print "Problem figuring out how many CPUs to use, setting to 1"
        number_cpus=1
        pass





################################################################################
#
#   ##
#   ##    Important!
#   ##    you will break the code!
#   ##
# 
#   ##    Do not put any options after this block!
#
################################################################################
if __name__ == "__main__":
    print
    print "This file defines some site properties, mostly paths."
    print
    print "Please run podi_testinstall to change settings in here"
    print "rather than editing the file by hand!"
    print
    
