#!/usr/bin/env python

from podi_definitions import *

wcs_ref_dir = "/datax/2mass_fits/"
wcs_ref_type = "2mass_nir"

number_cpus = "auto"
max_cpu_count = 6
if (cmdline_arg_isset("-ncpus")):
    number_cpus = int(cmdline_arg_set_or_default("-ncpus", number_cpus))
    print "Using user-defined CPU count of",number_cpus


#wcs_ref_dir = "/datax/2mass_fits/"
#wcs_ref_type = "2mass_opt"


if __name__ == "__main__":
    print "This file defines some site properties, mostly paths."
    print "Feel free to edit it to adapt directories to your site"
    
