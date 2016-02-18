#! /usr/bin/env python
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

Small tool to display the associations table in a collectcells-reduced frame.


Run as
-----------

``podi_showassociations.py (-full) input.fits``

By default, only the filenames of the associated input files are shown. With the
-full option, show_associations displays the full filenames of these files.

"""

import sys
import os
import pyfits
import numpy
from podi_definitions import *
from podi_commandline import *

import podi_logging
import logging

def read_associations(hdulist):

    logger = logging.getLogger("ReadAssocData")

    try:
        assoctable = hdulist['ASSOCIATIONS']
    except:
        logger.error("Couldn't find association table")
        return None

    reduction_steps = assoctable.data.field('correction')
    filenames = assoctable.data.field('filename_full')
    
    assocs = {}
    for idx in range(reduction_steps.shape[0]):
        if (not reduction_steps[idx] in assocs):
            assocs[reduction_steps[idx]] = []

        # We already have that step
        assocs[reduction_steps[idx]].append(filenames[idx])

    return assocs
    


def show_associations(inputfile, show_full_file):

    try:
        hdulist = pyfits.open(inputfile)
    except IOError:
        print "Can't open the input file: %s" % (inputfile)
        sys.exit(0)

    
    try:
        assoctable = hdulist['ASSOCIATIONS']
        # print assoctable
    except:
        stdout_write("Couldn't find association table, quitting\n\n")
        sys.exit(0)

    reduction_steps = assoctable.data.field('correction')
    if (show_full_file):
        filenames = assoctable.data.field('filename_full')
    else:
        filenames = assoctable.data.field('filename')

    # print reduction_steps

    stdout_write("\nAssociations for file %s ...\n" % (inputfile))
    stdout_write("-------------------------------\n")
    for step in set(reduction_steps):

        assoc_filenames = filenames[reduction_steps == step]

        for i in range(len(assoc_filenames)):
            if (i==0):
                stdout_write("% 15s : %s\n" % (step, assoc_filenames[i]))
            else:
                stdout_write("% 15s : %s\n" % (' ', assoc_filenames[i]))

    stdout_write("-------------------------------\n\n")

    return



def collect_reduction_files_used(masterlist, files_this_frame):
    """
    Keeps track of all files used during the reduction. This function maintains 
    a list of files and their corresponding reduction step, and properly handles
    individual files as well as list of filenames.
    
    This routine is called during reduction to keep track of all file 
    associations.

    Parameters
    ----------
    masterlist : dictionary

        Dictionary of all reduction steps that have external files associated. 
        For each step it maintains a list of files.

    files_this_frame : dictionary

        Like the master_list, just for only one step. 


    Returns
    -------
    the updated master_list

    """

    for key, value in files_this_frame.iteritems():
        if (key in masterlist):
            existing_keys = masterlist[key]
            if (type(value) == list):
                for val1 in value:
                    masterlist[key].append(val1)
            else:
                masterlist[key].append(value)
        else:
            # This is a new key, so just copy it
            if (type(value) == list):
                masterlist[key] = value
            else:
                masterlist[key] = [value]

    # print masterlist
    return masterlist










def create_association_table(master, verbose=False):
    """

    Convert the association dictionary maintained and updated by 
    :proc:collect_reduction_files_used and creates a FITS-compatible 
    associations table that will be stored in the resulting output FITS file.

    Parameters
    ----------
    master : dictionary

        the master associations dictionary 

    verbose : Bool

        Activate some debugging output

    Returns
    -------
    tbhdu : TableHDU

        A FITS TableHDU containing the assocations table. Each entry in this 
        table contains the reduction step, the name of the associated file, and 
        the full path of this file.

    """

    reduction_step = []
    full_filename = []
    short_filename = []

    for key, value in master.iteritems():
        # print key,":",value
        for filename in set(value):
            reduction_step.append(key)
            if (filename == None):
                continue
            full_filename.append(os.path.abspath(filename))
            
            dirname, filebase = os.path.split(filename)
            short_filename.append(filebase)

            if (verbose):
                print "% 15s : %s" % (key, filename)

    # print reduction_step
    # print short_filename

    columns = [\
        pyfits.Column(name='correction',    format='A25',  array=reduction_step),
        pyfits.Column(name='filename_full', format='A375', array=full_filename),
        pyfits.Column(name='filename',      format='A100', array=short_filename),
        ]

    coldefs = pyfits.ColDefs(columns)
    try:
        tbhdu = pyfits.BinTableHDU.from_columns(coldefs)
    except:
        tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    tbhdu.name = "ASSOCIATIONS"

    return tbhdu


if __name__ == "__main__":

    show_full_file = cmdline_arg_isset("-full")

    for inputfile in get_clean_cmdline()[1:]:
        show_associations(inputfile, show_full_file)
    

