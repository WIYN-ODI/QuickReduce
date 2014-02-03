#!/usr/bin/env python


import os
import sys
import pyfits
import subprocess

from podi_collectcells import *


swarp_exec = "swarp"
tmp_dir = "/scratch/"
single_dir = "."


if __name__ == "__main__":
    
    # Figure out the config path
    abspath = os.path.abspath(sys.argv[0])
    dirname, filename = os.path.split(abspath)
    swarp_default = "%s/.config/swarp.default" % (dirname)

    stacked_output = get_clean_cmdline()[1]

    inputfiles = get_clean_cmdline()[2:]

    pixelscale = float(cmdline_arg_set_or_default("-pixelscale", 0))
    subtract_back = cmdline_arg_isset("-bgsub")

    reuse_singles = cmdline_arg_isset("-reusesingles")

    add_only = cmdline_arg_isset("-add") and os.path.isfile(stacked_output)
    if (add_only):
        print "Activating ADD mode"

    if (stacked_output.endswith(".fits")):
        stacked_output = stacked_output[:-5]

    header_only_file = "%s/preswarp.fits" % (tmp_dir)

    reference_file = cmdline_arg_set_or_default("-reference", None)
    if (not add_only):
        #
        # This is the regular start-from-scratch mode
        #

        # Set some Swarp options
        swarp_opts = """ \
               -IMAGEOUT_NAME %(imageout)s \
               -WEIGHTOUT_NAME %(weightout)s \
               -COMBINE_TYPE %(combine_type)s \
              """ % {
                  'imageout': header_only_file,
                  'weightout': "%s/preswarp.weight.fits" % (tmp_dir),
                  'combine_type': 'AVERAGE',
              }

        if (pixelscale > 0):
            swarp_opts += " -PIXELSCALE_TYPE MANUAL -PIXEL_SCALE %.4f " % (pixelscale)

        swarp_opts += " -SUBTRACT_BACK %s " % ("Y" if subtract_back else "N")

        print swarp_opts

        # 
        # First create only the output header so we can pass some information 
        # to the user
        #
        swarp_cmd = "%(swarp)s %(opts)s -HEADER_ONLY Y %(files)s" % {
            'swarp': swarp_exec,
            'opts': swarp_opts,
            'files': " ".join(inputfiles),
        }
        print swarp_cmd

        try:
            retcode = subprocess.Popen(swarp_cmd.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE) #, shell=True)
            # if retcode < 0:
            #     print >>sys.stderr, "Child was terminated by signal", -retcode
            # else:
            #     print >>sys.stderr, "Child returned", retcode
            #print retcode.stdout.readlines()
            #print retcode.stderr.readlines()
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e

        #
        # some information about the resulting stack is in the output-file
        #

        output_info = pyfits.open(header_only_file)
        print "Stack information..."
        print "   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header)

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']
    else:
        #
        # This is the simpler add-only mode
        #
        
        # Open the existing output header and get data from there
        output_info = pyfits.open(stacked_output+".fits")
        print "Stack information..."
        print "   Output-dimensions: %(NAXIS1)5d x %(NAXIS2)5d" % (output_info[0].header)

        out_crval1 = output_info[0].header['CRVAL1']
        out_crval2 = output_info[0].header['CRVAL2']
        out_naxis1 = output_info[0].header['NAXIS1']
        out_naxis2 = output_info[0].header['NAXIS2']

        

    
    #
    # Prepare the individual frames, rectified and re-projected 
    # to the final grid
    #
    single_prepared_files = []
    for singlefile in inputfiles:
        hdulist = pyfits.open(singlefile)
        obsid = hdulist[0].header['OBSID']

        # assemble all swarp options for that run
        dic = {'singledir': single_dir,
               'obsid': obsid,
               'pixelscale': pixelscale,
               'pixelscale_type': "MANUAL" if pixelscale > 0 else "MEDIAN",
               'center_ra': out_crval1,
               'center_dec': out_crval2,
               'imgsizex': out_naxis1,
               'imgsizey': out_naxis2,
               'resample_dir': tmp_dir,
               'inputfile': singlefile,
               'swarp_default': swarp_default,
           }

        swarp_opts = """\
-c $(swarp_default)s \
-IMAGEOUT_NAME %(singledir)s/%(obsid)s.fits \
-WEIGHTOUT_NAME %(singledir)s/%(obsid)s.weight.fits \
-PIXEL_SCALE %(pixelscale)f \
-PIXELSCALE_TYPE %(pixelscale_type)s \
-COMBINE Y \
-COMBINE_TYPE AVERAGE \
-CENTER_TYPE MANUAL \
-CENTER %(center_ra)f,%(center_dec)f \
-IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
-RESAMPLE_DIR %(resample_dir)s \
-SUBTRACT_BACK N \
%(inputfile)s \
""" % dic

        single_file = "%(singledir)s/%(obsid)s.fits" % dic

        # print swarp_opts
        swarp_cmd = "%s %s" % (swarp_exec, swarp_opts)
        print "Preparing file %s, please wait ..." % (singlefile)
        print swarp_cmd

        if ((add_only or reuse_singles) and os.path.isfile(single_file)):
            print "This single-swarped file (%s) exist, skipping it" % (single_file)
        else:
#        if (not reuse_singles or not os.path.isfile(single_file)):
            try:
                ret = subprocess.Popen(swarp_cmd.split(), 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)
                (swarp_stdout, swarp_stderr) = ret.communicate()
                print swarp_stdout
                print swarp_stderr
                #print "\n".join(swarp_stderr)
                single_prepared_files.append(single_file)
            except OSError as e:
                print >>sys.stderr, "Execution failed:", e
#        else:


    # If in "add" mode, rename the previous output file and add it to the list of input files
    if (add_only):

        if (len(single_prepared_files) < 1):
            print "No new files were added, so there's nothing to do."
            sys.exit(0)

        prev = 1
        while (True):
            filename = "%s.prev%02d.fits" % (stacked_output, prev)
            if (not os.path.isfile(filename)):
                break
            prev += 1
            continue
                
        # Rename the current output file and its weights
        old_stacked = "%s.prev%02d.fits" % (stacked_output, prev)
        old_weight = "%s.prev%02d.weight.fits" % (stacked_output, prev)
        os.rename(stacked_output+".fits", old_stacked)
        os.rename(stacked_output+".weight.fits", old_weight)

        # Also add the new re-named old stacked file to list of input files
        single_prepared_files.append(old_stacked)

    #
    # Now all single files are prepared, go ahead and produce the actual stack
    #
    dic['combine_type'] = "AVERAGE"
    dic['imageout'] = stacked_output+".fits"
    dic['weightout'] = stacked_output+".weight.fits"
    dic['prepared_files'] = " ".join(single_prepared_files)
    dic['bgsub'] = "Y" if subtract_back else "N"

    swarp_opts = """\
-c %(swarp_default)s \
-IMAGEOUT_NAME %(imageout)s \
-WEIGHTOUT_NAME %(weightout)s \
-COMBINE_TYPE %(combine_type)s \
-PIXEL_SCALE %(pixelscale)f \
-PIXELSCALE_TYPE %(pixelscale_type)s \
-COMBINE Y \
-COMBINE_TYPE %(combine_type)s \
-CENTER_TYPE MANUAL \
-CENTER %(center_ra)f,%(center_dec)f \
-IMAGE_SIZE %(imgsizex)d,%(imgsizey)d \
-RESAMPLE_DIR %(singledir)s \
-SUBTRACT_BACK %(bgsub)s \
-WEIGHT_TYPE MAP_WEIGHT \
%(prepared_files)s \
""" % dic

    print "\n\n\nStarting final stacking...\n\n\n"
    # print swarp_opts

    swarp_cmd = "%s %s" % (swarp_exec, swarp_opts)
    print swarp_cmd
    try:
        ret = subprocess.Popen(swarp_cmd.split(), 
                                   stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
        (swarp_stdout, swarp_stderr) = ret.communicate()
        print swarp_stdout
        print swarp_stderr
        #print "\n".join(swarp_stderr)
    except OSError as e:
        print >>sys.stderr, "Execution failed:", e

 
