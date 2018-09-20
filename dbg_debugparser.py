#!/usr/bin/env python

import sys, os

if __name__ == "__main__":

    fulldebug_file = sys.argv[1]
    output_dir = sys.argv[2]

    df = open(fulldebug_file, "r")
    debug = df.readlines()
    
    start_found = False
    start = 0

    sublog_id = 0
    sublog_name = None

    n_success = 0
    n_failed = 0

    for i in range(len(debug)):

        if (debug[i].find("Starting --collectcells--")>=0 and not start_found):
            start_found = True
            start = i
        elif (debug[i].find("done with logging, closing file") >= 0 and start_found):
            sublog_id += 1

            print "Creating sublog",sublog_id, ";   current line =",i
            fn = "%s/debug.sublog.%04d" % (output_dir, sublog_id)
            if (sublog_name is not None):
                fn += ".%s" % (sublog_name)

            # Determine if the reduction was completed successfully
            success = False
            for dl in range(i,start,-1):
                if (debug[dl].find("Computing Ra/Dec for sky-samples in OTA") >= 0 or
                    debug[dl].find("All work completed successfully, output written to") >= 0):
                    success = True
                    break

            fn += ".%s" % ("SUCCESS" if success else "FAIL")
            if (success): n_success += 1
            if (not success): n_failed += 1

            f = open(fn, "w")
            f.write("".join(debug[start:i+1]))
            f.close()

            #start_found = True
            start = i+1
            sublog_name = None

            # Check next line and try to find what file we are dealing with
            try:
                next_line = debug[i+1]
                items = next_line.split()
                for it in range(len(items)):
                    if (items[it].find("collectcells.py")):
                        # This is the command, next file is the input file
                        inputfile = items[it+1]
                        _,exp = os.path.split(inputfile)
                        sublog_name = exp
                        break
            except IndexError:
                pass

    #print debug[:5]
    df.close()

    print "Total summary: %d successes, %d fails" % (n_success, n_failed)

    
