#! /usr/bin/env python

import sys
import os
import pyfits
import numpy
import scipy
import math
import subprocess
import time

pattern = [(  64,   48), 
           ( -34,  165), 
           ( -51,  -45), 
           (-131,  -70), 
           (-161, -113), 
           ( 240, -150),
           ( 189,   62), 
           ( 156,  168), 
           (-272,  -65),]



cmd_file = "/tmp/podi.cmd"
status_file = "/tmp/podi.status"

PODI_CMD_PAUSE = 1
PODI_CMD_RESUME = 2
PODI_CMD_ABORT =3

def check_for_command():
    if (os.path.isfile(cmd_file)):
        # We have a file that holds commands, so let's see what we have to do
        #print "We found the cmd file"
        file = open(cmd_file, "r")
        line = file.readline()
        if (len(line)>2):
            cmd = line.split()[0]
            return cmd
        
    return None

def reset_command():
    if (os.path.isfile(cmd_file)):
        file = open(cmd_file, "w")
        file.close()
    return

def handle_command(sum_dx, sum_dy):
    # Check if we received a command
    cmd = check_for_command()
    if (cmd == "pause"):
        reset_command()
        print """\
################
#
# Received pause command, waiting for resume or abort!
#
################
"""

        # Wait 1 second and check again
        counter = 0
        while( cmd != "abort" and cmd!= "resume" ):
            time.sleep(1)
            if (counter == 60):
                for i in range(counter):
                    sys.stdout.write("\b \b")
                counter = 0
            counter += 1
            sys.stdout.write("*")
            sys.stdout.flush()
            cmd = check_for_command()
            
        print """\
################
#
# Received resume command, continuing dither sequence!
#
################"""
        reset_command()

    if (cmd == "abort"):
        reset_command()
        # Return telescope to initial position and abort the dither sequence
        print """\
################
#
# Received abort command, undoing offsets and exiting script!
#
################"""
        offset_scope(-sum_dx, -sum_dy)
        sys.exit(0)

    return

        

def offset_scope(dx,dy):
    #
    # DISARM HERE
    #
    offset_cmd = 'echo target offset adjust %d\\",%d\\"' % (dx,dy)
    os.system(offset_cmd)
    return

if __name__ == "__main__":

#    if (len(sys.argv)<=2):
#        # Check if this is a command

    dither_pos_start = 1
    dither_pos_end = len(pattern)

    videxptime = 0.1
    preexptime = 1.0

    cmd = sys.argv[1]
    #print "reading command",cmd
    if (cmd == "pause" or
        cmd == "abort" or
        cmd == "resume"):
        print """\
#
# Sending '%s' command
#""" % (cmd)
        file = open(cmd_file, "w")
        print >>file, cmd
        file.close
        sys.exit(0)
    elif (cmd == "expose" or 
          cmd == "guide"):
        # check for parameters
        start = 2
        while(True):
            next = sys.argv[start]
            if (next[0] == "-"):
                param = next[1:]
                if (param == "start"):
                    dither_pos_start = int(sys.argv[start+1])
                    start += 2
                elif (param == "end"):
                    dither_pos_end = int(sys.argv[start+1])
                    start += 2
                elif (param == "vidtime"):
                    videxptime = float(sys.argv[start+1])
                    start += 2
                elif (param == "pretime"):
                    preexptime = float(sys.argv[start+1])
                    start += 2
                else:
                    print "unknown parameter:",next
            else:
                break
    else:
        print "Unknown command:",cmd
        sys.exit(0)

    reset_command()

    # Make sure dither positions are not negative
    dither_pos_start = numpy.max([dither_pos_start, 1])
    dither_pos_end   = numpy.max([dither_pos_end,   1])

    # Make sure dither positions do not exceed the number of pointings
    dither_pos_start = numpy.min([dither_pos_start, len(pattern)])
    dither_pos_end   = numpy.min([dither_pos_end,   len(pattern)])

    # Also make sure that start is before end
    if (dither_pos_start > dither_pos_end):
        dither_pos_start, dither_pos_end = dither_pos_end, dither_pos_start

    print "start-argv=",start
    if (len(sys.argv)-start<3):
        # Not enough parameters
        print """
To run this dither-sequence, run:
~/bin/podi_dither.py command "target name" exptime declination

command can be:
- expose - take image without guiding
- guide  - take image with guiding

The two commands support flags:
  -start x  -- start with dither position x instead of 1
  -end y    -- only execute up to dither position y
  -vidtime t -- change the video exposure time (default: 0.1)
  -pretime t -- change the pre-exposure exposure time (default: 1)

Additional commands:
- pause  - pause the script before the next commmand
- resume - resume script after a pause command
- abort  - terminate script and return telescope to origin

"""
        sys.exit(0)
    
    target_name = sys.argv[start]

    exptime = float(sys.argv[start+1])

    declination = float(sys.argv[start+2])

    print "Starting dither sequence"
    print "Object Name:  ", target_name
    print "Exposure time:", exptime
    print "Declination:  ", declination

    
    sum_dx, sum_dy = 0,0
#    for i in range(len(pattern)):
    for i in range(dither_pos_start, dither_pos_end+1):

        target_x, target_y = 0, 0
        for tmp in range(i):
            d_ra, d_dec = pattern[tmp]
            target_x += int(d_ra)
            target_y += int(d_dec)
        target_x *= math.cos(math.radians(declination))

        dx = target_x - sum_dx
        dy = target_y - sum_dy
        
        handle_command(sum_dx, sum_dy)
        print "#"
        print "#  Applying dither offset - position %d:" % (i)
        print "#     relative to last position: %+ 5d %+ 5d" % (dx, dy)
        print "#     relative to origin:        %+ 5d %+ 5d" % (target_x, target_y)
        print "#"
        offset_scope(dx, dy)

        sum_dx += dx
        sum_dy += dy

        handle_command(sum_dx, sum_dy)

        #
        # DISARM HERE
        #
        if (cmd == "guide"):
            ant_cmd = 'echo ant preimage -D_COMMENTS="%s" -D_EXPTIME=%f -D_PREEXPTIME=%f -D_VIDEXPTIME=%f' % (target_name, exptime, preexptime, videxptime)
        else:
            ant_cmd = 'echo ant expose_index -D_COMMENTS="%s" -D_EXPTIME=%f' % (target_name, exptime)
            
        #print ant_cmd
        sys.stdout.write("exposing")
        sys.stdout.flush()
        os.system(ant_cmd)
        time.sleep(exptime)
        sys.stdout.write("\n")
        sys.stdout.flush()

        
    # Return telescope to starting position
    offset_scope(-1*sum_dx, -1*sum_dy)


#os.system("ssh 192.168.1.10 afplay orch.wav")
