#!/usr/bin/env python

"""

This module listens to SAMP message announcing new ODI frames. Each frame is
this automatically quick-reduced, either locally or on a remote-machine via
ssh. 

"""


try:
#    import sampy
    import xsampy as sampy
except ImportError:
    print "For this to work you need the SAMPy package installed"
    raise

import os
import sys
import time
import multiprocessing
import datetime

from podi_definitions import *
from podi_commandline import *
import podi_collectcells
import podi_focus

import podi_logging

import podi_SAMPsetup as setup

import subprocess
import logging

m = multiprocessing.Manager()
process_tracker = m.Queue()

worker_queue = multiprocessing.JoinableQueue()

metadata = {"samp.name":
                "QR_listener",
            "samp.description.text":
                "QuickReduce SAMP Listener",
            "samp.icon.url": 
                "file:///work/podi_devel/test/qr.jpg",
            "samp.documentation.url": 
                "http://members.galev.org/rkotulla/research/podi-pipeline/",
            "author.name": 
                "Ralf Kotulla",
            "author.email": 
                "kotulla@uwm.edu",
            "author.affiliation": 
                "University of Wisconsin - Milwaukee",
            "home.page": 
                "http://members.galev.org/rkotulla/research/podi-pipeline",
            "cli1.version":"0.01",
}


def worker_slave(queue):
    """

    This function handles all work, either running collectcells locally or 
    remotely via ssh. Files to reduce are read from a queue.

    """

    print "Worker process started, ready for action..."

    if (not setup.use_ssh):
        # If we reduce frames locally, prepare the QR logging.
        options = podi_collectcells.read_options_from_commandline()
        options = podi_logging.setup_logging(options)
        options['clobber'] = False

    while (True):
        try:
            # print "\n\nWaiting for stuff to do\n\n"
            task = queue.get()
        except KeyboardInterrupt, SystemExit:
            # print "worker received termination notice"
            # Ignore the shut-down command here, and wait for the official 
            # shutdown command from main task
            continue

        if (task == None):
            print "Shutting down worker"
            queue.task_done()
            break

        filename, object_name, obsid = task

        print "starting work on file",filename

        ccopts = ""
        if (len(sys.argv) > 2):
            # There are some parameters to be forwarded to collectcells
            ccopts = " ".join(sys.argv[1:])
        # print "ccopts=",ccopts

        if (cmdline_arg_isset("-dryrun")):
            print "Sending off file",filename,"for reduction"
            print "task done!"
            queue.task_done()
            continue


        if (object_name.lower().find("focus") >= 0):
            #
            # This is most likely a focus exposure
            #
            n_stars = int(cmdline_arg_set_or_default("-nstars", 7))

            if (setup.use_ssh):

                remote_inputfile = setup.translate_filename_local2remote(filename)
                kw = {
                    'user': setup.ssh_user,
                    'host': setup.ssh_host,
                    'filename': remote_inputfile,
                    'podidir': setup.remote_podi_dir,
                    'outdir': setup.output_dir,
                    'nstars': n_stars,
                }
                ssh_command = "ssh %(user)s@%(host)s %(podidir)s/podi_focus.py -nstars=%(nstars)d %(filename)s %(outdir)s" % kw

                process = subprocess.Popen(ssh_command.split(), 
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                _stdout, _stderr = process.communicate()

                print "*"*80
                if (not _stdout == "" and not _stdout == None):
                    print _stdout
                    print "*"*80
                if (not _stderr == "" and not _stderr == None):
                    print _stderr
                    print "*"*80

            else:
                # Run locally
                print "\nAnalyzing focus sequence (%s)\n" % (filename)
                podi_focus.get_focus_measurement(filename, n_stars=n_stars, output_dir=setup.output_dir)

            # Now check if we are supposed to open/display the focus plot
            if (not setup.focus_display == None):
                
                remote_filename = "%s/%s_focus.png" % (setup.output_dir, obsid)
                local_filename = setup.translate_filename_remote2local(filename, remote_filename)

                cmd = "%s %s &" % (setup.focus_display, local_filename)
                os.system(cmd)

        else:
            #
            # This is NOT a focus exposure
            #

            if (setup.use_ssh):

                # This is not a focus exposure, to treat it as a normal science exposure
                remote_inputfile = setup.translate_filename_local2remote(filename)
                kw = {
                    'user': setup.ssh_user,
                    'host': setup.ssh_host,
                    'collectcells': setup.ssh_executable,
                    'options': ccopts,
                    'filename': remote_inputfile,
                    'outputfile': setup.output_format, 
                }

                ssh_command = "ssh %(user)s@%(host)s %(collectcells)s %(filename)s %(outputfile)s %(options)s -noclobber" % kw

                process = subprocess.Popen(ssh_command.split(), 
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)
                _stdout, _stderr = process.communicate()

                print "*"*80
                if (not _stdout == "" and not _stdout == None):
                    print _stdout
                    print "*"*80
                if (not _stderr == "" and not _stderr == None):
                    print _stderr
                    print "*"*80

            else:
                print "\nRunning collectcells (%s)\n" % (filename)
                podi_collectcells.collectcells_with_timeout(input=filename, 
                                                            outputfile=setup.output_format,
                                                            options=options,
                                                            timeout=300,
                                                            process_tracker=process_tracker)

            #
            # If requested, also send the command to ds9
            #
            if (cmdline_arg_isset("-forward2ds9")):
                local_filename = setup.translate_filename_remote2local(filename, setup.output_format)
                forward2ds9_option = cmdline_arg_set_or_default("-forward2ds9", "image")
                if (forward2ds9_option == "irafmosaic"):
                    cmd = "mosaicimage iraf %s" % (local_filename)
                else:
                    cmd = "fits %s" % (local_filename)

                print "filename", filename
                print "remote file", remote_inputfile
                print "local file", local_filename

                try:
                    cli1 = sampy.SAMPIntegratedClient(metadata = metadata)
                    cli1.connect()
                    cli1.enotifyAll(mtype='ds9.set', cmd=cmd)
                    cli1.disconnect()
                except:
                    print "Problems sending message to ds9"
                    pass


        #
        # Once the file is reduced, mark the current task as done.
        #
        print "task done!"
        queue.task_done()

    if (not setup.use_ssh):
        print "Shutting down QuickReduce logging"
        podi_logging.shutdown_logging(options)

    print "Terminating worker process..."

    return
        


def get_filename_from_input(input):
    """

    Convert the input string, which can be either a FITS filename or a directory,
    into a valid FITS filename of one OTA of the exposure.

    """

    if (os.path.isfile(input)):
        return input
    elif (os.path.isdir(input)):
        if (input.endswith("/")):
            input = input[:-1]
        dirname, base = os.path.split(input)
        filename = "%s/%s.33.fits" % (input, base)
        if (not os.path.isfile(filename)):
            filename += ".fz"
            if (not os.path.isfile(filename)):
                return None
            return filename
        return filename

    return input
    
def check_obstype(filename):
    obstype, object_name, obsid = "", "", ""
    with pyfits.open(filename) as hdulist:
        obstype = hdulist[0].header['OBSTYPE']
        object_name = hdulist[0].header['OBJECT'] \
                      if 'OBJECT' in hdulist[0].header else ""
        obsid = hdulist[0].header['OBSID'] \
                if 'OBSID' in hdulist[0].header else ""

    return obstype, object_name, obsid


def receive_msg(private_key, sender_id, msg_id, mtype, params, extra):
    """
    
    This function is a callbakc handler that is called everytime a message
    is received from the SAMP hub.

    """


    #print "\n"*5,"new file received!\n"
    #print private_key, sender_id, msg_id, mtype, params, extra
    #cli1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
    #                    "samp.result": {"result": "ok guys"}})

    filename = params['filename']
    print "\nReceived command to reduce %s (at %s) ..." % (
        filename, datetime.datetime.now().strftime("%H:%M:%S.%f")
        )
    if (not os.path.isdir(filename)):
        print "filename %s is not a valid directory" % (filename)
        return
        
    fits_file = get_filename_from_input(filename)
    obstype, object_name, obsid = check_obstype(fits_file)

    if (obstype == "FOCUS" or 
        obstype == "OBJECT" or 
        not cmdline_arg_isset("-onlyscienceframes")):

        worker_queue.put( (filename, object_name, obsid) )

    else:
        
        print """

Received input %s
   (translated to %s) ...
This is not a OBJECT or FOCUS frame.
I was told to ignore these kind of frames.
\
""" % (filename, fits_file)
            return
            
#    worker_queue.put( (filename, object_name, obsid) )

    print "Done with this one, hungry for more!"
    return


def create_client(metadata, wait=0):

    # Create client, connect to Hub, and install message listener
    cli1 = sampy.SAMPIntegratedClient(metadata = metadata)

    try:
        cli1.connect()
    except sampy.SAMPHubError, sampy.SAMPClientError:
        if (wait>0): time.sleep(wait)
        return None
    except :
        print "some other problem with connecting"
        raise

    try:
        cli1.bindReceiveMessage(setup.message_queue, receive_msg)
    except:
        print "Problem with bindReceiveMessage"

    return cli1

def SAMPListener():
    print """

   *******************************************************************
   * SAMPListener for automatic image reduction (locally/remote)     *
   * Part of the QuickReduce package for the WIYN One Degree Imager  *
   * Author: Ralf Kotulla, kotulla@uwm.edu                           *
   *******************************************************************
"""

    # Create a client
    print "Starting receiver ..."

    print "Trying to connect to Hub ..."
    try:
        while (True):
            cli1 = create_client(metadata)
            if (cli1 == None):
                time.sleep(1)
            else:
                break

    except KeyboardInterrupt, SystemExit:
        print "\rAborting and shutting down SAMPListener ..."
        sys.exit(0)

    print "Starting execution process..."
    worker_process = multiprocessing.Process(target=worker_slave,
                                      kwargs={
                                          'queue': worker_queue,
                                      }
    )
    worker_process.start()

    print "Setup complete, waiting for messages..."
    quiet = cmdline_arg_isset("-quiet")
    quiet_string = "|/-\\"
    quiet_pos = 0

    try:
        while (True):

            # Ping the Hub
            try:
                ret = cli1.ecallAndWait("hub", "samp.app.ping", "5")
                # if successful, wait a little and check again
                time.sleep(1)
            except KeyboardInterrupt, SystemExit:
                #print "Received abort command, forwarding exception"
                raise
            except Exception as e:
                # If timeout expires than a SAMPProxyError is returned

                # Make sure to terminate the client to terminate the background 
                # thread. Otherwise the program won't close as it's waiting for 
                # the thread to finish.
                cli1.client.stop() 
                del cli1

                # This means that most likely there's no Hub (anymore)
                # Try to connect
                print "\nLost connection to hub, trying to re-establish it ..."
                while (True):
                    cli1 = create_client(metadata)
                    if (cli1 == None):
                        time.sleep(1)
                    else:
                        break

            if (not quiet):
                sys.stdout.write("\rCurrent system-time: %s (press Crtl-C to quit)" % (datetime.datetime.now()))
            else:
                quiet_pos += 1
                sys.stdout.write("%s\b" % (quiet_string[quiet_pos % 4]))
            sys.stdout.flush()
                
    except KeyboardInterrupt, SystemExit:
        #print "Got termination notice"
        pass

        
    print # to get off the "current system time" line
    try:
        # Disconnect from hub
        print "Disconnecting from SAMP Hub ..."
        cli1.disconnect()
    except:
        pass

    # Finish up left over work
    try:
        print "Finishing up work (%d jobs), please wait ..." % (worker_queue.qsize())
    except:
        print "Finishing up work, please wait ..."
    worker_queue.put(None)
    worker_process.join(1)

    print "All done, goodbye!"


if __name__ == "__main__":


    if (len(sys.argv) > 1 and sys.argv[1] == "-testconnect"):

        try:
            cli1 = sampy.SAMPIntegratedClient()
            cli1.connect()
            cli1.bindReceiveMessage(setup.message_queue, receive_msg)
            cli1.disconnect()

            print "\nConnection successful!\n"
        except:
            print "\nProblem connecting\n"
            pass
        
        sys.exit(0)


    elif (cmdline_arg_isset("-yappi")):
        print "Running with yappi profiler"
        import yappi
        yappi.start()
        SAMPListener()
        profiler = open("profile.data", "w")
        yappi.get_func_stats().debug_print() #print_all()
        yappi.get_thread_stats().print_all(out=profiler)
        profiler.close()

    else:

        
        SAMPListener()

    
