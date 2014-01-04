#!/usr/bin/env python

"""

This module listens to SAMP message announcing new ODI frames. Each frame is
this automatically quick-reduced, either locally or on a remote-machine via
ssh. 

"""


try:
    import sampy
except ImportError:
    print "For this to work you need the SAMPy package installed"
    raise

import os
import sys
import time
import multiprocessing
import datetime

from podi_definitions import *
import podi_collectcells
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
        options = podi_collectcells.setup_logging(options)
        options['clobber'] = False

    while (True):
        try:
            # print "\n\nWaiting for stuff to do\n\n"
            task = queue.get()
        except KeyboardInterrupt, SystemExit:
            break
            #return

        if (task == None):
            queue.task_done()
            break

        filename = task

        print "starting work on file",filename

        ccopts = ""
        if (len(sys.argv) > 2):
            # There are some parameters to be forwarded to collectcells
            ccopts = " ".join(sys.argv[1:])
        # print "ccopts=",ccopts

        if (cmdline_arg_isset("-dryrun")):
            print "Sending off file",filename,"for reduction"

        elif (setup.use_ssh):
            #
            # Run collectcells on a different machine
            #

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

            cmd_items = ssh_command.split()
            # print "\nExecuting:\n%s\n" % (ssh_command)

            # Run ssh via a subprocess
            process = subprocess.Popen(cmd_items, stdout=subprocess.PIPE)
            _stdout, _stderr = process.communicate()

        else:
            #
            # Run collectcells locally
            #
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
        podi_logging.podi_log_master_quit(options['log_master_info'])

    print "Terminating worker process..."

    return
        


def get_filename_from_input(input):

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
    with pyfits.open(filename) as hdulist:
        obstype = hdulist[0].header['OBSTYPE']
    return obstype

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
    print "\nReceived command to reduce %s ..." % (filename)
    if (not os.path.isdir(filename)):
        print "filename %s is not a valid directory" % (filename)
        return
        
    if (cmdline_arg_isset("-onlyscienceframes")):
        fits_file = get_filename_from_input(filename)
        obstype = check_obstype(fits_file)
        if (obstype != "OBJECT"):
            print """

Received input %s
   (translated to %s) ...
This is not a OBJECT frame.
I was told to ignore non-OBJECT frames
\
""" % (filename, fits_file)
            return

    worker_queue.put( (filename) )

    print "Done with this one, hungry for more!"
    return



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

    # Create client, connect to Hub, and install message listener
    cli1 = sampy.SAMPIntegratedClient(metadata = metadata)
    cli1.connect()
    cli1.bindReceiveMessage(setup.message_queue, receive_msg)

    print "Starting execution process..."
    worker_process = multiprocessing.Process(target=worker_slave,
                                      kwargs={
                                          'queue': worker_queue,
                                      }
    )
    worker_process.start()

    print "Setup complete, waiting for messages..."

    try:
        while (True):
            time.sleep(2)
            sys.stdout.write("\rCurrent system-time: %s (press Crtl-C to quit)" % (datetime.datetime.now()))
            sys.stdout.flush()

            # try:
            #     cli1.connect()
            #     print "We were disconnected, connection re-established!"
            # except sampy.SAMPClientError:
            #     pass
            #     #print "We are connected!"
            # #print 

    except KeyboardInterrupt, SystemExit:
        pass

        
    try:
        # Disconnect from hub
        cli1.disconnect()
    except:
        pass

    # Finish up left over work
    print "\nFinishing up work (%d jobs)" % (worker_queue.qsize())
    worker_queue.put(None)
    worker_process.join()


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
    else:

        SAMPListener()

    
