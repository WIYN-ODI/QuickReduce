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

import podi_definitions
import podi_collectcells
import podi_SAMPsetup as setup

import subprocess

m = multiprocessing.Manager()
process_tracker = m.Queue()

worker_queue = multiprocessing.JoinableQueue()



def worker_slave(queue):
    """

    This function handles all work, either running collectcells locally or 
    remotely via ssh. Files to reduce are read from a queue.

    """

    print "Worker process started, ready for action..."

    while (True):
        try:
            task = queue.get()
        except KeyboardInterrupt, SystemExit:
            return

        if (task == None):
            queue.task_done()
            break

        filename = task

        # print "starting work on file",filename

        ccopts = ""
        if (len(sys.argv) > 2):
            # There are some parameters to be forwarded to collectcells
            ccopts = " ".join(sys.argv[1:])
        # print "ccopts=",ccopts

        if (setup.use_ssh):
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
            #print kw

            ssh_command = "ssh %(user)s@%(host)s %(collectcells)s %(filename)s %(outputfile)s %(options)s -noclobber" % kw
            #print ssh_command

            cmd_items = ssh_command.split()
            print "\nExecuting:\n%s\n" % (ssh_command)

            # Run ssh via a subprocess
            process = subprocess.Popen(cmd_items, stdout=subprocess.PIPE)
            _stdout, _stderr = process.communicate()
            #print "\n\nstdout=\n",_stdout
            #print "\n\nstderr=\n",_stderr
            #print ssh_command
            #print "done running ssh"

        else:
            #
            # Run collectcells locally
            #
            # print "Reading configuration for CollectCells"
            print "\nRunning collectcells (%s)\n" % (filename)
            options = podi_collectcells.read_options_from_commandline()
            options = podi_collectcells.setup_logging(options)
            options['clobber'] = False

            podi_collectcells.collectcells_with_timeout(input=filename, 
                                                        outputfile=setup.output_format,
                                                        options=options,
                                                        timeout=300,
                                                        process_tracker=process_tracker)

            podi_logging.podi_log_master_quit(options['log_master_info'])

        #
        # Once the file is reduced, mark the current task as done.
        #
        queue.task_done()

    return
        



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
        
    worker_queue.put( (filename) )

    print "Done with this one, hungry for more!"
    return



if __name__ == "__main__":

    print """

   *******************************************************************
   * SAMPListener for automatic image reduction (locally/remote)     *
   * Part of the QuickReduce package for the WIYN One Degree Imager  *
   * Author: Ralf Kotulla, kotulla@uwm.edu                           *
   *******************************************************************
"""

    # Create a client
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

    print "Starting receiver ..."

    # Create client, connect to Hub, and install message listener
    cli1 = sampy.SAMPIntegratedClient(metadata = metadata)
    cli1.connect()
    cli1.bindReceiveMessage("odi.image.load", receive_msg)

    print "Starting execution process..."
    worker_process = multiprocessing.Process(target=worker_slave,
                                      kwargs={
                                          'queue': worker_queue,
                                      }
    )
    worker_process.start()

    print "Setup complete, waiting for messages (press Crtl-C to quit)"

    try:
        while (True):
            time.sleep(2)
            sys.stdout.write("\rCurrent system-time: %s" % (datetime.datetime.now()))
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

#    
