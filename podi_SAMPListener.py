#!/usr/bin/env python

print "importing sampy"
import sampy
import os
import sys
import time
import multiprocessing

import podi_definitions
import podi_collectcells
import podi_SAMPsetup as setup


m = multiprocessing.Manager()
process_tracker = m.Queue()

print "Reading configuration for CollectCells"
options = podi_collectcells.read_options_from_commandline()

def receive_msg(private_key, sender_id, msg_id, mtype, params, extra):

    print private_key, sender_id, msg_id, mtype, params, extra
    #cli1.reply(msg_id, {"samp.status": SAMP_STATUS_OK,
    #                    "samp.result": {"result": "ok guys"}})


    filename = params['filename']
    print "Received command to reduce %s ..." % (filename)
    if (not os.path.isdir(filename)):
        print "filename %s is not a valid directory" % (filename)
        return

    print "All read, starting work"

    podi_collectcells.collectcells_with_timeout(input=filename, 
                                                outputfile=setup.output_format,
                                                options=options,
                                                timeout=300,
                                                process_tracker=process_tracker)

    print "Done with this one, hungry for more!"
    return

if __name__ == "__main__":
    print "Starting SAMPY message sender"

    # Create a client
    metadata = {"samp.name":"QR_listener",
                "samp.description.text":"QuickReduce SAMP Listener",
                "samp.icon.url": "file:///work/podi_devel/test/qr.jpg",
                "cli1.version":"0.01"}

    cli1 = sampy.SAMPIntegratedClient(metadata = metadata)

    cli1.connect()

    # Construct the message
    # make sure this is compatible with ODIFileBrowser

    # private_key = cli1.getPrivateKey()
    # msg = {"samp.mtype": "odi.image.load",
    #        "samp.params": {"filename": sys.argv[1],
    #                        },
    # }
    # cli1.notifyAll(msg)
    
    # Start the logging for collectcells
    global options
    options = podi_collectcells.setup_logging(options)

    print "Starting receiver"
    cli1.bindReceiveMessage("odi.image.load", receive_msg)
    print "listener started"

    try:
        while (True):
            time.sleep(5)
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
        cli1.disconnect()
    except:
        pass


#    
