#!/usr/bin/env python

print "importing sampy"
import sampy
import os
import sys
import time

message_queue = "odi.image.load"

if __name__ == "__main__":
    print "Starting SAMPY message sender"

    # Create a client
    metadata = {"samp.name":"QR_sender",
                "samp.description.text":"QuickReduce SAMP Sender",
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

    print "Sending message"
    cli1.enotifyAll(mtype=message_queue, filename=sys.argv[1])
    print "message sent"

    # time.sleep(5)

    cli1.disconnect()

#    
