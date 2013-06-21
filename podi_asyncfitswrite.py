#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#
import sys
import os
import pyfits
import Queue

import threading
import time
import multiprocessing
from podi_definitions import stdout_write, clobberfile

verbose = False

#
# This is the threaded class that does the actual work
#
class async_fits_writer_thread (threading.Thread):
    # Some initialization stuff to remember all variables
    def __init__(self, queue, queue_lock):
        threading.Thread.__init__(self)
        self.queue = queue
        self.queue_lock = queue_lock

    # This is the end-less loop checking the Queue and doing the work
    def run(self):
        if (verbose): print "Starting " + self.name

        while True:
            hdulist, filename, exit_cmd = self.queue.get(block=True)
            if (exit_cmd):
                if (verbose): print "exiting worker!"
                self.queue.task_done()
                break

            
            if (verbose): print "Doing some work here!"
            #time.sleep(2)
            clobberfile(filename)
            hdulist.writeto(filename, clobber=True)
            stdout_write("File %s finished writing to disk\n" % (filename))
            self.queue.task_done()

        if (verbose): print "All done, going home!" 



class async_fits_writer():
    queue_lock = None
    fits_queue = None
    threads = []

    def __init__(self, number_threads=1):
        # Here, create all tools for the threaded fits writer
        self.queue_lock = threading.Lock()
        self.fits_queue = multiprocessing.JoinableQueue() #Queue.Queue()

        # Create new threads
        for i in range(number_threads):
            thread = async_fits_writer_thread(self.fits_queue, self.queue_lock)
            thread.deamon = True
            thread.start()
            self.threads.append(thread)
       
    def write(self,hdulist, filename):
        stdout_write("Queued file %s for writing to disk.\n" % (filename))
        self.fits_queue.put((hdulist, filename, False))

    def finish(self):
        if (verbose): print "Sending shutdown commands"
        for t in self.threads:
            self.fits_queue.put((None, None, True))
        for t in self.threads:
            t.join()
        if (verbose): print "Finishing up work"

    def __del__(self):
        self.finish()


if __name__ == "__main__":

    import numpy
    zeros = numpy.zeros(shape=(5000,5000), dtype=numpy.float32)
    hdu = pyfits.PrimaryHDU(data=zeros)
    hdulist = pyfits.HDUList([hdu])
    
    afw = async_fits_writer(3)
    afw.write(hdulist, "deleteme.test1.fits")
    afw.write(hdulist, "deleteme.test2.fits")
    afw.write(hdulist, "deleteme.test3.fits")
    afw.write(hdulist, "deleteme.test4.fits")
    afw.write(hdulist, "deleteme.test5.fits")
    afw.write(hdulist, "deleteme.test6.fits")
    afw.write(hdulist, "deleteme.test7.fits")
    afw.write(hdulist, "deleteme.test8.fits")
    print "Done queueing all files for output!"
    del afw
    #time.sleep(5)


