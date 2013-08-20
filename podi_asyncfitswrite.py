#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#
import sys
import os
import pyfits
import Queue
import traceback

import threading
import time
import multiprocessing
from podi_definitions import stdout_write, clobberfile

verbose = False
verbose = True

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
            #self.queue_lock.acquire()
            try:
                print "\n\nTrying to get some work",self.name
                hdulist, filename, exit_cmd = self.queue.get(block=True)

                print "\n\n",self.name,": ",filename,"\n\n"
                print "\n\n",self.name,": ",hdulist,"\n\n"
                
                #self.queue_lock.release()
            except:
                stdout_write("\n\n##############################\n#\n# Something terrible happened!\n")
                etype, error, stackpos = sys.exc_info()
                stdout_write("# Exception report:\n")
                stdout_write("#  ==> %s\n" % (error))
                print traceback.format_exc()
                stdout_write("#\n##############################\n")
                print "\n\n\nCaught problem, continuing\n\n\n"
                #self.queue_lock.release()
                continue
            
            if (exit_cmd or hdulist==None or filename==None):
                if (verbose): print "\n\nexiting worker!",self.name,"\n\n"
                self.queue.task_done()
                break

            
            if (verbose): print "Doing some work here!",filename
            #hdulist.fileinfo()
            #time.sleep(2)
            clobberfile(filename)
            hdulist.writeto(filename, clobber=True)
            if (True): stdout_write("File %s finished writing to disk\n" % (filename))
            self.queue.task_done()

        if (verbose): print "All done, going home!",self.name



class async_fits_writer():
    queue_lock = None
    fits_queue = None
    threads = []

    def __init__(self, number_threads=1):
        # Here, create all tools for the threaded fits writer
        self.queue_lock = threading.Lock()
        self.fits_queue = multiprocessing.JoinableQueue() #Queue.Queue()
        self.number_threads = number_threads

        self.start_threads()
       
    def write(self,hdulist, filename):
        if (verbose): stdout_write("Queued file %s for writing to disk.\n" % (filename))
        #self.queue_lock.acquire()
        self.fits_queue.put((hdulist, filename, False), False)
        #self.queue_lock.release()

    def start_threads(self):
        # Create new threads
        for i in range(1): #self.number_threads):
            thread = async_fits_writer_thread(self.fits_queue, self.queue_lock)
            thread.deamon = True
            thread.start()
            self.threads.append(thread)
       
    def wait(self):
        self.finish()
        self.start_threads()

    def finish(self, userinfo=False):
        if (userinfo): stdout_write("Waiting for asynchronous I/O to complete ...")
        if (verbose): print "Sending shutdown commands"
        for t in self.threads:
            self.fits_queue.put((None, None, True))
#        for t in self.threads:
#            print "Joining thread"
#            t.join()

        print "Joinging Queue"
        self.fits_queue.join()
        self.threads = []

        if (verbose): print "Finishing up work (in async_fits_writer.finish)"
        if (userinfo): stdout_write(" done with writing files!\n")

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
    print "Waiting"
    afw.wait()
    print "Starting new round"
    afw.write(hdulist, "deleteme.test7.fits")
    afw.write(hdulist, "deleteme.test8.fits")
    print "Done queueing all files for output!"
    del afw
    #time.sleep(5)


