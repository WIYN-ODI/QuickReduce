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

class async_fits_writer_thread (threading.Thread):
    def __init__(self, queue, queue_lock):
        threading.Thread.__init__(self)
        self.queue = queue
        self.queue_lock = queue_lock
    def run(self):
        print "Starting " + self.name

        while True:
            #self.queue_lock.acquire()
            #if not self.queue.empty():
            hdulist, filename, exit_cmd = self.queue.get(block=True)
            if (exit_cmd):
                #self.queue_lock.release()
                print "exiting worker!"
                self.queue.task_done()
                break
            print "Doing some work here!"
            time.sleep(2)
            self.queue.task_done()
            #self.queue_lock.release()
            #time.sleep(0.1)

        print "All done, going home!" 



class async_fits_writer():
    queue_lock = None
    fits_queue = None
    threads = []

    def __init__(self, number_threads=1):
        # Here, create all tools for the threaded fits writer
        self.queue_lock = threading.Lock()
        self.fits_queue = multiprocessing.JoinableQueue() #Queue.Queue()
        #self.threads = []

        # Create new threads
        for i in range(number_threads):
            thread = async_fits_writer_thread(self.fits_queue, self.queue_lock)
            thread.deamon = True
            thread.start()
            self.threads.append(thread)
       
    def write(self,hdulist, filename):
        #self.queue_lock.acquire()
        self.fits_queue.put((hdulist, filename, False))
        #self.queue_lock.release()

    def finish(self):
        print "Sending shutdown commands"
        for t in self.threads:
            #self.queue_lock.acquire()
            self.fits_queue.put((None, None, True))
            #self.queue_lock.release()
        for t in self.threads:
            t.join()
        print "Finishing up work"

    def __del__(self):
        self.finish()


if __name__ == "__main__":

    afw = async_fits_writer(3)
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    afw.write(None, "test1")
    del afw
    time.sleep(5)


