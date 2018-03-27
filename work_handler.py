#!/usr/bin/env python

import os, sys
import multiprocessing
import queue
import threading
import time
import errno
import psutil

import podi_sitesetup as sitesetup
from podi_definitions import *
from podi_collectcells import *
import podi_logging
from podi_commandline import *

class reduce_collect_otas (object):

    def __init__(self, options, number_cpus):

        self.options = options
        self.number_cpus = number_cpus
        self.quit = False

        self.info = {}

        self.queue = multiprocessing.JoinableQueue()
        self.return_queue = multiprocessing.Queue()
        self.intermediate_queue = multiprocessing.Queue()
        self.shmem_dims = (4096, 4096)
        self.kw_worker_args= {
            'queue': self.queue,
            'return_queue': self.return_queue,
            'intermediate_queue': self.intermediate_queue,
            'options': self.options,
            'shmem': None,
            'shmem_dim': self.shmem_dims,
        }

        self.intermediate_results_done = multiprocessing.Lock()
        self.intermediate_results_done.acquire()


        self.files_to_reduce = []
        
        self.active_workers = 0

        #
        # Start feeding the workers
        #
        self.feed_worker_thread = threading.Thread(
            target=self.feed_workers,
            )

        #
        # Also start a worker to collect intermediate data
        #
        self.collect_intermediate_results_thread = threading.Thread(
            target=self.collect_intermediate_results
            )
        self.intermediate_results_complete = False

        self.intermediate_data_back_to_workers = None

        pass


    def start(self):
        self.feed_worker_thread.start()
        self.collect_intermediate_results_thread.start()

    def feed_workers(self):
        self.workers_started = 0
        should_be_working = []
        process_ids = []

        print "Starting to feed workers"
        x = 0
        while (not self.quit):
            #
            # Make sure all workers that should be working are doing so
            #
            workers_alive = 0
            for fn in self.info:
                process = self.info[fn]['process']
                if (process == None):
                    continue

                pid = process.pid
                if (pid == None):
                    # not started yet
                    continue

                ps = psutil.Process(pid)
                if (ps.status() in [psutil.STATUS_ZOMBIE,
                                    psutil.STATUS_DEAD]):
                    print "Found dead process: %d" % (pid)
                    process.terminate()
                    process.join(timeout=0.01)
                    self.info[fn]['process'] = None
                    self.active_workers -= 1
                    continue

                workers_alive += 1
            print "%d workers still alive" % (workers_alive)

            #
            # Start new workers if we have CPUs available
            #
            if (self.active_workers < self.number_cpus):
                print "we have some capacity to start new workers"
                # Check all workers, and start one if we find one that's not alive
                started_new_process = False
                for fn in self.files_to_reduce: #self.info:
                    if (self.info[fn]['process'] == None):
                        
                        #if (not self.info[fn]['process'].is_alive()):
                        print "starting worker for %s" % (fn)

                        p = multiprocessing.Process(target=parallel_collect_reduce_ota, 
                                                    kwargs=self.info[fn]['args'])

                        self.info[fn]['process'] = p
                        self.info[fn]['process'].start()

                        if (not self.info[fn]['intermediate_data'] == None):
                            # this process is being started after intermediate data
                            # has already been sent
                            # --> re-queue one more intermediate to allow completion
                            self.intermediate_queue.put(self.info[fn]['intermediate_data'])

                        self.active_workers += 1
                        #process_ids.append(self.info[fn]['process'].pid)
                        started_new_process = True
                        break
                if (started_new_process):
                    continue
            x += 1
            if (x%10 == 0): 
                print "still feeding workers"
                print ",".join(["%d" % (p) for p in process_ids])
            time.sleep(1)
        
    def collect_intermediate_results(self):
        print "Starting to collect intermediate results", len(self.info)
        self.intermediate_results_collected = 0

        while (self.intermediate_results_collected < len(self.info) and
               not self.quit):

            try:
                results = self.return_queue.get(timeout=0.1)
            except Queue.Empty:
                continue

            print "received some results!"
            self.intermediate_results_collected += 1
            self.active_workers -= 1

        print "***\n"*5,"All intermediate progress data received","\n***"*5
        self.intermediate_results_done.release()
        self.intermediate_results_complete = True

    def wait_for_intermediate_results(self):
        self.intermediate_results_done.acquire()
        self.intermediate_results_done.release()
        # Now we have all results

    def abort(self):
        self.quit = True
        if (not self.intermediate_results_complete):
            self.intermediate_results_done.release()
        print "Terminating feeder"
        #self.feed_worker_thread.terminate()
        for fn in self.info:
            try:
                print "terminating process for %s" % (fn)
                p = self.info[fn]['process']
                p.terminate()
                p.join()
                print "done!"
            except:
                pass

    def reduce_file(self, filename, id):

        if (not filename in self.info):
            self.info[filename] = {}
        else:
            print "we are already working on %s" % (filename)

        #
        # Setup a new process for this file
        #
        self.info[filename]['args'] = self.kw_worker_args
        self.info[filename]['args']['filename'] = filename
        self.info[filename]['args']['ota_id'] = id
        self.info[filename]['args']['shmem'] = multiprocessing.RawArray(ctypes.c_float, 4096*4096)
        self.info[filename]['args']['shmem_id'] = len(self.info)
        self.info[filename]['intermediate_data'] = None

        print "Setting up reduction for %s" % (filename)
        self.info[filename]['process'] = None

        self.files_to_reduce.append(filename)
        pass

        #, directory, filebase, otas):

        # print "Reducing these files:"
        # for ota in otas:
        #     filename = "%s/%s.%02d.fits" % (directory, filebase, ota)
        #     print filename


    
if __name__ == "__main__":

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)

    fn = sys.argv[1]
    directory, filebase = os.path.split(fn)
    filebase = filebase[:18]

    worker = reduce_collect_otas(options, 5)

    for ox,oy in itertools.product(range(2,5),repeat=2):
        ota = ox*10 + oy
        fn = "%s/%s.%02d.fits.fz" % (directory, filebase, ota)
        worker.reduce_file(fn, ota)

    worker.start()

    try:
        worker.wait_for_intermediate_results()
    except (SystemExit, KeyboardInterrupt):
        worker.abort()

    print "All data received!"
    worker.abort()

    podi_logging.shutdown_logging(options)
