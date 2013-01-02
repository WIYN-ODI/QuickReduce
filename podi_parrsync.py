#! /usr/bin/env python

import sys
import os
import subprocess as sp
import time

import Queue
import multiprocessing

def rsync(queue, target):
    while (True):
        shutdown, file = queue.get()
            
        if (shutdown):
            queue.task_done()
            return

        rsync_cmd = 'rsync -au "%s" %s' % (file, target)
            
        print rsync_cmd
        os.system(rsync_cmd)

        queue.task_done()


if __name__ == "__main__":
    if (len(sys.argv) <= 3):
        print "Usage:"
        print "./par_run.py (executable) (#CPUs) list of config files ..."
        sys.exit(0)

    target = sys.argv[1]
    number_cpus = int(sys.argv[2])

    queue = multiprocessing.JoinableQueue()
    processes = []
    for i in range(number_cpus):
        p = multiprocessing.Process(target=rsync, args=(queue,target))
        p.start()
        processes.append(p)

    file_list_starts_at = 3
    for file in sys.argv[file_list_starts_at:]:
        #if (not os.path.isfile()):
        #    continue

        queue.put((False,file))
        time.sleep(0.1)

    for i in range(number_cpus):
        queue.put((True,None))

    try:
        queue.join()
    except KeyboardInterrupt:
        for p in processes:
            p.terminate()


