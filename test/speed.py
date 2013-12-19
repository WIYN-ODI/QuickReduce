#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#


import sys
import os
import pyfits
import numpy
import scipy
from multiprocessing import Process, Pipe, JoinableQueue

def multiprocess_f(conn):
    data = numpy.zeros(shape=(4096,40960))
    #conn.send([42, None, 'hello'])
    conn.send(data)
    conn.close()



def pipes():
    parent_conn, child_conn = Pipe()
    p = Process(target=multiprocess_f, args=(child_conn,))
    p.start()
    data = parent_conn.recv()   # prints "[42, None, 'hello']"
    print data.shape
    p.join()



import pickle
def ospipe_f(r,w):
    data = numpy.zeros(shape=(4096,40960))
    #conn.send([42, None, 'hello'])
    os.close(r)
    w = os.fdopen(w, 'w')
    #conn.send(data)
    #conn.close()
    pickle.dump(data, w)
    #w.write(data)
    w.close()

def os_pipes():
    r, w = os.pipe()
    p = Process(target=ospipe_f, args=(r,w,))
    p.start()
    os.close(w)
    r = os.fdopen(r)
    #data = r.read()   # prints "[42, None, 'hello']"
    data = pickle.load(r)
    r.close()
    print data.shape
    p.join()



def mpqueue_f(queue):
    data = numpy.zeros(shape=(4096,40960))
    #conn.send([42, None, 'hello'])
    queue.put(data)
    #conn.send(data)
    #conn.close()
    queue.task_done()

def mpqueue():
    queue = JoinableQueue()
    #parent_conn, child_conn = os.pipe()
    p = Process(target=mpqueue_f, args=(queue,))
    p.start()
    data = queue.get()   # prints "[42, None, 'hello']"
    print data.shape
    queue.join()
    p.join()


if __name__ == '__main__':
    import cProfile, pstats
    cProfile.run('pipes()', "profiler_pipes")
    p = pstats.Stats("profiler_pipes")
    p.sort_stats('time').print_stats()

#    cProfile.run('os_pipes()', "profiler_ospipes")
#    p = pstats.Stats("profiler_ospipes")
#    p.sort_stats('time').print_stats()

    cProfile.run('mpqueue()', "profiler_ospipes")
    p = pstats.Stats("profiler_ospipes")
    p.sort_stats('time').print_stats()





    #p.strip_dirs().sort_stats('time').print_stats()
  
