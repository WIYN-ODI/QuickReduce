#!/usr/bin/env python


import sys
import numpy
import os
import pyfits
import datetime
import scipy
import scipy.stats
import math
import scipy.spatial
import itertools
import logging

import time
import Queue
import threading
import multiprocessing

from  podi_definitions import *

# logging.basicConfig(filename='example.log',level=logging.DEBUG)








logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='/tmp/myapp.log',
                        filemode='w')


class ProgressConsoleHandler(logging.StreamHandler):
    """
    A handler class which allows the cursor to stay on
    one line for selected messages
    """
    #on_same_line = False
    def emit(self, record):
        try:
            msg = self.format(record)
            stream = self.stream
            levelname = record.levelname
            #print levelname
            # same_line = hasattr(record, 'same_line')
            # if self.on_same_line and not same_line:
            #     stream.write(self.terminator)
            #print stream
            stream.write(msg+"\n")
            # if same_line:
            #     stream.write('... ')
            #     self.on_same_line = True
            # else:
            #     stream.write(self.terminator)
            #     self.on_same_line = False
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)




class PikaMessageHandler(logging.Handler):
    """
    A handler class which allows the cursor to stay on
    one line for selected messages
    """
    #on_same_line = False
    def __init__(self, queue="hello"):
        logging.Handler.__init__(self)
        import pika
        self.connection = pika.BlockingConnection(pika.ConnectionParameters(
               'localhost'))
        self.channel = self.connection.channel()
        self.channel.queue_declare(queue='hello')
        print "Pika Handler all setup, ready for action"

    def emit(self, record):
        try:
            msg = self.format(record)
            #stream = self.stream
            levelname = record.levelname
            #print levelname
            # same_line = hasattr(record, 'same_line')
            # if self.on_same_line and not same_line:
            #     stream.write(self.terminator)
            #print stream
            #stream.write(msg+"\n")
            self.channel.basic_publish(exchange='',
                      routing_key='hello',
                      body=msg)
            print " [x] Sent ",msg

            # if same_line:
            #     stream.write('... ')
            #     self.on_same_line = True
            # else:
            #     stream.write(self.terminator)
            #     self.on_same_line = False
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)



log_queue = multiprocessing.Queue()


def log_reader(queue,x):

    msg_counter = 0

    while (True):
        msg = queue.get()
        # if (msg == None):
        #     break
        print msg_counter, msg
        msg_counter += 1

    return


class podi_logging():

    def __init__(self):
        p = multiprocessing.Process(target=log_reader, args=(log_queue,1))
        p.start()
        self.log_reader = p

    def __del__(self):
        self.log_reader.terminate()


# Start logging
# x = podi_logging()


def dump_to_log(id, idx):

    time.sleep(1)
    for i in range(100):
        txt = "dumping log-txt from id=%d --> %d\n" % (id,i)
        
        #print "dumping log1 from id=",id,"-->",i
        #print "dumping log2 from id=",id,"-->",i
        #print "dumping log3 from id=",id,"-->",i
        # print txt*5
#        time.sleep(numpy.random.random(1)*0.05)

        log_queue.put(txt)

    return


def podi_logger(name):

    logger = logging.getLogger() #multiprocessing.get_logger()
    logger.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')

    progress = ProgressConsoleHandler(stream=sys.stdout)
    progress.setFormatter(formatter)

    use_pika = False
    if (use_pika):
        pika = PikaMessageHandler(queue='hello')
        pika.setFormatter(formatter)
        logging.getLogger('').addHandler(pika)

    logger.addHandler(progress)

    logger.info("test")

    return logger

from random import choice, random
import time

LEVELS = [logging.DEBUG, logging.INFO, logging.WARNING,
          logging.ERROR, logging.CRITICAL]

LOGGERS = ['a.b.c', 'd.e.f']

MESSAGES = [
    'Random message #1',
    'Random message #2',
    'Random message #3',
]


class QueueHandler(logging.Handler):
    """
    This is a logging handler which sends events to a multiprocessing queue.
    
    The plan is to add it to Python 3.2, but this can be copy pasted into
    user code for use with earlier Python versions.
    """

    def __init__(self, queue):
        """
        Initialise an instance, using the passed queue.
        """
        logging.Handler.__init__(self)
        self.queue = queue
        
    def emit(self, record):
        """
        Emit a record.

        Writes the LogRecord to the queue.
        """
        try:
            ei = record.exc_info
            if ei:
                dummy = self.format(record) # just to get traceback text into record.exc_text
                record.exc_info = None  # not needed any more
            self.queue.put_nowait(record)
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)



# This is the worker process top-level loop, which just logs ten events with
# random intervening delays before terminating.
# The print messages are just so you know it's doing something!
def worker_process(queue, configurer):
    configurer(queue)
    name = multiprocessing.current_process().name
    print('Worker started: %s' % name)
    for i in range(10):
        time.sleep(random())
        logger = logging.getLogger(choice(LOGGERS))
        level = choice(LEVELS)
        message = choice(MESSAGES)
        logger.log(level, message)
    print('Worker finished: %s' % name)




# The worker configuration is done at the start of the worker process run.
# Note that on Windows you can't rely on fork semantics, so each process
# will run the logging configuration code when it starts.
def log_slave_setup(queue):
    h = QueueHandler(queue) # Just the one handler needed
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.DEBUG) # send all messages, for demo; no other level or filter logic applied.




def log_master_setup():
    root = logging.getLogger()
    # h = logging.handlers.RotatingFileHandler('/tmp/mptest.log', 'a', 300, 10)
    h = logging.StreamHandler(stream=sys.stdout)
    f = logging.Formatter('%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    h.setFormatter(f)
    root.addHandler(h)



# This is the listener thread top-level loop: wait for logging events
# (LogRecords)on the queue and handle them, quit when you get a None for a 
# LogRecord.
def log_master(queue, configurer):
    configurer()
    while True:
        try:
            record = queue.get()
            if record is None: # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record) # No level or filter logic applied - just do it!
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import sys, traceback
            print >> sys.stderr, 'Whoops! Problem:'
            traceback.print_exc(file=sys.stderr)


def podi_log_master_start():

    queue = multiprocessing.Queue(-1)
    listener = threading.Thread(target=log_master,
                                kwargs={"queue": queue,
                                        "configurer": log_master_setup}
                            )
    listener.start()

    worker_setup = {"queue": queue,
                  "configurer": log_slave_setup}
    
    log_master_info = {"queue": queue,
                       "listener": listener
                   }

    return log_master_info, worker_setup


def podi_log_master_quit(log_master_info):
    
    
    log_master_info['queue'].put_nowait(None)
    log_master_info['listener'].join()

    return


if __name__ == "__main__":

    log_master_info, log_setup = podi_log_master_start()

    workers = []
    worker_log = {"queue": log_queue,
                  "configurer": log_slave_setup}

    for i in range(10):
        worker = multiprocessing.Process(target=worker_process,
                                       kwargs=worker_log)
        workers.append(worker)
        worker.start()
    for w in workers:
        w.join()


    podi_log_master_quit(log_master_info)
#    queue.put_nowait(None)
#    listener.join()
