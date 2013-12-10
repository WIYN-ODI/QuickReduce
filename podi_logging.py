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
import podi_sitesetup as sitesetup


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

################################################################################
#
# Testing code goes here
#
################################################################################



# This is the worker process top-level loop, which just logs ten events with
# random intervening delays before terminating.
# The print messages are just so you know it's doing something!
def test_worker_process(log_setup):

#    configurer(queue)
#    log_setup
    name = multiprocessing.current_process().name
    print('Worker started xxx: %s' % name)
    podi_logger_setup(log_setup)

    #logger = podi_getlogger(name, log_setup)
    for i in range(10):
        time.sleep(random())
        print "in worker ."
#        logger = podi_getlogger(name, log_setup)
        logger = logging.getLogger(choice(LOGGERS))
        level = logging.DEBUG #choice(LEVELS)
        message = "msg %d: %s" % (i+1, choice(MESSAGES))
        # print message
        logger.log(level, message)
        del logger 
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






################################################################################
#
# Real code below
#
################################################################################



class QueueHandler(logging.Handler):
    """
    This is a logging handler which sends events to a multiprocessing queue.
    
    """

    def __init__(self, queue):
        """
        Initialise an instance, using the passed queue.
        """
        logging.Handler.__init__(self)

        #import os
        #print "\n\n\n\n Setting up logging, in Process ",os.getpid(),"XXX\n\n\n\n"

        self.queue = queue
        self.msgcount = 0

    def flush(self):
        pass

    def emit(self, record):

        """
        Emit a record.

        Writes the LogRecord to the queue.
        """
        self.msgcount += 1
        #sys.stdout.write("Current msg count: %d\n" % (self.msgcount))
        #sys.stdout.flush()

        # print "emitting 1 entry,",self.msgcount,"so far"
        try:
            #print "before adding one to queue",self.queue.qsize()
            #print "adding log entry to queue",self.msgcount, self.format(record)
            self.queue.put_nowait(record)
            #print "after adding one queue",self.queue.qsize()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            sys.stdout.write("OOppsie!\n")
            sys.stdout.flush()
            self.handleError(record)






def log_master(queue, options):
    """

    This is the main process that handles all log output. 

    Each log-entry is received via the queue that's being fed by all
    sub-processes, and then forwarded to other log-handlers.

    """

    # print "starting logging!"

    import sys
    root = logging.getLogger()
    h = logging.NullHandler() #StreamHandler(stream=sys.stdout)
    f = logging.Formatter('ROOTHANDLER %(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    h.setFormatter(f)
    root.addHandler(h)
    root.propagete = False

    enable_debug = False
    if (cmdline_arg_isset("-debugfile") or sitesetup.debug_log_filename != None):

        debug_filename = cmdline_arg_set_or_default("-debugfile", sitesetup.debug_log_filename)
        try:
            debugfile = open(debug_filename, "w")
            enable_debug = True
            # print 'activating debug output'

            debug_logger = logging.getLogger('debug')
            # debug_logger = logging.getLogger()
            h = logging.StreamHandler(stream=debugfile)
            f = logging.Formatter('DEBUGHANDLER %(processName)-10s %(asctime)s %(name)s %(levelname)-8s %(message)s')
            h.setFormatter(f)
            debug_logger.addHandler(h)
            debug_logger.propagate=False
        except:
            pass
    else:
        debug_logger = root

    info = logging.getLogger('info')
    h = logging.StreamHandler(stream=sys.stdout)
    f = logging.Formatter('%(name)s: %(message)s')
    h.setFormatter(f)
    info.addHandler(h)
    
    infolog_file = open("quickreduce.log", "w")
    h = logging.StreamHandler(stream=infolog_file)
    f = logging.Formatter('INFOHANDLER %(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    h.setFormatter(f)
    info.addHandler(h)
    info.propagate = False
            
    #
    # Check if we can connect to a RabbitMQ server
    #
    enable_pika = False
    try:
        debug_logger.debug("Trying to establish PIKA connection")
        import pika
        import podi_pikasetup
        credentials = pika.PlainCredentials(podi_pikasetup.user, podi_pikasetup.password)
        connection = pika.BlockingConnection(
            pika.ConnectionParameters( 
                credentials=credentials, 
                host=podi_pikasetup.host, 
                virtual_host=podi_pikasetup.vhost)
        )
        channel = connection.channel()
        channel.queue_declare(queue=podi_pikasetup.jobstatus_queue, durable=True)

        debug_logger.debug("PIKA connection established!")
        enable_pika = True
    except:
        debug_logger.debug("No PIKA connection available")
        pass
        
    
    msg_received = 0
    while True:
        try:
            record = queue.get()
            # print "received log entry",queue.qsize(),msg_received
            if record is None: # We send this as a sentinel to tell the listener to quit.
                break

            msg_received += 1
            # Add some logic here

            #print "record-level:",record.levelno, record.levelname, msg_received

            if (enable_debug):
                #print "handling at debug level", record.msg
                debug_logger.handle(record)

            if ((record.levelno > logging.DEBUG) and
                (record.levelno <= logging.INFO) ):
                info.handle(record)
                #print "handling at info level"
            elif (record.levelno > logging.INFO):
                # logger = logging.getLogger(record.name)
                # print "msg",msg_received," --> ",record.msg
                
                #print "handling at root level"
                info.handle(record) # No level or filter logic applied - just do it!

            #print "done with record.\n"

            if (enable_pika and 
                (record.levelno >= logging.INFO)):
                # only sent select message via Pika

                pika_msg = "%s: %s" % (record.name, record.msg)
                channel.basic_publish(exchange='',
                                      routing_key=podi_pikasetup.jobstatus_queue,
                                      properties=pika.BasicProperties(delivery_mode = 2),#persistent
                                      body=str(record.msg)
                )
             
            queue.task_done()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import sys, traceback
            print >> sys.stderr, 'Whoops! Problem:'
            traceback.print_exc(file=sys.stderr)

    if (enable_debug):
        print >>debugfile, "done with logging, closing file"
        debugfile.close()

    if (enable_pika):
        connection.close()





def podi_log_master_start(options):
    """
    
    This function creates the logging sub-process that handles all log output.

    This function also prepares the necessary information so we can activate the
    multiprocessing-safe logging in all sub-processes

    """

    queue = multiprocessing.JoinableQueue()
    listener = multiprocessing.Process(target=log_master,
                                kwargs={"queue": queue,
                                        "options": options}
                            )
    listener.start()

    worker_setup = {"queue": queue,
                  "configurer": log_slave_setup}
    
    log_master_info = {"queue": queue,
                       "listener": listener
                   }

    # Also start a logger for the main process
    podi_logger_setup(worker_setup)

    return log_master_info, worker_setup


def podi_log_master_quit(log_master_info):
    """
    Shutdown the logging process
    """

    log_master_info['queue'].put_nowait(None)
    try:
        log_master_info['listener'].join()
    except (KeyboardInterrupt, SystemExit):
        pass

    return


def podi_logger_setup(setup):
    """
    This function re-directs all logging output to the logging queue that feeds
    the logging subprocess.
    """
    
    if (setup == None):
        return
        # handler = logging.StreamHandler(sys.stdout)
    else:
        handler = QueueHandler(setup['queue'])

    # import sys
    # handler = logging.StreamHandler(stream=sys.stdout)
    logger = logging.getLogger()

    for h in logger.handlers:
        logger.removeHandler(h)

    logger.setLevel(logging.DEBUG)

    f = logging.Formatter('MYLOGGER = %(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    handler.setFormatter(f)
    
    logger.addHandler(handler)
    logger.propagate = True

    logger.debug("Started logging for process %s" % (multiprocessing.current_process().name))

    return



# def podi_getlogger(name, setup):

#     if (setup == None):
#         handler = logging.StreamHandler(sys.stdout)
#     else:
#         handler = QueueHandler(setup['queue'])

#     logger = logging.getLogger(name)
#     logger.addHandler(handler)
#     logger.setLevel(logging.DEBUG) # send all messages, for demo; no other level or filter logic applied.
#     logger.propagate = False

#     return logger


if __name__ == "__main__":

    if (cmdline_arg_isset("-pikalisten")):
        
        import pika
        import podi_pikasetup
        credentials = pika.PlainCredentials(podi_pikasetup.user, podi_pikasetup.password)
        connection = pika.BlockingConnection(
            pika.ConnectionParameters( 
                credentials=credentials, 
                host=podi_pikasetup.host, 
                virtual_host=podi_pikasetup.vhost)
        )
        channel = connection.channel()
        channel.queue_declare(queue=podi_pikasetup.jobstatus_queue, durable=True)

        print "PIKA connection established!"

        def callback(ch, method, properties, body):
            print " [o] %r" % (body,)

        channel.basic_consume(callback,
                              queue=podi_pikasetup.jobstatus_queue,
                              no_ack=True)

        print ' [*] Waiting for messages. To exit press CTRL+C'
        channel.start_consuming()

        sys.exit(0)




    import podi_collectcells
    options = podi_collectcells.read_options_from_commandline()

    log_master_info, log_setup = podi_log_master_start(options)

    workers = []
    # for i in range(10):
    #     worker = multiprocessing.Process(target=worker_process, kwargs=worker_log)
    #     workers.append(worker)
    #     worker.start()
    # for w in workers:
    #     w.join()

    print log_setup

    for i in range(1):
        worker = multiprocessing.Process(target=test_worker_process, 
                                         kwargs={"log_setup": log_setup})
#                                         args=(worker_log))
        workers.append(worker)
        worker.start()
    for w in workers:
        w.join()

    logger = podi_getlogger("main process", log_setup)

    logger.info('test info')
    logger.debug('test debug')
    logger.critical('test critical')
    logger.error('test error')


    podi_log_master_quit(log_master_info)
#    queue.put_nowait(None)
#    listener.join()
