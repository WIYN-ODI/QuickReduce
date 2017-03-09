#!/usr/bin/env python

"""

This module contains all routines to configure and start the multi-processing
safe logging. All log output is queued, and handled in sequence by a separate
logging process.

"""


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
import traceback

from podi_definitions import *
from podi_commandline import *
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

    name = multiprocessing.current_process().name
    print('Worker started xxx: %s' % name)
    podi_logger_setup(log_setup)

    #logger = podi_getlogger(name, log_setup)
    for i in range(10):
        time.sleep(random())
        # print "in worker ."
        logger = logging.getLogger(choice(LOGGERS))
        level = choice(LEVELS)
        message = "msg %d: %s" % (i+1, choice(MESSAGES))
        # print message
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
    try:
        h = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        h = logging.StreamHandler(strm=sys.stdout)
    except:
        raise
    f = logging.Formatter('%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
    h.setFormatter(f)
    root.addHandler(h)




PPA_PROGRESS_LEVEL = 19

def ppa_update_progress(progress, message):

    jobid = int(cmdline_arg_set_or_default("-jobid", 0))
    logical_id = cmdline_arg_set_or_default("-logicalid", "???")

    logger = logging.getLogger("PPAUpdate(%s)" % logical_id)

    msg = "<Request><JobID>%(job_id)s</JobID><LogicalID>%(logical_id)s</LogicalID><Status>RUNNING:%(progress)02d</Status><Note>%(note)s</Note></Request>" % {
        "job_id": jobid,
        "logical_id": logical_id,
        "progress": progress,
        "note": message,
    }

    logger.log(PPA_PROGRESS_LEVEL, msg) 

    return


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
        except AssertionError:
            pass
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
    try:
        h = logging.NullHandler() #StreamHandler(stream=sys.stdout)
        f = logging.Formatter('%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
        h.setFormatter(f)
        root.addHandler(h)
    except AttributeError:
        # This happens in older Python versions that don't have a NULLHandler
        pass
    except:
        raise

    root.propagate = False

    enable_debug = False
    debug_logger = root
    if (cmdline_arg_isset("-debugfile") or sitesetup.debug_log_filename != None):

        debug_filename = cmdline_arg_set_or_default("-debugfile", sitesetup.debug_log_filename)
        try:
            open_mode = "a" if (sitesetup.debug_log_append) else "w"
            debugfile = open(debug_filename, open_mode)
            enable_debug = True
            print >>debugfile, " ".join(sys.argv)
            
            # print 'activating debug output'

            debug_logger = logging.getLogger('debug')
            # debug_logger = logging.getLogger()
            try:
                h = logging.StreamHandler(stream=debugfile)
            except TypeError:
                h = logging.StreamHandler(strm=debugfile)
            except:
                raise
            f = logging.Formatter('%(asctime)s -- %(levelname)-8s [ %(filename)30s : %(lineno)4s - %(funcName)30s() in %(processName)-12s] %(name)30s :: %(message)s')
            h.setFormatter(f)
            debug_logger.addHandler(h)
            debug_logger.propagate=False
        except:
            print "#@#@#@#@#@# Unable to write to debug file: %s" % (debug_filename)
            print "#@#@#@#@#@# Routing all debug output to stderr"
            debug_logger = logging.getLogger('debug')
            try:
                h = logging.StreamHandler(stream=sys.stderr)
            except TypeError:
                h = logging.StreamHandler(strm=sys.stderr)
            except:
                raise
            f = logging.Formatter('%(asctime)s -- %(levelname)-8s [ %(filename)30s : %(lineno)4s - %(funcName)30s() in %(processName)-12s] %(name)30s :: %(message)s')
            h.setFormatter(f)
            debug_logger.addHandler(h)
            debug_logger.propagate=False
            
            pass

    #
    # Create a handler for all output that also goes into the display
    #
    info = logging.getLogger('info')
    
    # set format for the terminal output
    try:
        h = logging.StreamHandler(stream=sys.stdout)
    except TypeError:
        h = logging.StreamHandler(strm=sys.stdout)
    except:
        raise
    # Add some specials to make sure we are always writing to a clean line
    # f = logging.Formatter('\r\x1b[2K%(name)s: %(message)s')
    f = logging.Formatter('%(name)s: %(message)s')
    h.setFormatter(f)
    info.addHandler(h)
    
    # 
    # Also write all info/warning/error messages to the logfile
    #
    infolog_file = open("quickreduce.log", "w")
    try:
        h = logging.StreamHandler(stream=infolog_file)
    except TypeError:
        h = logging.StreamHandler(strm=infolog_file)
    except:
        raise
    f = logging.Formatter('%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s')
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

        amqp_uri = os.getenv('AMQP_URI')
        if (not amqp_uri == None):
            debug_logger.debug("Using credentials from AMQP_URI")
            params = pika.connection.URLParameters(amqp_uri) 
            connection = pika.BlockingConnection(params) 
        else:
            debug_logger.debug("Using credentials from pika_setup file")
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
            try:
                record = queue.get(timeout=1.)
            except (KeyboardInterrupt, SystemExit):
                record = None
            except Queue.Empty:
                pass
                # print >>debugfile, "LogHandler: still running, but no message during the last second!"
                # print "."
                continue
            except:
                raise

            if (record == None): 
                break

            msg_received += 1
            # Add some logic here

            #print "record-level:",record.levelno, record.levelname, msg_received

            if (enable_debug):
                #print "handling at debug level", record.msg
                debug_logger.handle(record)

            if (record.levelno == PPA_PROGRESS_LEVEL): # and enable_pika):
                # print "\n"*10, "sending update to PPA:\n"+str(record.msg),"\n"*5
                # This message is sent via the update_ppa_progress command
                try:
                    channel.basic_publish(
                        exchange='',
                        routing_key=podi_pikasetup.jobstatus_queue,
                        properties=pika.BasicProperties(delivery_mode = 2),#persistent                                                           
                        body=str(record.msg)
                    )
                    # pass
                except:
                    pass
            elif ((record.levelno > logging.DEBUG) and
                (record.levelno <= logging.INFO) ):
                info.handle(record)
                #print "handling at info level"
            elif (record.levelno > logging.INFO):
                # logger = logging.getLogger(record.name)
                # print "msg",msg_received," --> ",record.msg
                
                #print "handling at root level"
                info.handle(record) # No level or filter logic applied - just do it!

                if ('SCA_PROGRESS_URL' in os.environ):
                    import requests
                    requests.post(os.environ['SCA_PROGRESS_URL']+".info",
                                  json={"msg": record.msg})

            #print "done with record.\n"

            #
            # only sent select message via Pika, and only if Pika has been 
            # initialized successfully
            #
            if (enable_pika and 
                (record.levelno >= logging.INFO)):

                pika_msg = "%s: %s" % (record.name, record.msg)
                try:
                    # use the optional formating routine to make the message
                    # compatible with e.g. the IU-PPA system
                    pika_msg = podi_pikasetup.format_msg(record)
                except:
                    pass
                    
                channel.basic_publish(exchange='',
                                      routing_key=podi_pikasetup.jobstatus_queue,
                                      properties=pika.BasicProperties(delivery_mode = 2),#persistent
                                      body=str(pika_msg)
                )
             
            queue.task_done()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            import sys, traceback
            print >> sys.stderr, 'Whoops! Problem:'
            traceback.print_exc(file=sys.stderr)

    if (enable_pika):
        try:
            connection.close()
        except:
            pass

    if (enable_debug):
        print >>debugfile, "done with logging, closing file"
        debugfile.close()





def podi_log_master_start(options):
    """
    
    This function creates the logging sub-process that handles all log output.

    This function also prepares the necessary information so we can activate the
    multiprocessing-safe logging in all sub-processes

    """

    queue = multiprocessing.JoinableQueue()

    #
    # Rename the thread name to make a more useful stacktrace
    #
    queue._start_thread()
    queue._thread.name = "QueueFeederThread___ParallelLogging"

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

    # print_stacktrace()

    return log_master_info, worker_setup


def podi_log_master_quit(log_master_info):
    """
    Shutdown the logging process
    """

    log_master_info['queue'].put(None)
    try:
        # print "joining log listener"
        log_master_info['listener'].join()
        # print "done joining log listener"
    except (KeyboardInterrupt, SystemExit):
        pass

    log_master_info['queue'].close()
    log_master_info['queue'].join_thread()
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


def log_exception(name=None):

    etype, error, stackpos = sys.exc_info()

    exception_string = ["\n",
                        "=========== EXCEPTION ==============",
                        "etype: %s" % (str(etype)),
                        "error: %s" % (str(error)),
                        "stackpos: %s" % (str(stackpos)),
                        "---\n",
                        traceback.format_exc(),
                        "--- end\n"
    ]
    logger = logging.getLogger(name)
    logger.critical("\n".join(exception_string))
    return


def log_platform_debug_data():
    logger = logging.getLogger("PLATFORM")
    try:
        import platform
        logger.debug("Python version: %s" % (str(platform.python_version())))
        logger.debug("Python compiler: %s" % (str(platform.python_compiler())))
        logger.debug("Python build: %s" % (str(platform.python_build())))

        logger.debug("OS version: %s" % (str(platform.platform())))

        logger.debug("OS uname: %s" % (" ".join(platform.uname())))
        logger.debug("OS system: %s" % (str(platform.system())))
        logger.debug("OS node: %s" % (str(platform.node())))
        logger.debug("OS release: %s" % (str(platform.release())))
        logger.debug("OS version: %s" % (str(platform.version())))
        logger.debug("OS machine: %s" % (str(platform.machine())))
        logger.debug("OS processor: %s" % (str(platform.processor())))

        logger.debug("interpreter: %s" % (" ".join(platform.architecture())))
    except:
        logger.debug("OS info not available, missing package platform")
        pass

    try:
        import socket
        logger.debug("Socket hostname: %s" % (socket.gethostname()))
    except:
        logger.debug("socket info not available, missing package socket")
        pass

    try:
        import getpass
        logger.debug("username: %s" % (getpass.getuser()))
    except:
        logger.debug("username not available, missing package getpass")
        pass
        
    return


def setup_logging(options=None):

    if (options == None):
        options = {}

    # Setup everything we need for logging
    log_master_info, log_setup = podi_log_master_start(options)
    options['log_setup'] = log_setup
    options['log_master_info'] = log_master_info

    log_platform_debug_data()

    return options
    
def shutdown_logging(options):
    podi_log_master_quit(options['log_master_info'])
    return


class fakefile (object):
    def __init__(self):
        self.text = ""
    def write(self, t):
        self.text += t
    def get(self):
        return self.text

def print_stacktrace(sleep=0, logger=None, info=True, stdout=False):

    time.sleep(sleep)

    ff = fakefile()
    print >>ff, "========================================================"
    print >>ff, "==   STACK TRACE -- BEGIN                             =="
    print >>ff, "========================================================"

    print >>ff, "\nCurrently running threads:\n -- %s" % ("\n -- ".join([str(x) for x in threading.enumerate()]))

    for thread_id, frame in sys._current_frames().iteritems():
        name = thread_id
        #print name, frame

        for thread in threading.enumerate():
            if thread.ident == thread_id:
                name = thread.name
        print >>ff,"\nSTACK-TRACE for %s" % (name)
        traceback.print_stack(frame, file=ff)

    try:
        import psutil
        print >>ff, "\nList of subprocess of current process:"
        this_process = psutil.Process()
        kids = this_process.children(recursive=True)
        if (len(kids) <= 0):
            print >>ff, "  This process does not have any child-processes"
        else:
            now_time = time.time()
            print >>ff, " -- %s" % ("\n -- ".join(["ProcID: %5d - %8.3f seconds" % (
                p.pid, now_time-p.create_time()) for p in kids]))
    except:

        print >>ff, "\nList of subprocesses not available, needs psutil package!"
        pass

    print >>ff, ""
    print >>ff, "========================================================"
    print >>ff, "==   STACK TRACE -- END                               =="
    print >>ff, "========================================================"

    if (stdout):
        print ff.get()
    else:
        if (logger == None):
            logger = logging.getLogger("StackTrace")
        if (info):
            logger.info("\n%s" % (ff.get()))
        else:
            logger.debug("\n%s" % (ff.get()))

    time.sleep(sleep)
    
    return ff.get()


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

        try:
            import pika
            print "Found PIKA module!"
        except ImportError:
            print "There's no pika package installed, quitting"
            sys.exit(0)
        
        try:
            import podi_pikasetup
            print "Found podi-setup for PIKA!"
        except ImportError:
            print "No podi-setup found."
            print "Maybe you need to run `cp podi_pikasetup.py.example podi_pikasetup.py ?"
            sys.exit(0)

        print "Trying to connect to AMQP server"
        try:
            credentials = pika.PlainCredentials(podi_pikasetup.user, podi_pikasetup.password)
            connection = pika.BlockingConnection(
                pika.ConnectionParameters( 
                    credentials=credentials, 
                    host=podi_pikasetup.host, 
                    virtual_host=podi_pikasetup.vhost)
            )
            channel = connection.channel()
            channel.queue_declare(queue=podi_pikasetup.jobstatus_queue, durable=True)
        except:
            print "Connection failed!"
            sys.exit(0)

        print "PIKA connection established!"

        def callback(ch, method, properties, body):
            print " [o] %r" % (body,)

        channel.basic_consume(callback,
                              queue=podi_pikasetup.jobstatus_queue,
                              no_ack=True)

        print ' [*] Waiting for messages. To exit press CTRL+C'
        try:
            channel.start_consuming()
        except (KeyboardInterrupt, SystemExit):
            print "\nShutting down listener, good bye!"
        except:
            pass

        sys.exit(0)

    else:

        import podi_collectcells
        options = podi_collectcells.read_options_from_commandline()

        log_master_info, log_setup = podi_log_master_start(options)

        # Setup the multi-processing-safe logging
        podi_logger_setup(options['log_setup'])

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

        logger = logging.getLogger("main process")

        logger.info('test info')
        logger.debug('test debug')
        logger.critical('test critical')
        logger.error('test error')


        podi_log_master_quit(log_master_info)
    #    queue.put_nowait(None)
    #    listener.join()
