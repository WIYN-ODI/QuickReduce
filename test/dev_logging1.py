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



if __name__ == "__main__":

    print "starting program"

    x = podi_logger('test')
    x.info("podi-logging works")

    logging.debug('some debug test')

    logging.info('some info test')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    #console = logging.StreamHandler()
    #console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    #console.setFormatter(formatter)
    # add the handler to the root logger
    #logging.getLogger('').addHandler(console)

    progress = ProgressConsoleHandler(stream=sys.stdout)
    progress.setFormatter(formatter)

    use_pika = False
    if (use_pika):
        pika = PikaMessageHandler(queue='hello')
        pika.setFormatter(formatter)
        logging.getLogger('').addHandler(pika)

#    console  = logging.StreamHandler()  
#    logger = logging.getLogger('test')
#    logger.setLevel(logging.DEBUG) 
#    logging.getLogger('').addHandler(progress)
    
    log = multiprocessing.get_logger()
    log.addHandler(progress)
    log.info("test")


    # processes = []
    # for i in range(4):
    #     worker_args = (i,i)
    #     p = multiprocessing.Process(target=dump_to_log, args=worker_args)
    #     p.start()
    #     processes.append(p)


    # Now, we can log to the root logger, or any other logger. First the root...
    logging.info('Jackdaws love my big sphinx of quartz.')

    # Now, define a couple of other loggers which might represent areas in your
    # application:

    logger1 = logging.getLogger('myapp.area1')
    logger2 = logging.getLogger('myapp.area2')

    logger1.debug('Quick zephyrs blow, vexing daft Jim.')
    logger1.info('How quickly daft jumping zebras vex.')
    logger2.warning('Jail zesty vixen who grabbed pay from quack.')
    logger2.error('The five boxing wizards jump quickly.')

    #global x
    #del x

    sys.exit(0)

