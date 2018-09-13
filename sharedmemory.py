#!/usr/bin/env python

# from podi_definitions import *
import multiprocessing
import time
import ctypes
import sys

from types import *   


try:
    import numpy
    _ctypes_to_numpy = {
        ctypes.c_char : numpy.int8,
        ctypes.c_wchar : numpy.int16,
        ctypes.c_byte : numpy.int8,
        ctypes.c_ubyte : numpy.uint8,
        ctypes.c_short : numpy.int16,
        ctypes.c_ushort : numpy.uint16,
        ctypes.c_int : numpy.int32,
        ctypes.c_uint : numpy.int32,
        ctypes.c_long : numpy.int32,
        ctypes.c_ulong : numpy.int32,
        ctypes.c_float : numpy.float32,
        ctypes.c_double : numpy.float64
    }
    numpy_available = True
except:
    numpy_available = False

def rebuild_after_pickle(obj):
    obj.rebuild()
    return obj

import ctypes
import multiprocessing.heap
#from multiprocessing.forking import assert_spawning, ForkingPickler
import mmap
import _multiprocessing
import logging

def address_of_buffer(buf):
    return ctypes.addressof(ctypes.c_char.from_buffer(buf))

class SharedMemory( object ):

    def __init__(self, _type, shape, log=False):
        # print "\n"*10, self

        self.logger = logging.getLogger("SharedMemory")
        self.log = log

        self._shape = shape
        self._type = _type
        self.allocated = False

        n_vars = 1
        for d in shape:
            n_vars *= d

        
        #
        # Recreate some shorter version of the workflow in multiprocessing.RawArray()
        #


        # typecode_or_type = _type
        # type_ = typecode_to_type.get(typecode_or_type, typecode_or_type)

        self.type_and_size = _type * n_vars
        self.unit_size = ctypes.sizeof(_type)
        self.n_bytes = self.unit_size * n_vars
        # create mmap object to allocate memory
        if (self.log):
            self.logger.debug("Allocating %.1f MBytes shared memory" % (self.n_bytes/1024/1024.))
        self.mm = mmap.mmap(-1, self.n_bytes)
        self.allocated = True

        # now convert the memory-map into the array
        self.rebuild()

    def rebuild(self):
        #ForkingPickler.register(self.type_and_size, rebuild_after_pickle)
        #(mm_address, mm_size) = _multiprocessing.address_of_buffer(self.mm)
        #(mm_address, mm_size) = address_of_buffer(self.mm)
        mm_address = address_of_buffer(self.mm)
        self._array = self.type_and_size.from_address(mm_address)
        
# # obj = _new_value(type_)
# size = ctypes.sizeof(type_)
# wrapper = heap.BufferWrapper(size)

# def rebuild_ctype(type_, wrapper, length):
#     if length is not None:
#         type_ = type_ * length
#     ForkingPickler.register(type_, reduce_ctype)
#     obj = type_.from_address(wrapper.get_address())
#     obj._wrapper = wrapper
#     return obj


# obj = rebuild_ctype(type_, wrapper, None)


        # ctypes.memset(ctypes.addressof(obj), 0, ctypes.sizeof(obj))

        # self._heap = multiprocessing.heap.Heap()

        # print "allocating %d xxx" % (n_vars)
        # self.shmem = multiprocessing.RawArray(self._type, n_vars)
        # print "memory allocated"

        self.allocated = True

    def to_ndarray(self):

        if (not numpy_available):
            return None


        """

        Helper function needed to allocate shared memory numpy-arrays.

        """
        mm_address = address_of_buffer(self.mm)
        #(mm_address, mm_size) = _multiprocessing.address_of_buffer(self.mm)
        #print("mm-size:", mm_size, self.n_bytes, self.n_bytes/4)
        #print("mm-address", mm_address, address_of_buffer(self.mm))

        #address = self.shmem._wrapper.get_address()
        #size = self.shmem._wrapper.get_size()
        dtype = _ctypes_to_numpy[self._type]
        class Dummy(object): pass
        d = Dummy()
        d.__array_interface__ = {
             'data' : (mm_address, False),
             'typestr' : ">f4", #FloatType, #"uint8", #numpy.uint8.str,
             'descr' : "", #"UINT8", #numpy.uint8.descr,
             'shape' : (self.n_bytes//4,), #(self.n_bytes/4), #
             'strides' : None,
             'version' : 3
        }
        return numpy.asarray(d).reshape(self._shape)#.view( dtype=numpy.float32 )

    def free(self):

        if (not self.allocated):
            return

        #
        # Now free up the shared memory in a way that actually releases the memory
        #
        if (self.log):
            self.logger.debug("Freeing up %.1f MBytes shared memory" % (self.n_bytes/1024/1024.))
        #print "freeing up memory"

        #blocks = shmem._wrapper._heap._allocated_blocks
        #for block in blocks:
        #    (arena, start, stop) = block

        # block,size = self.shmem._wrapper._state

        # del self.shmem._wrapper

        # print block
        # print size

        # # try:
        # #     self.shmem._wrapper._heap.free(block)
        # # except KeyError:
        # #     print "HEAP should be free already!"

        # (arena, start, stop) = block
        # # real memory is in arena.buffer
        # print arena

        # try:
        #     arena.buffer.close()
        # except:
        #     print "This should have succeeded"
        # print "memory should be freed now"

        self.mm.close()

        self.allocated = False

    def is_allocated(self):
        return self.allocated


    # def __del__(self):
    #     print "running __del__"
    #     self.free()



if __name__ == "__main__":

    shape = (1024, 1024, 1024/4, int(sys.argv[1]))

    

    time.sleep(2)
    print("Allocating shared memory")
    shmem = SharedMemory(ctypes.c_float, shape)
    print(shmem)

    time.sleep(2)


    #
    # Start new process, try to access the memory there
    #
    def mp(shmem):
        nd = shmem.to_ndarray()
        print(nd.shape)

    print("starting subprocess")
    p = multiprocessing.Process(target=mp, args=(shmem,))
    p.start()

    time.sleep(2)

    print("deleting shmam")
    #shmem.free()
    del shmem
    time.sleep(2)

    def run():
        print("running run")
        x = SharedMemory(ctypes.c_float, (1024,1024,1024))
        print("done!")

    run()
    time.sleep(2)

    # print "converting to numpy!"

    # nd = shmem.to_ndarray()
    # print nd.shape


    # time.sleep(3)


    # print "freeing memory"
    # del shmem
    # print "freed!"

    # time.sleep(3)

    # shmem2 = SharedMemory(ctypes.c_float, shape)
    # del shmem2

    # print "done!"



