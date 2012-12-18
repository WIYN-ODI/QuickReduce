
#
# (c) Ralf Kotulla for WIYN/pODI
#



"""
Some general definitions useful in a number of podi scripts

"""

import sys
import os


available_otas = [00, 16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]
available_ota_coords = [(0,0), (1,6), 
                        (2,2), (2,3), (2,4),
                        (3,2), (3,3), (3,4),
                        (4,2), (4,3), (4,4),
                        (5,5), (6,1),
                        ]

broken_ota_cells = [ [(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7)],  #00
                     [],#1.6
                     [],#2.2
                     [(6,6)],#2,3
                     [],#2,4
                     [],#3,2
                     [],#3,3
                     [(1,7),(3,1),],#3,4
                     [],#4,2
                     [],#4,3
                     [],#4,4
                     [],#55
                     [(1,3),(3,1)],#61
                     ]


all_otas = [00, 16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]
central_2x2 = [22,23,32,33]

which_otas_to_use = {"odi_g": all_otas,
		     "odi_r": all_otas,
                     "odi_i": all_otas,
                     "odi_z": all_otas,
                     "CTIO_Ha": central_2x2,
	             "sdss_u": central_2x2,
	             "CTIO_OIII": central_2x2,
                     "BATC_420": central_2x2,
	             "BATC_390": central_2x2,
		    }	


def stdout_write(str):
    sys.stdout.write(str)
    sys.stdout.flush()
    return

def clobberfile(filename):
    if (os.path.isfile(filename)):
        os.remove(filename)
    return



from types import *   
def shmem_as_ndarray( raw_array ):
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
    address = raw_array._wrapper.get_address()
    size = raw_array._wrapper.get_size()
    dtype = _ctypes_to_numpy[raw_array._type_]
    class Dummy(object): pass
    d = Dummy()
    d.__array_interface__ = {
         'data' : (address, False),
         'typestr' : ">f4", #FloatType, #"uint8", #numpy.uint8.str,
         'descr' : "", #"UINT8", #numpy.uint8.descr,
         'shape' : (size/4,),
         'strides' : None,
         'version' : 3
    }
    return numpy.asarray(d)#.view( dtype=numpy.float32 )



def cmdline_arg_isset(arg):
    # Go through all command line arguments and check
    # if the requested argument is one of them
    for cur in sys.argv[1:]:
        name,sep,value = cur.partition("=")
        if (name == arg):
            return True
    return False


def get_cmdline_arg(arg):
    # Check all arguments if 
    for cur in sys.argv[1:]:
        name,sep,value = cur.partition("=")
        if (name == arg):
            return value
    return None


def get_clean_cmdline():
    list = []
    for cur in sys.argv:
        if (cur[0] != "-"):
            list.append(cur)
    return list


def cmdline_arg_set_or_default(name, defvalue):
    if (cmdline_arg_isset(name)):
        return get_cmdline_arg(name)
    return defvalue
