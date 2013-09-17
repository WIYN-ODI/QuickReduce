#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#



"""
Some general definitions useful in a number of podi scripts

"""

import sys
import os
import numpy
import ctypes
import math
import numpy
import pyfits
import subprocess
import scipy
import scipy.ndimage
from bottleneck import nanmean, nanmedian


############
#                  | | |
#  Experimental    | | |
#                  V V V
############

#pyfits.USE_MEMMAP = False

############
#                  A A A
#  Experimental    | | |
#                  | | |
############

available_ota_coords = [
    (3,3),
    (3,4),
    (4,4),
    (4,3),
    (4,2),
    (3,2),
    (2,2),
    (2,3),
    (2,4),
    (5,5),
    (6,1),
    (1,6), 
    (0,0),
    ]

central_array_ota_coords = [
    (3,3),
    (3,4),
    (4,4),
    (4,3),
    (4,2),
    (3,2),
    (2,2),
    (2,3),
    (2,4),
    ]


broken_ota_cells = {
    "OTA00": [(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7)],  #00
    "OTA16": [],#1.6
    "OTA22": [],#2.2
    "OTA23": [(6,6)],#2,3
    "OTA24": [],#2,4
    "OTA32": [],#3,2
    "OTA33": [],#3,3
    "OTA34": [(1,7),(3,1),],#3,4
    "OTA42": [],#4,2
    "OTA43": [],#4,3
    "OTA44": [],#4,4
    "OTA55": [],#55
    "OTA61": [(1,3),(3,1)]
}

all_otas = [00, 16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]
central_3x3 = [22, 23,24, 32, 33, 34, 42, 43, 44]
central_2x2 = [22,23,32,33]


#
# Enter here the OTAs with full coverage in each of the filters. These OTAs
# will then be used to compute the median illumination level to correct/normalize the flat-fields
#
otas_to_normalize_ff = {"odi_g": [22,23,24,32,42],
                        "odi_r": all_otas,
                        "odi_i": all_otas,
                        "odi_z": all_otas,
                        "CTIO_Ha": central_2x2,
                        "sdss_u": central_2x2,
                        "CTIO_OIII": central_2x2,
                        "BATC_420": central_2x2,
                        "BATC_390": central_2x2,
                        "OPEN": all_otas,
                        "mosaic_u": central_2x2,
                        "s2_SII": central_2x2,
                        "823_v2": central_2x2,
                        "918R_v1": central_2x2,
                        "Us_solid": central_2x2,
                        "unknown": central_2x2,
                        }	

#
# These are the OTAs that have at least partial coverage in each filter
#
otas_for_photometry = {"odi_g": all_otas,
                       "odi_r": all_otas,
                       "odi_i": all_otas,
                       "odi_z": all_otas,
                       "CTIO_Ha": central_3x3,
                       "sdss_u": central_3x3,
                       "CTIO_OIII": central_3x3,
                       "BATC_420": central_3x3,
                       "BATC_390": central_3x3,
                       "OPEN": all_otas,
                       "mosaic_u": central_3x3,
                       "s2_SII": central_3x3,
                       "823_v2": central_3x3,
                       "918R_v1": central_3x3,
                       "unknown": central_2x2,
                       }	


#
# This list is used for photometric calibration.
# Enter the SDSS equivalent name for each of the filters
#
sdss_equivalents = {"odi_g": 'g',
                    "odi_r": 'r',
                    "odi_i": 'i',
                    "odi_z": 'z',
                    "CTIO_Ha": None,
                    "sdss_u": 'u',
                    "CTIO_OIII": None,
                    "BATC_420": None,
                    "BATC_390": None,
                    "OPEN": None,
                    "mosaic_u": 'u',
                    "s2_SII": None,
                    "823_v2": None,
                    "918R_v1": None,
                    "Us_solid": 'u',
                    "unknown": None,
		    }	

sdss_photometric_column = {"u":  2,
                           "g":  4,
                           "r":  6,
                           "i":  8,
                           "z": 10,
                           }

cellmode_ids = {
    "S": 0,
    "V": 1,
    # Add other cellmodes and some numerical representation here
}

list_of_valid_filter_ids = {
    "4": "odi_g",
    "w104": "odi_g",
    "w105": "odi_r",
    "w103": "odi_i",
    "w101": "odi_z",
    "c6009": "CTIO_Ha",
    "c6022": "sdss_u",
    "c6014": "CTIO_OIII",
    "k1053": "BATC_420",
    "k1052": "BATC_390",
    "clear": "OPEN",
#    "4": "mosaic_u",
#    "1": "s2_SII",
    "k1047": "823_v2",
    "k1028": "918R_v1",
#    "1": "Us_solid",
    }

def get_valid_filter_name(hdr):

    filter_id = hdr['FILTERID'].strip()
    if (filter_id in list_of_valid_filter_ids):
        filter_name = list_of_valid_filter_ids[filter_id]
        return filter_name

    return 'unknown'





def get_cellmode(primhdr, cellhdr):
    
    ota = int(primhdr['FPPOS'][2:4])

    ota_name = "OTA%02d" % ota
    extname = "OTA%02d.SCI" % ota

    cell = cellhdr['EXTNAME']

    # Check if this is one of the permanently broken cells
    wm_cellx, wm_celly = cellhdr['WN_CELLX'], cellhdr['WN_CELLY']
    broken = False
    list_of_broken_cells = broken_ota_cells[ota_name]
    for broken_cell in list_of_broken_cells:
        x,y = broken_cell
        if (wm_cellx == x and wm_celly == y):
            return -1
            break

    # It's not one of the broken cells, but it might still be a guide/video cell
    idx = wm_cellx + 8 * wm_celly
    cellmode = primhdr['CELLMODE']
    this_cellmode = cellmode[idx]
    cell_id = cellmode_ids[this_cellmode]

    return cell_id


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
        if (cur[0] != "-" or cur[1].isdigit()):
            list.append(cur)
    return list


def cmdline_arg_set_or_default(name, defvalue):
    if (cmdline_arg_isset(name)):
        return get_cmdline_arg(name)
    return defvalue



def sexa2deg(sexa):
    components = sexa.split(":")
    
    if (len(components) != 3):
        return 1e99
    
    deg = float(components[0])
    min = float(components[1])
    sec = float(components[2])
    
    return math.copysign(math.fabs(deg) + math.fabs(min/60.0) + math.fabs(sec/3600.0), deg)

def deg2sexa(deg, signed=False):

    unsigned = math.fabs(deg)

    degrees = math.floor(unsigned)
    rest = (unsigned - degrees) * 60.0

    minutes = math.floor(rest)
    rest = (rest - minutes) * 60.0

    seconds = math.floor(rest)

    num = [math.copysign(degrees, deg), minutes, seconds]

    if (signed):
        text = "%+03d:%02d:%04.1f" % (int(math.copysign(degrees, deg)), int(minutes), seconds)
    else:
        text = "%02d:%02d:%04.1f" % (int(math.copysign(degrees, deg)), int(minutes), seconds)

    return text, num



headers_to_inherit = [
    'RA', 'DEC', 'TARGRA', 'TARGDEC', 'TELRAOFF', 'TELDECOF', 
    'FILTER', 'FILTERID', 'FILTDSCR', 'EXPTIME',
    'OBSID', 'OBJECT', 'OBSTYPE',
    'WCSASTRM',
    'EXPMEAS',
    
    'ORIGIN', 'INSTRUME',
    'FILENAME', 
    'OBSLOGIN',
    'RADESYS', 'TIMESYS', 'LSTHDR',
    'OBSERVAT', 'TELESCOP',
    'OBSERVER', 'PROPOSER', 'PROPID', 'PROGID', 'TACID', 'PROPPERD',
    'DATE-OBS', 'TIME-OBS', 'MJD-OBS', 'DATE',
    'ZD', 'AIRMASS',
    'TELFOCUS',
    'TRACK',
    'ELMAP', 'AZMAP', 'ROTPORT',
    'FOLDPOS','OBSBLOCK',
    'ADCMODE', 'ADCANG1', 'ADCANG2', 'ADCJD',
    'ROTSTART', 'ROTEND', 'ROTOFF', 
    
    'TEMPSTAT', 'DEWAR', 'COOLHEAD', 'COLPLATE', 'FOCPLATE', 'DEWPRESS',
    'FLTARM1A', 'FLTARM1B', 'FLTARM1C',
    'FLTARM2A', 'FLTARM2B', 'FLTARM2C',
    'FLTARM3A', 'FLTARM3B', 'FLTARM3C',
    'SHUTDIR', 'SHUTOPEN', 'SHUTCLOS',
    'CONTROLR',
    'IMAGESWV',
    ]

headers_to_delete_from_otas = [
    'CELLGAP1', 'CELLGAP2',
    'CNAXIS1', 'CNAXIS2',
    'NAMPS', 'NEXTEND', 
    'PRESCAN1', 'PRESCAN2',
    'OVRSCAN1', 'OVRSCAN2',
    'IMNAXIS1', 'IMNAXIS2',
    'EXTEND'
    ]


# .13
pupilghost_centers = {"OTA33.SCI": (4190, 4150),
                      "OTA34.SCI": (4205, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }

for ext in pupilghost_centers:
    cx, cy = pupilghost_centers[ext]
    cx += -9
    cy += +4
    pupilghost_centers[ext] = (cx, cy)
#    pupilghost_centers[ext][0] += 1
#    pupilghost_centers[ext][1] += 1


# .14
pupilghost_centers = {"OTA33.SCI": (4181, 4154),
                      "OTA34.SCI": (4196, -161),
                      "OTA44.SCI": (   1, -181),
                      "OTA43.SCI": (   1, 4134),
                      }


# .15
pupilghost_centers = {"OTA33.SCI": (4181+2, 4154+2),
                      "OTA34.SCI": (4196+2, -161-2),
                      "OTA44.SCI": (   1-2, -181-2),
                      "OTA43.SCI": (   1, 4134),
                      }


# .16
pupilghost_centers = {"OTA33.SCI": (4181-8, 4154+2),
                      "OTA34.SCI": (4196-10, -161-2),
                      "OTA44.SCI": (   1+5, -181-7),
                      "OTA43.SCI": (   1, 4134),
                      }
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4186, -163),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   1, 4134),
                      }

# 17
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   1, 4134),
                      }

# 18
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4195, -172),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   1, 4134),
                      }

# 19
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   -4, 4139),
                      }

# 20
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   -2, 4137),
                      }

# 21
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   8, -186),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 22
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   4, -190),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 23
pupilghost_centers = {"OTA33.SCI": (4173, 4156),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   0, -194),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 24
pupilghost_centers = {"OTA33.SCI": (4175, 4158),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   6, -188),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 25
pupilghost_centers = {"OTA33.SCI": (4175, 4158),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -178),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 26
pupilghost_centers = {"OTA33.SCI": (4175, 4158),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  21, -173),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 27
pupilghost_centers = {"OTA33.SCI": (4170, 4153),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -178),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 28
pupilghost_centers = {"OTA33.SCI": (4170, 4153),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  21, -173),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 29
pupilghost_centers = {"OTA33.SCI": (4164, 4147),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  15, -179),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 30
pupilghost_centers = {"OTA33.SCI": (4167, 4150),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  18, -176),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 31
pupilghost_centers = {"OTA33.SCI": (4166, 4149),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  17, -177),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 32
pupilghost_centers = {"OTA33.SCI": (4166, 4149),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  18, -176),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 33
pupilghost_centers = {"OTA33.SCI": (4167, 4150),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  17, -177),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 34
pupilghost_centers = {"OTA33.SCI": (4167, 4150),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  12, -177),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 35
pupilghost_centers = {"OTA33.SCI": (4167, 4150),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  12, -172),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 36
pupilghost_centers = {"OTA33.SCI": (4170, 4153),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  15, -169),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 37
pupilghost_centers = {"OTA33.SCI": (4164, 4147),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   9, -175),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 38
pupilghost_centers = {"OTA33.SCI": (4165, 4148),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  10, -174),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 39
pupilghost_centers = {"OTA33.SCI": (4165, 4148),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (   9, -175),
                      "OTA43.SCI": (   -1, 4136),
                      }




# 40
pupilghost_centers = {"OTA33.SCI": (4164, 4147),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  10, -174),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 41
pupilghost_centers = {"OTA33.SCI": (4164, 4147),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  15, -169),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 42
pupilghost_centers = {"OTA33.SCI": (4167, 4150),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -168),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 43
pupilghost_centers = {"OTA33.SCI": (4161, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -168),
                      "OTA43.SCI": (   -1, 4136),
                      }

# 44
pupilghost_centers = {"OTA33.SCI": (4165, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -164),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 45
pupilghost_centers = {"OTA33.SCI": (4165, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  28, -152),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 46
pupilghost_centers = {"OTA33.SCI": (4165, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  20, -160),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 47
pupilghost_centers = {"OTA33.SCI": (4165, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  16, -164),
                      "OTA43.SCI": (  -1, 4136),
                      }



# 48
pupilghost_centers = {"OTA33.SCI": (4165, 4144),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  20, -160),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 49
pupilghost_centers = {"OTA33.SCI": (4165, 4140),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  20, -155),
                      "OTA43.SCI": (   -1, 4136),
                      }



# 50
pupilghost_centers = {"OTA33.SCI": (4165, 4140),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  10, -155),
                      "OTA43.SCI": (   -1, 4136),
                      }


# 51
pupilghost_centers = {"OTA33.SCI": (4165, 4140),
                      "OTA34.SCI": (4190, -167),
                      "OTA44.SCI": (  20, -145),
                      "OTA43.SCI": (   -1, 4136),
                      }



# .13
pupilghost_centers = {"OTA33.SCI": (4190, 4150),
                      "OTA34.SCI": (4205, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .52
pupilghost_centers = {"OTA33.SCI": (4185, 4150),
                      "OTA34.SCI": (4200, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .53
pupilghost_centers = {"OTA33.SCI": (4180, 4150),
                      "OTA34.SCI": (4195, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .54
pupilghost_centers = {"OTA33.SCI": (4175, 4150),
                      "OTA34.SCI": (4190, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .55
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4190, -165),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .56
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4195, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .57
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4185, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .58
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4175, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .59
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4175, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  15, 4130),
                      }


# .60
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4180, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .61
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4172, -173),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  13, 4133),
                      }



# .62
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4170, -175),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .63
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4165, -180),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }



# .64
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4185, -160),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .65
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4175, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  20, 4140),
                      }



# .66
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4165, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4130),
                      }


# .67
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4165, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4125), #0,-5
                      }

# .68
pupilghost_centers = {"OTA33.SCI": (4178, 4148), # +3 +3
                      "OTA34.SCI": (4165, -170), #
                      "OTA44.SCI": (   7, -188), # -3 -3
                      "OTA43.SCI": (  10, 4125), #
                      }

# .69 from 66
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4165, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4128), #0,-2
                      }

# .70 from 66
pupilghost_centers = {"OTA33.SCI": (4175, 4145),
                      "OTA34.SCI": (4165, -170),
                      "OTA44.SCI": (  10, -185),
                      "OTA43.SCI": (  10, 4133), #0,+3
                      }


# using swarp and d_pixel offsets
pupilghost_centers = {"OTA33.SCI": (4168, 4170),
                      "OTA34.SCI": (4209, -175),
                      "OTA44.SCI": (  22, -200),
                      "OTA43.SCI": (   7, 4144), #0,+3
                      }

pupilghost_centers = {"OTA33.SCI": (4158, 4170),
                     "OTA34.SCI": (4199, -175),
                     "OTA44.SCI": (  12, -200),
                     "OTA43.SCI": (  -3, 4144), #0,+3
                     }

#try_33.1
pupilghost_centers = {"OTA33.SCI": (4168, 4170),}

#try33.2
pupilghost_centers = {"OTA33.SCI": (4118, 4120),}
#try33.3
pupilghost_centers = {"OTA33.SCI": (4125, 4115),}
#try33.4
pupilghost_centers = {"OTA33.SCI": (4110, 4130),}
#try33.5
pupilghost_centers = {"OTA33.SCI": (4100, 4110),}
#try33.6
pupilghost_centers = {"OTA33.SCI": (4107, 4117),}
#try33.7
pupilghost_centers = {"OTA33.SCI": (4117, 4117),}
#try33.8
pupilghost_centers = {"OTA33.SCI": (4095, 4117),}
#try33.9
pupilghost_centers = {"OTA33.SCI": (4100, 4117),}
#try33.10
pupilghost_centers = {"OTA33.SCI": (4100, 4130),}
#try33.11
pupilghost_centers = {"OTA33.SCI": (4105, 4125),}
#try33.12
pupilghost_centers = {"OTA33.SCI": (4115, 4115),}
#try33.13
pupilghost_centers = {"OTA33.SCI": (4125, 4105),}
#try33.14
pupilghost_centers = {"OTA33.SCI": (4135, 4095),}
#try33.15/16/17
pupilghost_centers = {"OTA33.SCI": (4140, 4090),}
#try33.18
pupilghost_centers = {"OTA33.SCI": (4130, 4090),}
#try33.19
pupilghost_centers = {"OTA33.SCI": (4124, 4082),}

#try33.18 <- best fit for ota 33
pupilghost_centers = {"OTA33.SCI": (4130, 4090),}

# try34.01
pupilghost_centers = {"OTA34.SCI": (4199, -175)}
# try34.02
pupilghost_centers = {"OTA34.SCI": (4230, -225)}
# try34.02
pupilghost_centers = {"OTA34.SCI": (4230, -237)}


# try 43.01
pupilghost_centers = {"OTA43.SCI": (  -3, 4144)}
# try 43.02
pupilghost_centers = {"OTA43.SCI": (  17, 4127)}
# try 43.03
pupilghost_centers = {"OTA43.SCI": (   3, 4107)}
# try 43.04 <-- nailed it ;-)
pupilghost_centers = {"OTA43.SCI": (  -8, 4151)}


# try 44.01
pupilghost_centers = {"OTA44.SCI": (  12, -200)}
# try 44.02
pupilghost_centers = {"OTA44.SCI": (  12, -204)}


#try33.18 
pupilghost_centers = {"OTA33.SCI": (4130, 4090),}
#try33.20
pupilghost_centers = {"OTA33.SCI": (4176, 4182),}
#try33.21 <- good enough
pupilghost_centers = {"OTA33.SCI": (4172, 4182),}


# now all together
pupilghost_centers = {"OTA33.SCI": (4172, 4182),
                      "OTA43.SCI": (  -8, 4151),
                      "OTA34.SCI": (4230, -237),
                      "OTA44.SCI": (  12, -204)}

# try34.02
# This works well for creating a template from the g-band
pupilghost_centers = {"OTA33.SCI": (4182, 4155),
                      "OTA43.SCI": ( -23, 4147),
                      "OTA34.SCI": (4207, -189),
                      "OTA44.SCI": ( -12, -204)}

#

pupilghost_centers = {
    "odi_g": {"OTA33.SCI": (4182, 4155),
              "OTA43.SCI": ( -23, 4147),
              "OTA34.SCI": (4207, -189),
              "OTA44.SCI": ( -12, -204)},

    "odi_i": {"OTA33.SCI": (4064, 4148),
              "OTA34.SCI": (4084, -216),
              "OTA44.SCI": ( -84, -240),
              "OTA43.SCI": (-104, 4124)},

    # use the direct coordinates as a backup for compatibility for now
    "OTA33.SCI": (4182, 4155),
    "OTA43.SCI": ( -23, 4147),
    "OTA34.SCI": (4207, -189),
    "OTA44.SCI": ( -12, -204),
}
    


def inherit_headers(header, primary_header):
    for header in headers_to_inherit:
        if (not header in primary_header):
            print "Problem with header ",header
            continue

        card = primary_header.ascardlist()[header]
        header.update(card.key, card.value, card.comment)

        

def rebin_image(data, binfac):

    if (binfac < 1):
        stdout_write("Rebinning at the moment only supports binning to larger pixels with binfac>1\n")
        return None
    elif (binfac == 1):
        return data

    out_size_x, out_size_y = int(math.ceil(data.shape[0]*1.0/binfac)), int(math.ceil(data.shape[1]*1.0/binfac))

    if (out_size_x*binfac != data.shape[0] or out_size_y*binfac != data.shape[1]):
        # The input array size is not a multiple of the new binning
        # Create a slightly larger array to hold the data to be rebinned
        container = numpy.zeros(shape=(out_size_x*binfac, out_size_y*binfac))

        # And insert the original data
        container[0:data.shape[0], 0:data.shape[1]] = data[:,:]
    else:
        container = data 
        
#    rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)

    rb1 = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
    rb2 = nanmean(rb1, axis=-1)
    rebinned = nanmean(rb2, axis=1)

#    rb1 = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
#    rb2 = nanmedian(rb1, axis=-1)
#    rebinned = nanmedian(rb2, axis=1)

#    #.nanmean(axis=-1).nanmean(axis=1)

    return rebinned



def center_coords(hdr):
        
    try:
        centerx, centery = hdr['NAXIS1']/2, hdr['NAXIS2']/2
    except:
        centerx, centery = 2048., 2048.

    center_ra  = (centerx-hdr['CRPIX1'])*hdr['CD1_1'] + (centery-hdr['CRPIX2'])*hdr['CD1_2'] + hdr['CRVAL1']
    center_dec = (centerx-hdr['CRPIX1'])*hdr['CD2_1'] + (centery-hdr['CRPIX2'])*hdr['CD2_2'] + hdr['CRVAL2']

    return center_ra, center_dec

    

def break_region_string(str_region):
    reg = str_region[1:-1]
    x,dummy,y = reg.partition(",")
    x1,dummy,x2 = x.partition(":")
    y1,dummy,y2 = y.partition(":")
    return int(x1)-1, int(x2)-1, int(y1)-1, int(y2)-1

def extract_region(data, str_region):
    x1,x2,y1,y2 = break_region_string(str_region)
    return data[y1:y2+1, x1:x2+1]


def insert_into_array(data, from_region, target, target_region):

    fx1, fx2, fy1, fy2 = break_region_string(from_region)
    tx1, tx2, ty1, ty2 = break_region_string(target_region)

    if (fx2-fx1 != tx2-tx1 or fy2-fy1 != ty2-ty1):
        print "Dimensions do not match, doing nothing"
    else:
        target[ty1:ty2+1, tx1:tx2+1] = data[fy1:fy2+1, fx1:fx2+1]

    return 0

def mask_broken_regions(datablock, regionfile, verbose=False):

    counter = 0
    file = open(regionfile)
    for line in file:
        if (line[0:3] == "box"):
            coords = line[4:-2]
            coord_list = coords.split(",")
                        
            if (not datablock == None):
                x, y = int(float(coord_list[0])), int(float(coord_list[1]))
                dx, dy = int(0.5*float(coord_list[2])), int(0.5*float(coord_list[3]))
                #mask[y-dy:y+dy,x-dx:x+dx] = 1

                x1 = numpy.max([0, x-dx])
                x2 = numpy.min([datablock.shape[1], x+dx])
                y1 = numpy.max([0, y-dy])
                y2 = numpy.min([datablock.shape[0], y+dy])
                datablock[y1:y2, x1:x2] = numpy.NaN

                # print x,x+dx,y,y+dy
            counter += 1

    file.close()
    if (verbose):
        print "Marked",counter,"bad pixel regions"
    return datablock



def is_image_extension(hdr):
    try:
        extname = hdr['EXTNAME']
        if (extname[:3] == "OTA" and extname[-4:] == ".SCI"):
            return True
    except:
        pass

    return False



def get_svn_version():

    try:
        p = subprocess.Popen('svnversion -n', shell=True, stdout=subprocess.PIPE)
        svn_version, err = p.communicate()
        ret = p.wait()
        if (ret != 0):
            svn_version = "problem_with_svn"
    except:
        svn_version="no_svnversion_found"

    return svn_version

def log_svn_version(hdr):
    svn = get_svn_version()
    hdr.update("QPIPESVN", svn, "QuickReduce Revision")
    return



def rotate_around_center(data, angle, mask_limit = 0.1, verbose=True, safety=1, mask_nans=True, spline_order=3):

    if (verbose):
        stdout_write("Rotating data block by %.1f deg ..." % (angle))

    # Prepare mask so we can mask out undefined regions
    mask = numpy.zeros(shape=data.shape)
    mask[numpy.isnan(data)] = 1.

    # Replace NaN's with some numer to make interpolation work
    data[numpy.isnan(data)] = 0.


    # Now rotate both the image and its mask
    rotated = scipy.ndimage.interpolation.rotate(input=data, angle=angle, axes=(1, 0), 
                                                 reshape=False, order=spline_order,
                                                 mode='constant', cval=0, )

    if (mask_nans):
        rotated_mask = scipy.ndimage.interpolation.rotate(input=mask, angle=angle, axes=(1, 0), 
                                                          reshape=False, order=1,
                                                          mode='constant', cval=1, )

        # Blur out the NaN mask to make sure we get clean edges. This approach
        # is on the conservative side, rather clipping some pixels too many then to 
        # add artificial edges that later are hard to remove and/or deal with
        filter_gaussed = scipy.ndimage.filters.gaussian_filter(input=rotated_mask, order=0, sigma=safety)

        # And finally apply the mask
        # rotated[rotated_mask > mask_limit] = numpy.NaN
        rotated[filter_gaussed > mask_limit] = numpy.NaN
        
    # and return the results
    if (verbose): stdout_write(" done!\n")
    return rotated



def get_filter_level(header):
    
    filter_level = 0
    filter_count = 0
 
    for lvl in range(1,4):
        for arm in "ABC":
            keyword = "FLTARM%d%s" % (lvl, arm)
            # print keyword

            if (header[keyword].strip() == "IN"):
                filter_count += 1
                if (filter_count == 1):
                    filter_level = lvl

    if (filter_count > 1):
        return -1

    return filter_level
            




def cell2ota__extract_data_from_cell(data_in=None):

    if (data_in == None):
        return numpy.zeros(shape=(494,480))

    return data_in[0:494, 0:480]

def cell2ota__get_target_region(x, y):

    #taken from ODI OTA Technical Datasheet (det area 480x494, streets 11/28 px)
    _y = 7 - y
    y1 = (505*_y)  #was 503 
    y2 = y1 + 494
    x1 = 508 * x
    x2 = x1 + 480

    return x1, x2, y1, y2




def three_sigma_clip(input, ranges=[-1e9,1e9], nrep=3, return_mask=False):

    valid = (input > ranges[0]) & (input < ranges[1])
                
    for rep in range(nrep):
        lsig = scipy.stats.scoreatpercentile(input[valid], 16)
        hsig = scipy.stats.scoreatpercentile(input[valid], 84)
        median = numpy.median(input[valid])
        sigma = 0.5 * (hsig - lsig)

        mingood = numpy.max([median - 3*sigma, ranges[0]])
        maxgood = numpy.min([median + 3*sigma, ranges[1]])

            #print median, sigma
        valid = (input > mingood) & (input < maxgood)
        
    if (return_mask):
        return input[valid], valid

    return input[valid]


def derive_ota_outlines(otalist):
    from astLib import astWCS

    all_corners = []
    for ext in range(len(otalist)):
        if (type(otalist[ext]) == pyfits.hdu.image.ImageHDU):
            
            wcs = astWCS.WCS(otalist[ext].header, mode='pyfits')
            
            corner_coords = []
            corner_coords.append(wcs.pix2wcs(                            0,                             0))
            corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'],                             0))
            corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'], otalist[ext].header['NAXIS2']))
            corner_coords.append(wcs.pix2wcs(                            0, otalist[ext].header['NAXIS2']))

            all_corners.append(corner_coords)

    return all_corners
