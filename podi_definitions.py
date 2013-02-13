#
# (c) Ralf Kotulla for WIYN/pODI
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
otas_to_normalize_ff = {"odi_g": all_otas,
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


pupilghost_centers = {
    "odi_r": {"OTA33.SCI": (4080, 4180),
              "OTA34.SCI": (4115, -190),
              "OTA44.SCI": (-90, -230),
              "OTA43.SCI": (-120, 4150),
              },
    "odi_i": {"OTA33.SCI": (4080, 4180),
              "OTA34.SCI": (4115, -190),
              "OTA44.SCI": (-90, -230),
              "OTA43.SCI": (-120, 4150),
              }
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
        
    rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)

    return rebinned



def center_coords(hdr):
        
    centerx, centery = hdr['NAXIS1']/2, hdr['NAXIS2']/2
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

