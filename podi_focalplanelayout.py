#!/usr/bin/env python

import pyfits
import logging
import itertools
import math, sys

import podi_crosstalk
import podi_logging
import numpy
import os

class FocalPlaneLayout(object):


    def __init__(self, inp=None):

        self.logger = logging.getLogger("FocalPlaneLayout")
        self.valid = False
        self.logger.debug("Creating focal plane from type %s" % (str(type(inp))))

        if (type(inp) == str):
            # Assume its a filename
            hdulist = pyfits.open(inp)
            hdu = hdulist[0]
        elif (type(inp) == pyfits.hdu.hdulist.HDUList):
            hdu = inp[0]
        elif (type(inp) == pyfits.hdu.image.ImageHDU or
              type(inp) == pyfits.hdu.compressed.CompImageHDU or
              type(inp) == pyfits.hdu.image.PrimaryHDU):
            hdu = inp
        elif (inp == None):
            # This is a fall-back mode, creating a class that can not do 
            # everything it could do otherwise
            return
        else:
            self.logger.error("Unrecognized type: %s" % (type(inp)))
            return

        # Assume this focal plane layout is well determined
        self.valid = True
        self.logger.debug("HDU type is now %s" % (str(type(hdu))))

        #
        # Set internal values to a safe value
        #
        self.wcs = None

        self.setup_general()
        
        # Find date of exposure
        mjd_obs = hdu.header['MJD-OBS']

        if (mjd_obs < 57023):
            # This exposure was obtained BEFORE 01/01/2015
            # --> use original pODI layout (3x3 + 4)
            self.setup_podi_layout()
            self.fp_config = 'podi'
        else:
            # This is the refurbished, extended layout commissioned in 2015
            self.setup_5x6_layout()
            self.fp_config = "odi_5x6"
        return


    #############################################################################
    #
    # The information below was formerly stored in the podi_definitions.py file
    #
    #############################################################################
    def setup_podi_layout(self):

        self.layout = 'pODI'

        self.wcs_default = "2mass_distort5.fits"

        # this is for pODI 
        self.available_ota_coords = [
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

        # this is for pODI 
        self.central_array_ota_coords = [
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

        # this is for pODI 
        self.blocked_out_otas = {
            "odi_g": [],
            "odi_i": [],
            "odi_r": [],
            "odi_z": [],
            #
            # All other filters
            #
              "823_v2": [55, ],
             "918R_v1": [55, ],
            "BATC_390": [55, ],
            "BATC_420": [55, ],
             "CTIO_Ha": [55, 61],
         "CTIO_Ha_8nm": [55, ],
           "CTIO_OIII": [55, ],
     "Halpha_and_odiz": [55, ],
            "KPNO_815": [55, ],
            "mosaic_u": [55, ],
    "MosaicU_and_odir": [55, ],
                "OPEN": [55, ],
              "s2_SII": [55, ],
              "sdss_u": [55, ],
                 "UG5": [55, ],
             "unknown": [55, ],
             "UNKNOWN": [55, ],
            "Us_solid": [55, ],
         "windowGlass": [55, ],
                "WRC3": [55, ],
            "Harris_B": [55, ],
            }
        

        # this is for pODI 
        self.broken_ota_cells = {
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

        # this is for pODI 
        self.all_otas = [00, 16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]
        self.central_3x3 = [22, 23,24, 32, 33, 34, 42, 43, 44]
        self.central_2x2 = [22,23,32,33]
        self.non_vignetted = [16, 22, 23, 24, 32, 33, 34, 42, 43, 44, 55, 61]

        #
        # Enter here the OTAs with full coverage in each of the filters. These OTAs
        # will then be used to compute the median illumination level to correct/normalize the flat-fields
        #
        # this is for pODI 
        self.otas_to_normalize_ff = {
            #
            # Full ODI filters
            #
            "odi_g": self.non_vignetted,
            "odi_i": self.non_vignetted,
            "odi_r": self.non_vignetted,
            "odi_z": self.non_vignetted,
            #
            # All other filters
            #
            "823_v2": self.central_2x2,
            "918R_v1": self.central_2x2,
            "BATC_390": self.central_2x2,
            "BATC_420": self.central_2x2,
            "CTIO_Ha": self.central_2x2,
            "CTIO_Ha_8nm": self.central_2x2,
            "CTIO_OIII": self.central_2x2,
            "Halpha_and_odiz": self.central_2x2,
            "KPNO_815": self.central_2x2,
            "mosaic_u": self.central_2x2,
            "MosaicU_and_odir": self.central_2x2,
            "OPEN": self.non_vignetted,
            "s2_SII": self.central_2x2,
            "CTIO_SII": self.central_2x2,
            "sdss_u": self.central_2x2,
            "UG5": self.central_2x2,
            "unknown": self.central_2x2,
            "UNKNOWN": self.central_2x2,
            "Us_solid": self.central_2x2,
            "windowGlass": self.all_otas,
            "WRC3": self.central_2x2,
            "Harris_B": self.central_2x2,
            "KPNO_815_v2": self.central_2x2,
            }	
        
        #
        # These are the OTAs that have at least partial coverage in each filter
        #
        # this is for pODI 
        self.otas_for_photometry = {
            #
            # Full ODI filters
            #
            "odi_g": self.all_otas,
            "odi_r": self.all_otas,
            "odi_i": self.all_otas,
            "odi_z": self.all_otas,
            #
            # All other filters
            #
            "823_v2": self.central_3x3,
            "918R_v1": self.central_3x3,
            "BATC_390": self.central_3x3,
            "BATC_420": self.central_3x3,
            "CTIO_Ha": self.central_3x3,
            "CTIO_Ha_8nm": self.central_3x3,
            "CTIO_OIII": self.central_3x3,
            "Halpha_and_odiz": self.central_3x3,
            "KPNO_815": self.central_3x3,
            "mosaic_u": self.central_3x3,
            "MosaicU_and_odir": self.central_3x3,
            "OPEN": self.all_otas,
            "sdss_u": self.central_3x3,
            "s2_SII": self.central_3x3,
            "CTIO_SII": self.central_3x3,
            "UG5": self.central_3x3,
            "unknown": self.central_3x3,
            "UNKNOWN": self.central_3x3,
            "Us_solid": self.central_3x3,
            "windowGlass": self.all_otas,
            "WRC3": self.central_3x3,
            "Harris_B": self.central_3x3,
            }	

        # this is for pODI 
        self.pupilghost_centers = {
            "odi_g": {"OTA33.SCI": (4182, 4155),
                      "OTA43.SCI": ( -23, 4147),
                      "OTA34.SCI": (4207, -189),
                      "OTA44.SCI": ( -12, -204)},

            "odi_i": {"OTA33.SCI": (4225, 4145),
                      "OTA34.SCI": (4225, -143),
                      "OTA44.SCI": ( -31, -207),
                      "OTA43.SCI": ( -63, 4081)},
            
            # use the direct coordinates as a backup for compatibility for now
            "OTA33.SCI": (4182, 4155),
            "OTA43.SCI": ( -23, 4147),
            "OTA34.SCI": (4207, -189),
            "OTA44.SCI": ( -12, -204),
            }

        pass

    def setup_5x6_layout(self):

        #
        # For now default to the podi layout
        # Change this once we know the final focal plane configuration
        #
        self.setup_podi_layout()

        self.layout = 'ODI_5x6'
        self.wcs_default = "wcs_5x6_skeleton.fits"

        self.available_ota_coords = itertools.product(range(1,7), repeat=2)# [
            # (3,3),
            # (3,4),
            # (4,4),
            # (4,3),
            # (4,2),
            # (3,2),
            # (2,2),
            # (2,3),
            # (2,4),
            # (5,5),
            # (6,1),
            # (1,6), 
            # (0,0),
            # ]

        self.central_array_ota_coords = self.available_ota_coords

        pass

    
    def setup_general(self):

        self.broken_cells = {
            #
            # Original pODI OTAs
            #--------------------------------------------------------------pODI ODI5
            '13838': [(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7)],  #  00   55
            '13880': [],                                                 #  16   25
            '13835': [],                                                 #  22   26
            '13901': [(6,6)],                                            #  23   14
            '13968': [],                                                 #  24   11
            '13974': [],                                                 #  32   12
            '13879': [],                                                 #  33   45
            '13923': [(1,7),(3,1),],                                     #  34   13
            '13792': [],                                                 #  42   46
            '13902': [],                                                 #  43   36
            '13947': [],                                                 #  44   35
            '13946': [],                                                 #  55   16
            '13837': [(1,3),(3,1)],                                      #  61   15
            #
            # Additions for ODI 5x6
            #
            '17198': [(4,0),
                      (0,7),(1,7),(2,7),(3,7),(4,7),(5,7),(6,7),(7,7),], #  --   21
            '17187': [(6,6)],                                            #  --   31
            '17234': [(0,0),
                      (0,7),(1,7),(2,7),(3,7),(4,7),(5,7),(6,7),(7,7),], #  --   41
            '17253': [(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7),
                      (6,0),(6,1),(6,2),(6,3),(6,4),(6,5),(6,6),(6,7),], #  --   51
            '17297': [],                                                 #  --   22
            '17231': [],                                                 #  --   32
            '17277': [],                                                 #  --   42
            '17190': [(1,5),
                      (7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6),(7,7),], #  --   52
            '17144': [(0,0),(0,7)],                                      #  --   23
            '17121': [],                                                 #  --   33
            '17341': [],                                                 #  --   43
            '17278': [],                                                 #  --   53
            '17275': [(0,1),(0,2),(1,1)],                                #  --   24
            '17166': [(3,0)],                                            #  --   34
            '17167': [],                                                 #  --   44
            '17122': [(5,0),(6,0),(7,0)],                                #  --   54
            '8101' : [],                                                 #  --   56
            #
            # These are added, but should not exist
            #
            '17189': [],
            }

        #
        # For testing, disable all broken cells
        #
        for ota in self.broken_cells:
            self.broken_cells[ota] = []


    def is_cell_broken(self, ota_id, cx, cy):
        cell_xy = (cx,cy)
        # self.logger.debug("OTA-ID: %s, cell %d,%d -> %s" % (ota_id, cx, cy, str(cell_xy in self.broken_cells[ota_id])))
        return (cell_xy in self.broken_cells[ota_id])

    def get_otaid_from_position(self, xy):
        if (self.fp_config == 'podi'):
            _fp = {00: 13838,
                   16: 13880,
                   22: 13835, 23: 13901, 24: 13968,
                   32: 13974, 33: 13879, 34: 13923,
                   42: 13792, 43: 13902, 44: 13947,
                   55: 13946,
                   61: 13837,
              }
        else:
            _fp = {11:13968, 12:13974, 13:13923, 14:13901, 15:13837, 16:13946,
                   21:17198, 22:17297, 23:17144, 24:17275, 25:13880, 26:13835, 
                   31:17187, 32:17231, 33:17121, 34:17166, 35:13947, 36:13902, 
                   41:17234, 42:17277, 43:17341, 44:17167, 45:13879, 46:13792,
                   51:17253, 52:17190, 53:17278, 54:17122, 55:13838, 56:8101, 
            }
        if (xy in _fp):
            return _fp[xy]
        return -1

    def get_crosstalk_matrix(self, extname):

        if (extname in podi_crosstalk.xtalk_matrix):
            return podi_crosstalk.xtalk_matrix[extname]

        return podi_crosstalk.xtalk_matrix['OTA33.SCI']

    def crosstalk_saturation_limit(self, extname):
        return podi_crosstalk.xtalk_saturation_limit

    def crosstalk_saturation_correction(self, extname):
        return podi_crosstalk.xtalk_saturated_correction


    def apply_wcs_distortion(self, filename, hdu, binning):

        if (os.path.isdir(filename)):
            # We received a directory instead of a filename
            # ==> use the default filename in this directory
            filename = "%s/%s" % (filename, self.wcs_default)
            self.logger.debug("Choosing default WCS solution --> %s" % (filename))

        try:
            self.wcs = pyfits.open(filename)
        except:
            self.logger.error("Could not open WCS distortion model (%s)" % (filename))
            return

        extname = hdu.header['EXTNAME']

        try:
            wcs_header = self.wcs[extname].header
        except:
            self.logger.warning("Could not find distortion model for %s" % (extname))
            return

        try:
            for hdr_name in ('CRPIX1', 'CRPIX2'):
                wcs_header[hdr_name] /= binning
            for hdr_name in ('CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'):
                wcs_header[hdr_name] *= binning

            cards = wcs_header.cards
            for (keyword, value, comment) in cards:
                if (keyword == 'CRVAL1'):
                    hdu.header["CRVAL1"] -= wcs_header['CRVAL1'] / math.cos(math.radians(hdu.header['CRVAL2']))
                elif (keyword == "CRVAL2"):
                    hdu.header['CRVAL2'] -= wcs_header['CRVAL2']
                else:
                    hdu.header[keyword] = (value, comment)

            # Make sure to write RAs that are positive
            if (hdu.header["CRVAL1"] < 0):
                hdu.header['CRVAL1'] += 360.
        except:
            self.logger.critical("something went wrong while applying the WCS model")
            podi_logging.log_exception()

            
        return


    def get_science_area_otas(self, filtername, include_vignetted=True):
        
        if (include_vignetted):
            if (filtername in self.otas_for_photometry):
                return self.otas_for_photometry[filtername]
            else:
                return self.central_2x2

        if (filtername in self.otas_to_normalize_ff):
            return self.otas_to_normalize_ff[filtername]
        
        return central_2x2




    def create_radially_sorted_ota_list(self):
        
        # all_xy = []
        
        # for x,y in itertools.product(range(8),repeat=2):
        
        #     center_x = x * 4300 + 2000
        #     center_y = y * 4300 + 2000

        #     r = math.sqrt((17000 - center_x)**2 + (17000-center_y)**2)
        
        #     all_xy.append([x,y,r,center_x, center_y])

        # all_xy = numpy.array(all_xy)
    
        # # sort by r
        # si = numpy.argsort(all_xy[:,2])

        # sortedxy = all_xy[si]

        # # numpy.savetxt(sys.stdout, sortedxy, "% 2d % 2d %.0f % 6d % 6d")
        # return list(numpy.array(sortedxy[:,0]*10+sortedxy[:,1], dtype=numpy.int))

        #
        # Compact version
        #
        y,x = numpy.indices((8,8))
        r = numpy.hypot(17000-(4300*x[0:1,:]+2000), 17000-(4300*y[:,0:1]+2000))
        si = numpy.argsort(r.ravel())
        idx = x[0:1,:]*10+y[:,0:1]
        return list(idx.ravel()[si])
