#!/usr/bin/env python

import pyfits
import logging
import itertools

import podi_crosstalk


class FocalPlaneLayout(object):

    def __init__(self, inp):

        self.logger = logging.getLogger("FocalPlaneLayout")
        self.valid = False

        if (type(inp) == str):
            # Assume its a filename
            hdulist = pyfits.open(inp)
            hdu = hdulist[0]
        elif (type(inp) == pyfits.hdu.hdulist.HDUList):
            hdu = hdulist[0]
        elif (type(inp) == pyfits.hdu.image.ImageHDU or
              type(inp) == pyfits.hdu.compressed.CompImageHDU or
              type(inp) == pyfits.hdu.image.PrimaryHDU):
            hdu = inp
        else:
            self.logger.error("Unrecognized type: %s" % (type(inp)))
            return

        # Assume this focal plane layout is well determined
        self.valid = True

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
            }

    def is_cell_broken(self, ota_id, cx, cy):
        cell_xy = (cx,cy)
        return (cell_xy in self.broken_cells[ota_id])



    def get_crosstalk_matrix(self, extname):

        if (extname in podi_crosstalk.xtalk_matrix):
            return podi_crosstalk.xtalk_matrix[extname]

        return podi_crosstalk.xtalk_matrix['OTA33.SCI']

    def crosstalk_saturation_limit(self, extname):
        return podi_crosstalk.xtalk_saturation_limit

    def crosstalk_saturation_correction(self, extname):
        return podi_crosstalk.xtalk_saturated_correction
