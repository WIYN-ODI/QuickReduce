
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
