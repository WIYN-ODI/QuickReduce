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

This module contains a number of global constants required during reduction, as
well as functions that are widely used throughout the other modules.  Some
general definitions useful in a number of podi scripts

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
import scipy.special
from bottleneck import nanmean, nanmedian

from podi_commandline import *
import podi_sitesetup as sitesetup

#
# Update this variable below when changing versions
#
pipeline_plver = "QuickReduce 1.1"
pipeline_name = "QuickReduce"
pipeline_version = "1.1"
try:
    cmd = "svnversion -n %s" % (sitesetup.exec_dir)
    ret = subprocess.Popen(cmd.split(), 
                           stdout=subprocess.PIPE, 
                           stderr=subprocess.PIPE)
    out, err = ret.communicate()
    if (ret.returncode == 0 and not out.startswith("Un")):
        pipeline_version = "%s.%s" % (pipeline_version, out)
except:
    pass

#
#
#


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



#
# This list is used for photometric calibration.
# Enter the SDSS equivalent name for each of the filters
#
sdss_equivalents = {
    #
    # Full ODI filters
    #
    "odi_g": 'g',
    "odi_r": 'r',
    "odi_i": 'i',
    "odi_z": 'z',
    #
    # All other filters
    #
    "823_v2": "z",
    "918R_v1": "z",
    "BATC_390": "u",
    "BATC_420": "g",
    "CTIO_Ha": "r",
    "CTIO_Ha_8nm": "r",
    "CTIO_OII": "u",
    "CTIO_OIII": "g",
    "CTIO_SII": "r",   # same filter as s2_SII below
    "Halpha_and_odiz": None,
    "KPNO_815": "z",
    "KPNO_815_v2": "z",
    "mosaic_u": 'u',
    "MosaicU_and_odir": None,
    "OPEN": "r",
    "sdss_u": 'u',
    "s2_SII": "r",
    "UG5": None,
    "unknown": None,
    "UNKNOWN": None,
    "Us_solid": 'u',
    "windowGlass": None,
    "WRC3": "g",
    "Harris_B": "g",
    }	


filter_bandpass = {
    #           filter-name:      filter definition filename mean_pos, center,    fwhm,  left&right fwhm,     max,  mean,    area, left-5%, right-5%

    #
    # ODI filters
    #
                    "odi_g": (  		 "odi_g.txt", 4816.17, 4739.85,  1516.4,  4032.6,  5549.0, 0.3199, 0.2209,  412.56, 3989.12, 5596.64),
                    "odi_r": (  		 "odi_r.txt", 6231.93, 6214.50,  1383.1,  5526.7,  6909.8, 0.3910, 0.3109,  516.96, 5469.03, 6965.25),
                    "odi_i": (  		 "odi_i.txt", 7524.93, 7560.00,  1268.8,  6919.1,  8188.0, 0.3806, 0.2512,  409.39, 6845.65, 8277.06),
                    "odi_z": (  		 "odi_z.txt", 8605.45, 8643.00,   888.2,  8168.2,  9056.4, 0.2494, 0.1442,  197.98, 8069.76, 9261.17),
    #
    # mosaic filters - some may need to be renamed to match the name in the ODI filter system
    #
                  "CTIO_Ha": (         "c6009_ha_Halpha.txt", 6562.13, 6562.00,    80.1,  6522.4,  6602.5, 0.3035, 0.1199,   25.07, 6496.66, 6629.17),
                "CTIO_OIII": (  	 "c6014_o3_OIII.txt", 4999.06, 5001.80,    50.1,  4973.3,  5023.4, 0.2097, 0.0839,   10.90, 4954.99, 5046.71),
                   "sdss_u": (  	  "c6022_u_SDSS.txt", 3609.31, 3592.00,   491.3,  3363.0,  3854.4, 0.1639, 0.0861,   76.34, 3215.61, 3971.17),
                 "BATC_390": (  	"k1052_390_BATC.txt", 3903.53, 3926.50,   260.8,  3764.3,  4025.2, 0.1800, 0.0920,   44.01, 3720.52, 4091.28),
                 "BATC_420": (  	"k1053_420_BATC.txt", 4175.67, 4190.50,   324.2,  4009.4,  4333.6, 0.1960, 0.1217,   60.64, 3969.61, 4392.18),
                  "918R_v1": (      "k1028_918R_918R_v1.txt", 9203.83, 9170.21,    79.4,  9171.1,  9250.5, 0.1028, 0.0365,    8.70, 9140.66, 9277.68),
                   "823_v2": (        "k1047_823_823_v2.txt", 8238.68, 8239.00,    68.3,  8204.2,  8272.6, 0.2404, 0.1106,   16.40, 8187.34, 8289.58),
                     "WRC3": (      "wrc3_WR_CIII_k1021.txt", 4660.27, 4661.50,    51.6,  4634.5,  4686.1, 0.1996, 0.0923,   10.54, 4620.85, 4700.38),
                   "s2_SII": (  	  "c6013_s2_SII.txt", 6722.04, 6701.65,    80.5,  6679.5,  6760.0, 0.2962, 0.1131,   24.67, 6653.53, 6787.84),
                 "CTIO_SII": (  	  "c6013_s2_SII.txt", 6722.04, 6701.65,    80.5,  6679.5,  6760.0, 0.2962, 0.1131,   24.67, 6653.53, 6787.84),
                 "KPNO_815": (        "k1026_815_815_v1.txt", 8151.56, 8152.00,    70.5,  8116.5,  8187.1, 0.2561, 0.1163,   18.24, 8098.74, 8206.36),
              "CTIO_Ha_8nm": (    "c6011_ha8_Halpha+8nm.txt", 6652.38, 6652.25,    73.8,  6615.7,  6689.6, 0.3457, 0.1267,   26.28, 6589.44, 6716.45),
                 "Us_solid": (        "k1044_Us_solid_U.txt", 3698.69, 3642.00,   487.6,  3468.4,  3956.1, 0.2467, 0.1271,  114.21, 3264.72, 4015.57),
                 "Harris_B": (  	"c6002_B_Harris.txt", 4358.13, 4650.00,  1014.7,  3814.3,  4829.1, 0.1890, 0.0849,  191.36, 3605.61, 5211.97),
                 "CTIO_OII": (  	  "c6012_o2_OII.txt", 3726.11, 3727.25,    49.7,  3700.7,  3750.4, 0.1670, 0.0730,    8.57, 3682.05, 3771.56),


    #
    # more mosaic filters - mostly for completeness - many of which have not yet been used in ODI
    #
                   "quotaU": (  		"quotaU.dat", 3691.79, 3630.00,   482.8,  3455.6,  3938.4, 0.2443, 0.1257,  112.13, 3253.70, 3997.21),
                  "c6001_U": (  	       "c6001_U.txt", 3655.95, 3656.02,   528.5,  3388.9,  3917.4, 0.1910, 0.0980,   97.41, 3225.31, 4057.35),
           "c6004_R_Harris": (  	"c6004_R_Harris.txt", 6575.37, 7377.36,  1450.9,  5695.2,  7146.2, 0.3203, 0.1275,  474.31, 5588.92, 8331.47),
       "c6007_M_Washington": (      "c6007_M_Washington.txt", 5350.75, 5665.93,  1449.9,  4565.3,  6015.3, 0.2832, 0.1641,  393.86, 4500.81, 6526.56),
         "c6008_D51_DDO_51": (        "c6008_D51_DDO_51.txt", 5140.39, 5137.00,   164.3,  5054.3,  5218.6, 0.2490, 0.1113,   41.65, 5010.35, 5266.26),
             "c6017_g_SDSS": (  	  "c6017_g_SDSS.txt", 4802.37, 4523.57,  1405.4,  4081.2,  5486.7, 0.2779, 0.1711,  344.60, 3911.93, 5529.97),
             "c6018_r_SDSS": (  	  "c6018_r_SDSS.txt", 6256.66, 6234.01,  1388.9,  5564.8,  6953.8, 0.3715, 0.2844,  470.17, 5464.22, 6993.57),
             "c6019_i_SDSS": (  	  "c6019_i_SDSS.txt", 7715.44, 7684.00,  1419.1,  7052.1,  8471.3, 0.3127, 0.2036,  375.85, 6859.41, 8517.47),
             "c6020_z_SDSS": (  	  "c6020_z_SDSS.txt", 8902.49, 8927.00,   981.2,  8296.7,  9277.9, 0.1901, 0.0924,  197.50, 8044.05, 9998.00),
           "c6024_Bj_Tyson": (  	"c6024_Bj_Tyson.txt", 4406.86, 4627.38,  1516.8,  3661.7,  5178.6, 0.2432, 0.1314,  356.75, 3285.31, 5313.76),
         "c6025_It_Tyson_I": (        "c6025_It_Tyson_I.txt", 8518.51, 8636.37,  1351.6,  7705.3,  9056.9, 0.2549, 0.1323,  359.30, 7434.68, 9996.90),
           "c6026_V_Harris": (  	"c6026_V_Harris.txt", 5446.90, 5682.09,   991.8,  4906.7,  5898.5, 0.2742, 0.1401,  270.14, 4797.11, 6351.86),
                  "c6028_I": (  	       "c6028_I.txt", 7947.11, 8012.00,  1410.4,  7279.0,  8689.5, 0.3198, 0.1985,  385.47, 7178.43, 8805.38),
                  "k1001_U": (  	       "k1001_U.txt", 3648.13, 3639.00,   515.8,  3389.6,  3905.5, 0.2011, 0.1030,   99.40, 3225.02, 4041.07),
           "k1002_B_Harris": (  	"k1002_B_Harris.txt", 4389.83, 4709.00,   968.4,  3863.6,  4832.1, 0.2102, 0.0932,  206.61, 3637.00, 5242.87),
           "k1003_V_Harris": (  	"k1003_V_Harris.txt", 5459.44, 5734.00,  1003.4,  4910.5,  5914.0, 0.2738, 0.1452,  271.54, 4834.02, 6366.58),
           "k1004_R_Harris": (  	"k1004_R_Harris.txt", 6589.38, 7396.00,  1454.6,  5696.0,  7150.6, 0.3370, 0.1328,  498.01, 5592.21, 8606.07),
     "k1005_I_Nearly_Mould": (    "k1005_I_Nearly_Mould.txt", 8060.02, 8165.00,  1564.2,  7251.7,  8816.0, 0.3087, 0.1867,  436.96, 7078.62, 9249.18),
       "k1006_C_Washington": (      "k1006_C_Washington.txt", 3962.75, 4040.00,   923.8,  3479.0,  4402.9, 0.2282, 0.1244,  202.84, 3280.91, 4730.99),
       "k1007_M_Washington": (      "k1007_M_Washington.txt", 5247.02, 5462.01,  1193.1,  4607.1,  5800.3, 0.2649, 0.1523,  303.72, 4502.71, 6205.88),
         "k1008_D51_DDO_51": (        "k1008_D51_DDO_51.txt", 5145.55, 5147.00,   162.7,  5064.2,  5227.0, 0.2611, 0.1569,   41.45, 5040.83, 5250.47),
          "k1009_ha_Halpha": (         "k1009_ha_Halpha.txt", 6574.83, 6575.00,    80.5,  6534.4,  6615.0, 0.3781, 0.1567,   31.37, 6511.84, 6638.28),
     "k1010_ha4_Halpha+4nm": (    "k1010_ha4_Halpha+4nm.txt", 6620.84, 6623.00,    80.4,  6580.3,  6660.8, 0.3680, 0.1539,   30.66, 6559.02, 6683.46),
     "k1011_ha8_Halpha+8nm": (    "k1011_ha8_Halpha+8nm.txt", 6654.37, 6657.00,    81.2,  6613.5,  6694.8, 0.3629, 0.1516,   30.51, 6592.63, 6717.57),
   "k1012_ha12_Halpha+12nm": (  "k1012_ha12_Halpha+12nm.txt", 6698.75, 6700.00,    82.8,  6657.1,  6740.0, 0.3458, 0.1451,   29.78, 6635.13, 6762.87),
   "k1013_ha16_Halpha+16nm": (  "k1013_ha16_Halpha+16nm.txt", 6730.77, 6731.00,    81.0,  6690.2,  6771.3, 0.3645, 0.1519,   30.57, 6667.94, 6793.68),
         "k1014_O3_OIII_N2": (        "k1014_O3_OIII_N2.txt", 5025.16, 5026.00,    55.4,  4997.2,  5052.7, 0.2200, 0.1030,   12.48, 4984.23, 5066.47),
     "k1015_Ooff_OIII+30nm": (    "k1015_Ooff_OIII+30nm.txt", 5317.18, 5316.00,   242.1,  5197.7,  5439.9, 0.2960, 0.1560,   68.21, 5154.73, 5475.68),
             "k1016_wh_Bk7": (  	  "k1016_wh_Bk7.txt", 6184.92, 6182.45,  5162.9,  3498.0,  8660.9, 0.3878, 0.2833, 1735.51, 3209.40, 9500.00),
             "k1017_g_SDSS": (  	  "k1017_g_SDSS.txt", 4759.17, 4665.00,  1397.2,  4052.3,  5449.5, 0.2821, 0.2084,  353.25, 3873.91, 5477.29),
             "k1018_r_SDSS": (  	  "k1018_r_SDSS.txt", 6298.72, 6262.50,  1456.8,  5572.6,  7029.4, 0.3719, 0.2924,  497.08, 5460.01, 7062.19),
             "k1019_i_SDSS": (  	  "k1019_i_SDSS.txt", 7655.58, 7655.33,  1549.4,  6940.3,  8489.7, 0.3281, 0.2228,  424.57, 6784.81, 8524.69),
             "k1020_z_SDSS": (  	  "k1020_z_SDSS.txt", 8893.57, 8927.50,   985.7,  8290.0,  9275.8, 0.1928, 0.0939,  201.07, 8033.83, 10000.00),
       "k1021_wrc3_WR_CIII": (      "k1021_wrc3_WR_CIII.txt", 4660.27, 4661.50,    51.6,  4634.5,  4686.1, 0.1996, 0.0923,   10.54, 4620.85, 4700.38),
       "k1022_wr475_WR_475": (      "k1022_wr475_WR_475.txt", 4760.78, 4765.50,    51.1,  4734.7,  4785.9, 0.2250, 0.1006,   11.89, 4721.82, 4802.27),
      "k1023_wrhe2_WR_HeII": (     "k1023_wrhe2_WR_HeII.txt", 4695.00, 4695.50,    50.4,  4669.9,  4720.4, 0.2240, 0.0970,   11.07, 4654.52, 4735.27),
        "k1024_wrc4_WR_CIV": (       "k1024_wrc4_WR_CIV.txt", 5825.66, 5830.00,    41.8,  5804.8,  5846.7, 0.2717, 0.1065,   11.42, 5794.53, 5858.87),
           "k1025_Bw_NDWFS": (  	"k1025_Bw_NDWFS.txt", 4148.40, 4131.00,  1205.4,  3541.6,  4747.1, 0.2891, 0.2241,  312.02, 3453.62, 4775.17),
         "k1027_823_823_v1": (        "k1027_823_823_v1.txt", 8232.35, 8232.00,    69.9,  8197.3,  8267.3, 0.2596, 0.1087,   17.94, 8178.55, 8285.36),
       "k1040_VR_Bernstein": (      "k1040_VR_Bernstein.txt", 6007.82, 6090.00,  2060.0,  4943.4,  7003.5, 0.3823, 0.2519,  688.24, 4779.41, 7147.25),
         "k1041_Un_Steidel": (        "k1041_Un_Steidel.txt", 3569.63, 3591.03,   496.9,  3326.3,  3823.2, 0.1860, 0.1070,   85.38, 3197.71, 3938.25),
         "k1042_Gn_Steidel": (        "k1042_Gn_Steidel.txt", 4799.43, 4819.02,   839.0,  4365.8,  5204.8, 0.2567, 0.1768,  183.29, 4338.72, 5265.78),
         "k1043_Rs_Steidel": (        "k1043_Rs_Steidel.txt", 6942.58, 7085.38,  1047.2,  6446.2,  7493.4, 0.3605, 0.2096,  349.47, 6300.32, 7560.98),
    "k1045_Ud_Dey_custom_U": (   "k1045_Ud_Dey_custom_U.txt", 3514.97, 3482.00,   272.3,  3390.0,  3662.3, 0.1754, 0.0782,   46.84, 3243.32, 3717.05),
         "k1046_815_815_v2": (        "k1046_815_815_v2.txt", 8151.55, 8148.00,    68.9,  8117.7,  8186.7, 0.2248, 0.1031,   15.50, 8097.04, 8203.28),
           "k1051_337_BATC": (  	"k1051_337_BATC.txt", 3384.58, 3392.51,   222.8,  3274.9,  3497.8, 0.0757, 0.0370,   17.01, 3184.82, 3564.86),
           "k1054_454_BATC": (  	"k1054_454_BATC.txt", 4514.50, 4524.00,   364.8,  4328.0,  4692.8, 0.2358, 0.1478,   77.93, 4290.51, 4736.08),
           "k1055_493_BATC": (  	"k1055_493_BATC.txt", 4896.78, 4908.50,   323.6,  4729.7,  5053.4, 0.2462, 0.1549,   76.58, 4699.82, 5096.68),
           "k1056_527_BATC": (  	"k1056_527_BATC.txt", 5238.93, 5238.50,   340.9,  5075.8,  5416.8, 0.2600, 0.1755,   82.86, 5024.22, 5444.26),
           "k1057_579_BATC": (  	"k1057_579_BATC.txt", 5811.26, 5781.86,   254.8,  5686.0,  5940.9, 0.3281, 0.1792,   79.55, 5659.44, 5969.39),
           "k1058_607_BATC": (  	"k1058_607_BATC.txt", 6084.32, 6086.00,   272.2,  5948.1,  6220.3, 0.3438, 0.2206,   89.79, 5921.12, 6249.89),
           "k1059_666_BATC": (  	"k1059_666_BATC.txt", 6706.62, 6721.50,   429.6,  6498.3,  6928.0, 0.3493, 0.2047,  140.92, 6449.43, 6980.19),
           "k1060_705_BATC": (  	"k1060_705_BATC.txt", 7031.02, 7034.50,   173.4,  6943.9,  7117.4, 0.3321, 0.1379,   58.84, 6899.64, 7165.97),
           "k1061_755_BATC": (  	"k1061_755_BATC.txt", 7543.91, 7557.00,   224.2,  7441.8,  7666.1, 0.2662, 0.1577,   47.49, 7426.94, 7683.38),
           "k1062_802_BATC": (  	"k1062_802_BATC.txt", 8043.15, 8047.50,   231.4,  7925.6,  8157.1, 0.2511, 0.1473,   56.02, 7897.93, 8190.77),
           "k1063_848_BATC": (  	"k1063_848_BATC.txt", 8502.59, 8502.50,   175.1,  8416.8,  8592.0, 0.2206, 0.1376,   37.15, 8387.81, 8620.26),
           "k1064_918_BATC": (  	"k1064_918_BATC.txt", 9164.06, 9150.82,   254.6,  9040.7,  9295.4, 0.0962, 0.0562,   22.45, 9005.83, 9327.37),
           "k1065_973_BATC": (  	"k1065_973_BATC.txt", 9711.74, 9715.00,   261.4,  9587.1,  9848.6, 0.0532, 0.0332,   12.53, 9547.24, 9886.48),


}

sdss_photometric_column = {"u":  2,
                           "g":  4,
                           "r":  6,
                           "i":  8,
                           "z": 10,
                           }
ucac_photometric_column = {"UCAC_Red": 2,
                           "B":  4,
                           "V":  6,
                           "g":  8,
                           "r": 10,
                           "i": 12,
                           }

cellmode_ids = {
    "S": 0,
    "V": 1,
    "D": 1,
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



#
# Source Extractor flags
# 
# 0b00000001 =   1: near neighbor
# 0b00000010 =   2: deblended source
# 0b00000100 =   4: >1 pixel saturated
# 0b00001000 =   8: truncated / image boundary
# 0b00010000 =  16: aperture data incomplete/corrupted
# 0b00100000 =  32: isophotal data incomplete/corrupted
# 0b01000000 =  64: memory overflow during deblending
# 0b10000000 = 128: memory overflow during extraction
#
# The following flags define the FLAGS that exclude a source
#
# WCS can handle saturated/deblended/crowded sources
sexflag_wcs  = 0b11111000
#
# Photometric calibration needs perfect stars
sexflag_phot = 0b11111111



#
# Define some names for the columns in the source extractor catalog
#
SXapertures = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0] #arcsec
SXcolumn_names = [
    'ra', 'dec',
    'x', 'y', 
    'fwhm_image', 'fwhm_world',
    'background',
    'flags',
    'ota',
    'mag_aper_2.0',
    'mag_aper_3.0',
    'mag_aper_4.0',
    'mag_aper_5.0',
    'mag_aper_6.0',
    'mag_aper_8.0',
    'mag_aper_10.0',
    'mag_aper_12.0',
    'mag_err_2.0',
    'mag_err_3.0',
    'mag_err_4.0',
    'mag_err_5.0',
    'mag_err_6.0',
    'mag_err_8.0',
    'mag_err_10.0',
    'mag_err_12.0',
    'mag_auto',
    'mag_err_auto',
    'flux_max',
    'major_axis',
    'minor_axis',
    'position_angle',
    'elongation',
    'ellipticity',
]
# Convert columns into dictionary to make look-up easier
SXcolumn = {}
for name in SXcolumn_names:
    SXcolumn[name] = len(SXcolumn)

SXcolumn_descriptions = [
    "ALPHAWIN_J2000      Right ascension of barycenter (J2000)                      [deg]",
    "DELTAWIN_J2000      Declination of barycenter (J2000)                          [deg]",
    "XWIN_IMAGE          Object position along x (double precision)                 [pixel]",
    "YWIN_IMAGE          Object position along y (double precision)                 [pixel]",
    "FWHM_IMAGE          FWHM assuming a gaussian core                              [pixel]",
    "FWHM_WORLD          FWHM assuming a gaussian core                              [pixel]",
    "BACKGROUND          Background at centroid position [count]",
    "FLAGS               Extraction flags",
    "EXT_NUMBER          FITS extension number",
    "MAG_APER_2.0        magnitude in aperture with radius 2.0''                    [mag]",
    "MAG_APER_3.0        magnitude in aperture with radius 3.0''                    [mag]",
    "MAG_APER_4.0        magnitude in aperture with radius 4.0''                    [mag]",
    "MAG_APER_5.0        magnitude in aperture with radius 5.0''                    [mag]",
    "MAG_APER_6.0        magnitude in aperture with radius 6.0''                    [mag]",
    "MAG_APER_8.0        magnitude in aperture with radius 8.0''                    [mag]",
    "MAG_APER_10.0       magnitude in aperture with radius 10.0''                   [mag]",
    "MAG_APER_12.0       magnitude in aperture with radius 12.0''                   [mag]",
    "MAG_ERR_2.0         magnitude uncertainty in aperture with radius 2.0''        [mag]",
    "MAG_ERR_3.0         magnitude uncertainty in aperture with radius 3.0''        [mag]",
    "MAG_ERR_4.0         magnitude uncertainty in aperture with radius 4.0''        [mag]",
    "MAG_ERR_5.0         magnitude uncertainty in aperture with radius 5.0''        [mag]",
    "MAG_ERR_6.0         magnitude uncertainty in aperture with radius 6.0''        [mag]",
    "MAG_ERR_8.0         magnitude uncertainty in aperture with radius 8.0''        [mag]",
    "MAG_ERR_10.0        magnitude uncertainty in aperture with radius 10.0''       [mag]",
    "MAG_ERR_12.0        magnitude uncertainty in aperture with radius 12.0''       [mag]",
    "MAG_AUTO            Kron-like elliptical aperture magnitude                    [mag]",
    "MAGERR_AUTO         RMS error for AUTO magnitude                               [mag]",
    "FLUX_MAX            Peak flux above background                                 [count]",
    "AWIN_IMAGE          Profile RMS along major axis                               [pixel]",
    "BWIN_IMAGE          Profile RMS along minor axis                               [pixel]",
    "THETA_IMAGE         Position angle (CCW/x)                                     [deg]",
    "ELONGATION          A_IMAGE/B_IMAGE",
    "ELLIPTICITY         1 - B_IMAGE/A_IMAGE",
]

   
reference_zeropoint = {
    "odi_g": [26.0, 26.2],
    "odi_r": [26.1, 26.3],
    "odi_i": [25.6, 25.8],
    "odi_z": [24.6, 24.8],
}
atm_extinction = {
    "odi_g": 0.2,
    "odi_r": 0.12,
    "odi_i": 0.058,
    "odi_z": 0.04,
}
photzp_colorterms = {
    "odi_g": [ 0.1600, 'g', 'r'],
#   "odi_r": [ 0.0047, 'g', 'i'],
#   "odi_r": [ 0.0074, 'g', 'r'],
#   "odi_r": [-0.0001, 'r', 'i'],
#   "odi_i": [-0.0027, 'r', 'i'],
#   "odi_i": [ 0.0022, 'i', 'z'],
#   "odi_i": [-0.0006, 'r', 'z'],
    "odi_z": [-0.1277, 'i', 'z'],
}


#
# Header keyword names in the TECHDATA extension
# to get the actual keywords, add the cell-id string, 
# i.e. "%02d%d%d" % (ota, cellx, celly)
#
techdata_keywords = [
    "GN__", "GN_E",
    "RN__", "RN_E",
    "RNE_", "RNEE", 
]
backup_gain = 1.3
backup_readnoise = 6.5
backup_readnoise_electrons = 9.0

invalid_sky_level_value = -99999.99


def get_valid_filter_name(hdr):
    """

    Convert the header FILTERID entry into a valid filtername. This is necessary
    as some of the filter names change over time, but the filter-id should
    remain unchanged.

    """

    try:
        filter_id = hdr['FILTERID'].strip()
        if (filter_id in list_of_valid_filter_ids):
            filter_name = list_of_valid_filter_ids[filter_id]
            return filter_name
    except:
        try: 
            filter = hdr['FILTER']
            return filter
        except:
            return "unknown"
    return 'unknown'





def get_cellmode(primhdr, cellhdr, focalplanelayout):
    """

    Check if the specified cell, identified by OTA and CELL-ID, is either broken
    or marked as a guide/video cell.

    """

    ota = int(primhdr['FPPOS'][2:4])

    ota_name = "OTA%02d" % ota
    extname = "OTA%02d.SCI" % ota

    cell = cellhdr['EXTNAME']

    # Check if this is one of the permanently broken cells
    wm_cellx, wm_celly = cellhdr['WN_CELLX'], cellhdr['WN_CELLY']
    broken = focalplanelayout.is_cell_broken(primhdr['OTA_ID'], wm_cellx, wm_celly)
    if (broken):
        return -1

    # It's not one of the broken cells, but it might still be a guide/video cell
    idx = wm_cellx + 8 * wm_celly
    cellmode = primhdr['CELLMODE']
    this_cellmode = cellmode[idx]
    cell_id = cellmode_ids[this_cellmode]

    return cell_id


def stdout_write(str):
    """

    Write a given text to stdout and flush the terminal. Pretty much what print
    does, but flushing the output afterwards.

    """

    sys.stdout.write(str)
    sys.stdout.flush()
    return

def clobberfile(filename):
    """

    Delete a file if it already exists, otherwise do nothing.

    """

    if (os.path.isfile(filename)):
        os.remove(filename)
    return



from types import *   
def shmem_as_ndarray( raw_array ):
    """

    Helper function needed to allocate shared memory numpy-arrays.

    """
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




def sexa2deg(sexa):
    """

    Convert a sexa-decimal coordinate string (e.g. 12:30:00.0") into a float
    (12.5 in the example).

    """

    components = sexa.split(":")
    
    if (len(components) != 3):
        return 1e99
    
    deg = float(components[0])
    min = float(components[1])
    sec = float(components[2])
    
    return math.copysign(math.fabs(deg) + math.fabs(min/60.0) + math.fabs(sec/3600.0), deg)

def deg2sexa(deg, signed=False):
    """

    Convert a float coordinate into the more user-friednly sexa-decimal format

    """

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


    


def inherit_headers(header, primary_header):
    """Copy all headers from the primary header into the current header

    """

    for header in headers_to_inherit:
        if (not header in primary_header):
            print "Problem with header ",header
            continue

        card = primary_header.ascardlist()[header]
        header[card.key] = (card.value, card.comment)

        

def rebin_image(data, binfac, operation=numpy.mean):
    """

    Apply a binning factor to a data array.

    Parameters
    ----------
        
    data : ndarray

        Input data array. Only tested to work on two-dimensional data.

    binfac : int

        binning factor, e.g. 2 for 2xw binning. Only identical binning in both
        dimensions is supported at the present.

    operation : function (default: numpy.mean)

        What operation to use when combining the pixels. All functions operating
        on ndarrays are supported, but typical cases are likely one of the
        following:

        * numpy.mean
        * numpy.average
        * numpy.sum
        * numpy.median

        Or from the bottleneck package:

        * bottleneck.nanmean
        * bottleneck.nanmedian

    """
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
    rb2 = operation(rb1, axis=-1)
    rebinned = operation(rb2, axis=1)

#    rb1 = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
#    rb2 = nanmedian(rb1, axis=-1)
#    rebinned = nanmedian(rb2, axis=1)

#    #.nanmean(axis=-1).nanmean(axis=1)

    return rebinned



def center_coords(hdr):
    """

    Return the center coordinates of a given HDU, based on the WCS information
    (does not include distortion).

    """
    
    try:
        centerx, centery = hdr['NAXIS1']/2, hdr['NAXIS2']/2
    except:
        centerx, centery = 2048., 2048.

    center_ra  = (centerx-hdr['CRPIX1'])*hdr['CD1_1'] + (centery-hdr['CRPIX2'])*hdr['CD1_2'] + hdr['CRVAL1']
    center_dec = (centerx-hdr['CRPIX1'])*hdr['CD2_1'] + (centery-hdr['CRPIX2'])*hdr['CD2_2'] + hdr['CRVAL2']

    return center_ra, center_dec

    

def break_region_string(str_region):
    """

    Break down a IRAF-like string (e.g. [0:100,0:200]) into its four components
    and return them separately.

    """
    reg = str_region[1:-1]
    x,dummy,y = reg.partition(",")
    x1,dummy,x2 = x.partition(":")
    y1,dummy,y2 = y.partition(":")
    return int(x1)-1, int(x2)-1, int(y1)-1, int(y2)-1

def extract_region(data, str_region):
    """Extract a region based on the a IRAF-like region definition.

    See also
    --------
    `break_region_string`

    """
    x1,x2,y1,y2 = break_region_string(str_region)
    return data[y1:y2+1, x1:x2+1]


def insert_into_array(data, from_region, target, target_region):
    """

    Copy data from one array into another array, with source and target
    coordinates specified in the IRAF-style format.

    """
    fx1, fx2, fy1, fy2 = break_region_string(from_region)
    tx1, tx2, ty1, ty2 = break_region_string(target_region)

    if (fx2-fx1 != tx2-tx1 or fy2-fy1 != ty2-ty1):
        print "Dimensions do not match, doing nothing"
    else:
        target[ty1:ty2+1, tx1:tx2+1] = data[fy1:fy2+1, fx1:fx2+1]

    return 0

def mask_broken_regions(datablock, regionfile, verbose=False):
    """

    Mask out broken regions in a data array. Regions are defined via a ds9
    region file to allow for ay creation by the user.

    """



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


def is_image_extension(hdu):
    """

    Check if a given HDU is a Image extension

    """
    if (hdu == None):
        return False

    if (type(hdu) == pyfits.hdu.image.ImageHDU or
        type(hdu) == pyfits.hdu.compressed.CompImageHDU):

        try:
            if (hdu.data == None):
                return False
        except:
            # Can't access the .data block, hence this can't be an image extension
            return False

        return True

    elif (type(hdu) == pyfits.hdu.image.PrimaryHDU):
        # This might or might not be an image
        if ('NAXIS' in hdu.header and
            'NAXIS1' in hdu.header and
            'NAXIS2' in hdu.header and
            hdu.header['NAXIS'] == 2 and
            hdu.header['NAXIS1'] > 0 and
            hdu.header['NAXIS2'] > 0):
            # This looks like it might be an image
            return True
            
    return False



def get_svn_version():
    """

    Return current SVN version as given by the `svnversion` command.

    """

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
    """

    Add SVN version number to FITS header.

    """
    svn = get_svn_version()
    hdr["QPIPESVN"] = (svn, "QuickReduce Revision")
    return



def rotate_around_center(data, angle, mask_limit = 0.1, verbose=True, safety=1, mask_nans=True, spline_order=3):
    """

    Rotate a given data array. Rotation center is the center of the data frame.

    """
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
    """

    Return the level of the installed filter, based on the information in the
    FLTARM header keywords. If more than one filter is active, this function
    returns the level of the lowest populated filter arm.

    """
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
    """
    Don't use anymore!
    """
    if (data_in == None):
        return numpy.zeros(shape=(494,480))

    return data_in[0:494, 0:480]

def cell2ota__get_target_region(x, y, binning=1):
    """

    Get the location of a given cell in the monolithic OTA array, accounting for
    binning (only 1x1 and 2x2 supported).

    """
    #taken from ODI OTA Technical Datasheet (det area 480x494, streets 11/28 px)

    # Y-coordinates of cell numbers are counted top down,
    # but pixel coordinates are counted bottom up
    _y = 7 - y

    if (binning == 1):
        y1 = (505*_y)  #was 503 
        y2 = y1 + 494
        x1 = 508 * x
        x2 = x1 + 480
    elif (binning == 2):
        y1 = int(math.ceil(252.5 * _y))
        y2 = y1 + 247
        x1 = 254 * x
        x2 = x1 + 240
        
    return x1, x2, y1, y2




def three_sigma_clip(input, ranges=[-1e9,1e9], nrep=3, return_mask=False, nsigma=3):
    """

    Perfom an iterative 3-sigma clipping on the passed data array. 

    """

    old_valid = numpy.isfinite(input)
    valid = (input > ranges[0]) & (input < ranges[1])
                
    for rep in range(nrep):
        if (numpy.sum(valid) < 1):
            valid = old_valid
            break

        lsig = scipy.stats.scoreatpercentile(input[valid], 16)
        hsig = scipy.stats.scoreatpercentile(input[valid], 84)
        median = numpy.median(input[valid])
        sigma = 0.5 * (hsig - lsig)

        mingood = numpy.max([median - nsigma*sigma, ranges[0]])
        maxgood = numpy.min([median + nsigma*sigma, ranges[1]])

            #print median, sigma
        old_valid = valid
        valid = (input > mingood) & (input < maxgood)
        
    if (return_mask):
        return input[valid], valid

    return input[valid]


def derive_ota_outlines(otalist):
    """

    For each OTA (extension) in the pased list, derive the sky-position of all 4
    corners.

    """
    from astLib import astWCS

    all_corners = []
    for ext in range(len(otalist)):
        if (not is_image_extension(otalist[ext])):
            continue
#type(otalist[ext]) == pyfits.hdu.image.ImageHDU):
            
        wcs = astWCS.WCS(otalist[ext].header, mode='pyfits')
            
        corner_coords = []
        corner_coords.append(wcs.pix2wcs(                            0,                             0))
        corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'],                             0))
        corner_coords.append(wcs.pix2wcs(otalist[ext].header['NAXIS1'], otalist[ext].header['NAXIS2']))
        corner_coords.append(wcs.pix2wcs(                            0, otalist[ext].header['NAXIS2']))

        all_corners.append(corner_coords)

    return all_corners


def create_qa_filename(outputfile, plotname, options):
    """

    Return the filename for a given diagnostic plot, accounting for
    user-specified preferences.

    """
    if (options['structure_qa_subdirs']):
        dirname, basename = os.path.split(outputfile)
        if (dirname == None or dirname == ''): 
            dirname = "."
        qa_plotdir = "%s/%s/" % (dirname, options['structure_qa_subdir_name'])
        if (not os.path.isdir(qa_plotdir)):
            os.mkdir(qa_plotdir)
        qa_plotfile = "%s/%s" % (qa_plotdir, plotname)
    else:
        qa_plotfile = "%s.%s" % (outputfile[:-5], plotname)

    return qa_plotfile

def create_qa_otaplot_filename(plotname, ota, structure_qa_subdirs):
    """

    Return the filename for a given OTA-level diagnostic plot, accounting for
    user-specified preferences.

    """

    if (structure_qa_subdirs):
        # in this case, plotname is to be taken as directory
        if (not os.path.isdir(plotname)):
            os.mkdir(plotname)
        # The actual filenames are the directory and the OTA
        qa_plotfile = "%s/OTA%02d" % (plotname, ota)
    else:
        qa_plotfile = "%s_OTA%02d" % (plotname, ota)

    return qa_plotfile








#
# Additional functions to handle binned data more elegantly
#
def get_binning(hdr):
    """

    Get the binning factor of a given frame/cell based on its header.

    """
    # print "determining binning factor"
    
    # The CCDBIN1 keyword should be foung in primary header, so try this one first
    try:
        binfac = hdr['CCDBIN1']
        # print binfac
        return int(binfac)
    except:
        pass

    # Couldn't find the CCDBIN1 header, so this is likely not a primary HDU
    # Try the CCDSUM keyword next
    try:
        ccdsum = hdr['CCDSUM']
        # print ccdsum
        items = ccdsum.split()
        return int(items[0])
    except:
        pass

    #
    # If we still don't find a valid header, look at the size of the data section
    #

    # Still no luck finding a header? Assume it's 1 and continue
    return 1


def get_collected_image_dimensions(binning):
    """

    Return the dimension of the monolithic OTA frame, accounting for binning

    """
    if (binning == 1):
        sizex, sizey = 4096, 4096
    elif (binning == 2):
        sizex, sizey = 2048, 2048

    return sizex, sizey


def extract_datasec_from_cell(data, binning):
    """

    Return the science region of a given cell, accounting for binning

    """

    if (binning == 1):
        dx, dy = 480, 494
    elif (binning == 2):
        dx, dy = 240, 247

    # print "extracting datasec", dx, dy
    return data[0:dy, 0:dx]


def extract_biassec_from_cell(data, binning):
    """

    Return the overscan region of a given cell, accounting for binning

    """

    if (binning == 1):
        dx1, dx2, dy1, dy2 = 500, 550, 0, 494
    elif (binning == 2):
        dx1, dx2, dy1, dy2 = 260, 277, 0, 246

    # print dx1, dx2, dy1, dy2
    return data[dy1:dy2, dx1:dx2]


def match_catalog_areas(src, to_match, radius):
    """

    Match the area coverage of two catalogs, allowing for some extra coverage
    around the edges.

    """

    #
    # Account for the cos(dec) effect
    #
    max_dec = numpy.max(numpy.fabs(src[:,1]))
    if (max_dec > 85): max_dec = 85
    cos_dec = numpy.cos(numpy.radians(max_dec))

    #
    # Apply cos(dec) correction to the src catalog 
    # (make sure to copy the catalogs to prevent altering the 
    # actual input catalogs) ...
    #
    src_radec = src[:,0:2].copy()
    src_radec[:,0] *= cos_dec
    #src_radec[:,0] *= numpy.cos(numpy.radians(src_radec[:,1]))
    src_tree = scipy.spatial.cKDTree(src_radec)

    #
    # ... and to the match catalog
    #
    match_radec = to_match[:,0:2].copy()
    match_radec[:,0] *= cos_dec

    #
    # Now determine, for each match source (in 2MASS), the distance 
    # to the 1 (k=1) closest source in the src (ODI) catalog. If the distance is
    # less than 'radius', the source is close enough and is considered matched
    #
    d,l = src_tree.query(match_radec, k=1, eps=0.05, p=2, distance_upper_bound=radius)
    matched = to_match[d<radius]

    return matched






def collect_reduction_files_used(masterlist, files_this_frame):
    """
    Keeps track of all files used during the reduction. This function maintains 
    a list of files and their corresponding reduction step, and properly handles
    individual files as well as list of filenames.
    
    This routine is called during reduction to keep track of all file 
    associations.

    Parameters
    ----------
    masterlist : dictionary

        Dictionary of all reduction steps that have external files associated. 
        For each step it maintains a list of files.

    files_this_frame : dictionary

        Like the master_list, just for only one step. 


    Returns
    -------
    the updated master_list

    """

    for key, value in files_this_frame.iteritems():
        if (key in masterlist):
            existing_keys = masterlist[key]
            if (type(value) == list):
                for val1 in value:
                    masterlist[key].append(val1)
            else:
                masterlist[key].append(value)
        else:
            # This is a new key, so just copy it
            if (type(value) == list):
                masterlist[key] = value
            else:
                masterlist[key] = [value]

    # print masterlist
    return masterlist










def create_association_table(master, verbose=False):
    """

    Convert the association dictionary maintained and updated by 
    :proc:collect_reduction_files_used and creates a FITS-compatible 
    associations table that will be stored in the resulting output FITS file.

    Parameters
    ----------
    master : dictionary

        the master associations dictionary 

    verbose : Bool

        Activate some debugging output

    Returns
    -------
    tbhdu : TableHDU

        A FITS TableHDU containing the assocations table. Each entry in this 
        table contains the reduction step, the name of the associated file, and 
        the full path of this file.

    """

    reduction_step = []
    full_filename = []
    short_filename = []

    for key, value in master.iteritems():
        # print key,":",value
        for filename in set(value):
            reduction_step.append(key)
            if (filename == None):
                continue
            full_filename.append(os.path.abspath(filename))
            
            dirname, filebase = os.path.split(filename)
            short_filename.append(filebase)

            if (verbose):
                print "% 15s : %s" % (key, filename)

    # print reduction_step
    # print short_filename

    columns = [\
        pyfits.Column(name='correction',    format='A25',  array=reduction_step),
        pyfits.Column(name='filename_full', format='A375', array=full_filename),
        pyfits.Column(name='filename',      format='A100', array=short_filename),
        ]

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.BinTableHDU.from_columns(coldefs)

    tbhdu.name = "ASSOCIATIONS"

    return tbhdu



def add_fits_header_title(header, title, before_keyword):
    header.add_blank(before=before_keyword)
    header.add_blank(title, before=before_keyword)
    header.add_blank(before=before_keyword)
    return

def sigma_to_percentile(s):
    return 100 * 0.5*(1 + scipy.special.erf(s/math.sqrt(2)))




def format_filename(input_filename_or_header, outputfile):
    """

    Format the input string and replace special keywords with information from
    the FITS header.

    """

    if (type(input_filename_or_header) == str):
        # This is a string, interpret it as a filename
        if (os.path.isfile(input_filename_or_header)):
            filename = input_filename_or_header
        elif (os.path.isdir(input_filename_or_header)):
            dirname = input_filename_or_header
            if (dirname.endswith("/")): dirname = dirname[:-1]
            # For now, assume that OTA 33 always exists, and that there is no .fz ending
            directory, obsid = os.path.split(dirname)
            filename = "%s/%s.33.fits" % (dirname, obsid)
            if (not os.path.isfile(filename)):
                filename = filename+".fz"
                if (not os.path.isfile(filename)):
                    return outputfile

        hdulist = pyfits.open(filename)
        header = hdulist[0].header
    else:
        # This is a header
        header = input_filename_or_header

    while (outputfile.find("%") >= 0):
        start = outputfile.find("%") 
        if (outputfile[start:start+7] == "%FILTER"):
            outputfile = outputfile[:start] + header['FILTER'] + outputfile[start+7:]
        elif (outputfile[start:start+7] == "%OBJECT"):
            # The object name might contain spaces, replace them with underscores
            # Also replace all other special characters with underscores
            objectname = header['OBJECT'].replace(' ', '_')
            objectname = objectname.replace(',', '_')
            objectname = objectname.replace('(', '_')
            objectname = objectname.replace(')', '_')
            objectname = objectname.replace('/', '_')
            objectname = objectname.replace('\\', '_')
            objectname = objectname.replace('`', '_')
            objectname = objectname.replace('"', '_')
            objectname = objectname.replace('\'', '_')
            objectname = objectname.replace(')', '_')
            outputfile = outputfile[:start] + objectname  + outputfile[start+7:]
        elif (outputfile[start:start+6] == "%OBSID"):
            outputfile = outputfile[:start] + header['OBSID'] + outputfile[start+6:]
        elif (outputfile[start:start+7] == "%OBSSEQ"):
            obsid = header['OBSID']
            dot_position = obsid.find(".")
            obs_sequence = obsid[:dot_position]
            outputfile = outputfile[:start] + obs_sequence + outputfile[start+7:]
        elif (outputfile[start:start+8] == "%EXPTIME"):
            outputfile = "%s%.1f%s" % (outputfile[:start], header['EXPTIME'], outputfile[start+8:])
        else:
            stdout_write("found unknown tag in %s\n" % outputfile)
            break

    return outputfile



