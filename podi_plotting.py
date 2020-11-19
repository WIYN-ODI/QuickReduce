#!/usr/bin/env python3
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

This module contains all information and imports dealing with plotting

"""


#################################3
#
# Some definitions for colors 
#
#################################3
# Set plotting backend to Agg to make it work in terminals without X support

import matplotlib
matplotlib.use('Agg')
import matplotlib.colors
from matplotlib.colors import LinearSegmentedColormap

colorfade_bluewhite = {'red':   ((0.0, 1.0, 1.0),
                                 (1.0, 142./255., 1.)),

                       'green': ((0.0, 1.0, 1.0),
                                 (1.0, 163./255., 0.0)),

                       'blue':  ((0.0, 1.0, 1.0),
                                 (1.0, 218./255., 0.0))
                       }
cmap_bluewhite = matplotlib.colors.LinearSegmentedColormap('BlueWhite', colorfade_bluewhite)

import matplotlib.pyplot
matplotlib.pyplot.register_cmap(cmap=cmap_bluewhite)

import matplotlib.patches
