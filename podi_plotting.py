#!/usr/local/bin/python

#
# (c) Ralf Kotulla for WIYN/pODI
#


#
# This file contains all information, routines and imports dealing with plotting
#


#################################3
#
# Some definitions for colors 
#
#################################3
# Set plotting backend to Agg to make it work in terminals without X support

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
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
matplotlib.pyplot.register_cmap(cmap=cmap_bluewhite)
