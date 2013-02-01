#!/usr/bin/env python


import sys
import numpy
import os
import math

from podi_definitions import *
import pyfits
import matplotlib.pyplot as plot
import matplotlib.cm as cm


def view_shift_history(filename):

    # Load the shiftf file
    hdu = pyfits.open(filename)
    dir, name = os.path.split(filename)

    #print hdu.info()
    
    shift_history = hdu[1].data


    if (cmdline_arg_isset("-shiftdelay")):
    
        shift_complete = hdu[1].data.field('shift_complete')
        #numpy.savetxt("/Volumes/ODIScratch/rkotulla/shift.dump", shift_history)

        time_diff = (shift_history.field('shift_complete') - shift_history.field('message_received')) * 1000.
        print time_diff.shape
        print time_diff
        count, bins = numpy.histogram(time_diff, 400, (-1,1))
        #for i in range(count.shape[0]):
        #    print bins[i], bins[i+1], count[i]

        mean_time = 0.5*(bins[:-1] + bins[1:])
        plot.plot(mean_time, count)
        plot.ylabel("Count")
        plot.xlabel("time delay between start and shift [milli-secs]")
        plot.show()
    
    elif (cmdline_arg_isset("-xyhist")):
        print "PLotting xy-histogram of shift positions"
        count, xedges, yedges = numpy.histogram2d(shift_history.field('shift1'), shift_history.field('shift2'),
                                                  bins=[20,20], range=[[-10,10], [-10,10]])

        fig = plot.figure()
        # split filename into directory and filename
        fig.suptitle("Shift history for\n%s" % (name))

        # Draw the 2-d histogram on the right
        img = plot.imshow(count, interpolation='nearest', extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),origin='lower')
        fig.colorbar(img)

        maxcount = numpy.max(count)
        for x in range(count.shape[0]):
            for y in range(count.shape[1]):
                color = "white" if count[x,y] < 0.2*max_count else "black"
                print count[x,y], max_count, color
                plot.text(0.5*(xedges[x]+xedges[x+1]), 0.5*(yedges[y]+yedges[y+1]), "%d" % count[x,y], ha='center', va='center', color=color)

        plot.xlabel("shift_x")
        plot.ylabel("shift_y")
        
        plot.show()


        #if (cmdline_arg_isset("-xyscatter")):
    else:

        n_shifts = shift_history.field('shift1').shape[0]
        # Create plot
        fig = plot.figure(figsize=(12,6))
        fig.canvas.set_window_title("podi_viewshifthistory: Shift history for: %s" % (name))
        fig.suptitle("Shift history for: %s  --  #shifts = %d" % (name, n_shifts ))

        little_scatter_x = numpy.random.rand(n_shifts)/2-0.5
        little_scatter_y = numpy.random.rand(n_shifts)/2-0.5

        # Draw the two plots on the left, bototm one first
        plot.axes([0.05,0.1, 0.38,.35])
        plot.scatter(shift_history.field('framenum'), shift_history.field('shift2')+little_scatter_y)
        plot.grid(True)
        plot.xlabel("frame number")
        plot.ylabel("shift y")
        
        #plot.subplot(2,1,2)
        plot.axes([0.05,0.55, 0.38,0.35])
        plot.scatter(shift_history.field('framenum'), shift_history.field('shift1')+little_scatter_x)
        plot.grid(True)
        plot.xlabel("frame number")
        plot.ylabel("shift x")


        # Now draw the 2-d histogram pn the right

        count, xedges, yedges = numpy.histogram2d(shift_history.field('shift1'), shift_history.field('shift2'),
                                                  bins=[20,20], range=[[-10,10], [-10,10]])

        plot.axes([0.5,0.1,0.5,0.8])
        img = plot.imshow(count, interpolation='nearest', extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),origin='lower')
        fig.colorbar(img)

        maxcount = numpy.max(count)
        for x in range(count.shape[0]):
            for y in range(count.shape[1]):
                color = "#000000" #if count[x,y] < 0.3*maxcount else "black"
                plot.text(0.5*(xedges[x]+xedges[x+1]), 0.5*(yedges[y]+yedges[y+1]), "%d" % count[y,x], ha='center', va='center', color=color)

        plot.xlabel("shift_x")
        plot.ylabel("shift_y")
 

        
        plot.show()
        
    print 
    

    
if __name__ == "__main__":
    
    filename = sys.argv[1]
    
    view_shift_history(filename)
    
