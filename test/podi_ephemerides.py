#!/usr/bin/env python


# coding: utf-8
import telnetlib
import time
import os, sys
import numpy
import scipy, scipy.interpolate
import matplotlib, matplotlib.pyplot
import datetime

sys.path.append("/work/podi_devel/")
from podi_commandline import *
from podi_definitions import *


def get_ephemerides_for_object(object_name, 
                               start_datetime="2012-01-01",
                               end_datetime="2014-01-01",
                               time_interval="6h",
                               session_log_file=None,
                               verbose=True):

    tn = telnetlib.Telnet('horizons.jpl.nasa.gov', 6775)
    tn.write("vt100\n")

    session_log = ""

    return_string = tn.read_until("Horizons> ")
    session_log += return_string

    # NASA should return:
    # ---
    # JPL Horizons, version 3.91 
    # Type `?' for brief intro, `?!' for more details 
    # System news updated Apr 22, 2014
    # 
    # Horizons> 
    # ---

    tn.write("%s\n" % (object_name))
    # Now read the NASA return.
    # For a unnumbered object, NASA will ask for a confirmation

    done_reading = False
    return_string = ""
    while (not done_reading):
        return_string += tn.read_some()
        lines = return_string.split("\n")
        if (len(lines) >= 4):
            done_reading = True

    session_log += return_string
    if (verbose): print return_string

    if (return_string.find(">EXACT< designation search [CASE & SPACE sensitive]:") >= 0):
        # This was a search for the exact designation

        # Read the rest of the output
        return_string += tn.read_until(":")
        #
        # NASA returns
        # ---
        # Horizons> 2014 CR13
        #  
        # >EXACT< designation search [CASE & SPACE sensitive]:
        #   DES = 2014 CR13;
        # Continue [ <cr>=yes, n=no, ? ] : 
        # ---
        #

        # Send confirmation
        tn.write("\n")

    search_result = tn.read_until("<cr>:")
    session_log += search_result
    print session_log
    if (verbose): print search_result
    # For numbered objects, and for un-numbered objects after confirmation, we
    # should received something like this:
    # ---
    # *******************************************************************************
    # JPL/HORIZONS                 300163 (2006 VW139)           2014-May-21 11:43:05
    # Rec #:300163 (+COV)   Soln.date: 2013-Aug-10_15:34:30    # obs: 150 (2000-2013)
    #
    # FK5/J2000.0 helio. ecliptic osc. elements (au, days, deg., period=Julian yrs): 
    #
    #   EPOCH=  2454959.5 ! 2009-May-08.00 (CT)          Residual RMS= .39144        
    #    EC= .1998703078417539   QR= 2.440143900621482   TP= 2455761.497290819       
    #    OM= 83.22983169924017   W=  281.8850275006057   IN= 3.238832031574682       
    #    A= 3.049685475412755    MA= 211.5793264733121   ADIST= 3.659227050204028    
    #    PER= 5.32587            N= .185063804           ANGMOM= .029434476          
    #    DAN= 2.8121             DDN= 3.05355            L= 5.1333415                
    #    B= -3.1693289           MOID= 1.44190001        TP= 2011-Jul-18.9972908190  
    #
    # Asteroid physical parameters (km, seconds, rotational period in hours):        
    #    GM= n.a.                RAD= n.a.               ROTPER= n.a.                
    #    H= 16.1                 G= .150                 B-V= n.a.                   
    #                            ALBEDO= n.a.            STYP= n.a.                  
    #
    # ASTEROID comments: 
    # 1: soln ref.= JPL#2, OCC=0
    # 2: source=ORB
    # *******************************************************************************
    #  Select ... [A]pproaches, [E]phemeris, [F]tp,[M]ail,[R]edisplay, [S]PK,?,<cr>: 
    # ---

    if (search_result.find("No matches found.") >= 0) :
        # If there's no object found, we will instead read this:
        # ---
        # *******************************************************************************
        # JPL/DASTCOM3           Small-body Index Search Results     2014-May-21 11:44:24
        #
        #  Comet AND asteroid index search:
        #
        #     DES = 2014 cr13;
        #
        #  Matching small-bodies: 
        #     No matches found.
        # *******************************************************************************
        # ---

        tn.close()
        print "\n"*5,"No matches found","\n"*5
        return None


    # No request [E]phemeris
    tn.write("E\n")

    horizon_return = tn.read_until("Observe, Elements, Vectors  [o,e,v,?] :")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("o\n")

    horizon_return = tn.read_until(" Coordinate center [ <id>,coord,geo  ] :")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("695@399\n")

    horizon_return = tn.read_until(" Confirm selected station    [ y/n ] --> ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("y\n")

    horizon_return = tn.read_until(" Starting UT  [>=   1599-Dec-11 23:59] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("2012-01-01\n")

    horizon_return = tn.read_until(" Ending   UT  [<=   2500-Dec-30 23:58] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("2015-01-05\n")

    horizon_return = tn.read_until(" Output interval [ex: 10m, 1h, 1d, ? ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("1d\n")

    horizon_return = tn.read_until(" Accept default output [ cr=(y), n, ?] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("n\n")

    horizon_return = tn.read_until(" Select table quantities [ <#,#..>, ?] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("1,3,7,9\n")

    horizon_return = tn.read_until(" Output reference frame [J2000, B1950] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("J2000\n")

    horizon_return = tn.read_until(" Time-zone correction   [ UT=00:00,? ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("\n")

    horizon_return = tn.read_until(" Output UT time format   [JD,CAL,BOTH] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("JD\n")

    horizon_return = tn.read_until(" Output time digits  [MIN,SEC,FRACSEC] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("MIN\n")

    horizon_return = tn.read_until(" Output R.A. format       [ HMS, DEG ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("DEG\n")

    horizon_return = tn.read_until(" Output high precision RA/DEC [YES,NO] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("YES\n")

    horizon_return = tn.read_until(" Output APPARENT [ Airless,Refracted ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("Airless\n")

    horizon_return = tn.read_until(" Set units for RANGE output [ KM, AU ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("AU\n")

    horizon_return = tn.read_until(" Suppress RANGE_RATE output [ YES,NO ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("YES\n")

    horizon_return = tn.read_until(" Minimum elevation [ -90 <= elv <= 90] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("\n")

    horizon_return = tn.read_until(" Maximum air-mass  [ 1 <=   a  <= 38 ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("\n")

    horizon_return = tn.read_until(" Print rise-transit-set only [N,T,G,R] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("N\n")

    horizon_return = tn.read_until(" Skip printout during daylight [ Y,N ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("N\n")

    horizon_return = tn.read_until(" Solar elongation cut-off   [ 0, 180 ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("\n")

    horizon_return = tn.read_until(" Local Hour Angle cut-off       [0-12] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("\n")

    horizon_return = tn.read_until(" Spreadsheet CSV format        [ Y,N ] : ")
    session_log += horizon_return 
    if (verbose): print horizon_return
    tn.write("Y\n")

    horizon_return = tn.read_until("$$SOE\r\n")
    session_log += horizon_return 
    if (verbose): print horizon_return

    ephemdata = tn.read_until("$$EOE")
    session_log += ephemdata
    print ephemdata

    datafile = open("telnet.csv", "w")
    print >>datafile, ephemdata
    datafile.close()
    
    horizon_return = tn.read_until(" >>> Select... [A]gain, [N]ew-case, [F]tp, [K]ermit, [M]ail, [R]edisplay, ? : ")
    session_log += horizon_return 
    if (verbose): print horizon_return

    # send the quit command
    tn.write("q\n")

    horizon_return = tn.read_all()
    session_log += horizon_return 
    if (verbose): print horizon_return

    tn.close()

    if (not session_log_file == None):
        logfile = open(session_log_file, "w")
        print >>logfile, session_log
        logfile.close()

    #
    # Analyse the output 
    #
    jd2mjd = 2400000.5
    np = []
    for line in ephemdata.split("\n")[:-1]:
        items = line.split(',')
        mjd = float(items[0]) - jd2mjd
        ra = float(items[3])
        dec = float(items[4])
        rate_ra = float(items[5])
        rate_dec = float(items[6])
        np.append([mjd,ra,dec,rate_ra,rate_dec])
    data = numpy.array(np)
    numpy.savetxt("ephem.data", data)

    # Now create some interpolation vectors
    ra_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,1], kind='linear' )
    dec_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,2], kind='linear' )
    rate_ra_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,3], kind='linear' )
    rate_dec_vs_mjd = scipy.interpolate.interp1d( data[:,0], data[:,4], kind='linear' )
    results = {
        'data': data,
        'ra': ra_vs_mjd,
        'dec': dec_vs_mjd,
        'rate_ra': rate_ra_vs_mjd,
        'rate_ra': rate_dec_vs_mjd,
        }
    #
    #
    #

    return results


def find_timescale(filelist):

    mjd_list = []
    for fitsfile in filelist:
        if (not os.path.isfile(fitsfile)):
            continue

        hdulist = pyfits.open(fitsfile)
        for ext in hdulist:
            if ('MJD-OBS' in ext.header):
                obs_mjd = ext.header['MJD-OBS']
                break
        hdulist.close()

        mjd_list.append(obs_mjd)

    mjd_list = numpy.array(mjd_list)

    # Now determine the min and max times
    mjd_min = numpy.min(mjd_list)
    mjd_max = numpy.max(mjd_list)

    # print mjd_min, mjd_max

    return mjd_min, mjd_max


def mjd_to_datetime(mjd, fmt):

    mjd_to_y2k = 51544.

    deltatime = datetime.timedelta(days=mjd-mjd_to_y2k)

    mjd_ref_date = datetime.date(year=2000, month=1, day=1)
    mjd_ref_time = datetime.time(0,0)
    mjd_ref_datetime = datetime.datetime(2000,1,1,0,0)

    real_time = mjd_ref_datetime + deltatime

    # print mjd_ref_date.isoformat()
    # print mjd_ref_date.strftime(fmt)

    # print real_time.strftime(fmt)

    return real_time.strftime(fmt)


def get_ephemerides_for_object_from_filelist(object_name, 
                                             filelist,
                                             session_log_file=None,
                                             verbose=True):

    print object_name
    print "\n".join(filelist)+"\n"

    mjd_min, mjd_max = find_timescale(filelist)

    delta_mjd = mjd_max - mjd_min
    time_steps = "1h"

    if (delta_mjd > 100):
        timebuffer = 5.
        time_steps = "6h"
    elif (delta_mjd > 5):
        timebuffer = 1.
        time_steps = "4h"
    elif (delta_mjd >= 0.5):
        timebuffer = 0.1
        time_steps = "1h"
    else:
        timebuffer = 0.1
        time_steps = "30min"
    
    # times need to be in format 1599-Dec-11 23:59
    start_mjd = math.floor((mjd_min-timebuffer)/timebuffer) * timebuffer
    end_mjd   =  math.ceil((mjd_max+timebuffer)/timebuffer) * timebuffer

    print start_mjd, end_mjd

    start_datetime = mjd_to_datetime(start_mjd, "%Y-%m-%d %H:%M")
    end_datetime = mjd_to_datetime(end_mjd, "%Y-%m-%d %H:%M")

    print start_datetime, end_datetime

    results =  get_ephemerides_for_object(object_name, 
                                          start_datetime=start_datetime,
                                          end_datetime=end_datetime,
                                          time_interval=time_steps,
                                          session_log_file=session_log_file,
                                          verbose=verbose)

    return results
   
if __name__ == "__main__":

    if (cmdline_arg_isset("-times")):
        filelist = get_clean_cmdline()[1:]
        find_timescale(filelist)

    if (cmdline_arg_isset("-mjd2date")):
        mjd = float(get_clean_cmdline()[1])
        print mjd_to_datetime(mjd, "%Y-%m-%d %H:%M:%S")

    if (cmdline_arg_isset("-fromlist")):
        objname = get_clean_cmdline()[1]
        filelist = get_clean_cmdline()[2:]
        results = get_ephemerides_for_object_from_filelist(object_name=objname, 
                                                           filelist=filelist,
                                                           session_log_file="session.log",
                                                           verbose=False)
        print results

    else:
        object_name = get_clean_cmdline()[1]
        get_ephemerides_for_object(object_name, verbose=False)

