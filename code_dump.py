

improve_wcs = None #"CRVAL" # or "CRPIX" or "shift" or None

        # matched_cat contains the reference coordinates and some info on source pairing
        #if (improve_wcs != None):
        #    improve_wcs_solution(ota_list, matched_cat, fixwcs_odi_sourcecat)


def wcserr_crpix(p, xy, radec, wcs):

    #print "\n\n\n\n"
    #print "p=",p

    wcs.header['CRPIX1'] = p[0]
    wcs.header['CRPIX2'] = p[1]

    if (len(p) > 2):
        wcs.header['CD1_1'] = p[2]
        wcs.header['CD1_2'] = p[3]
        wcs.header['CD2_1'] = p[4]
        wcs.header['CD2_2'] = p[5]

    wcs.updateFromHeader()

    xy_comp = wcs2odi(radec, wcs)
    dxy = xy_comp - xy
    #print dxy[0:10]
    err2 = numpy.sqrt(numpy.sum(dxy**2,axis=1)) #* 3600.0
    return err2

def wcserr_crval(p, xy, radec, wcs):

    #print "\n\n\n\n"
    #print "p=",p

    wcs.header['CRVAL1'] = p[0]
    wcs.header['CRVAL2'] = p[1]

    wcs.updateFromHeader()

    radec_comp = odi2wcs(xy, wcs)
    dradec = radec_comp - radec

    dradec[:,0] *= numpy.cos(numpy.radians(radec[:,1]))

    #print "median_radec err (crval)", numpy.median(dradec, axis=0), "(",p,")"

    err2 = numpy.sqrt(numpy.sum(dradec**2,axis=1)) * 3600.0
    return err2

def wcserr_shift(p, ref, odi):
    odi_x = odi + p
    err2 = numpy.sqrt(numpy.sum((ref-odi_x)**2,axis=1))
    return err2

def improve_wcs_solution(ota_list, matched_cat, fixwcs_odi_sourcecat):

    print "Improving WCS solution..."
    for ext in range(1, len(ota_list)):
        ota = int(ota_list[ext].header['EXTNAME'][3:5])

        #print wcs.header

        # Now select all stars from this extension

        in_this_ota = fixwcs_odi_sourcecat[:,12] == ota
        
        matched_ota = matched_cat[in_this_ota]
        sources_ota = fixwcs_odi_sourcecat[in_this_ota]

        # Pick sources that are within a matching radius
        match_distance = numpy.sqrt((matched_ota[:,0] - matched_ota[:,2])**2 + (matched_ota[:,1] - matched_ota[:,3])**2)
        nearby_matches = match_distance < 5 * (1./3600.)

        # Merge both catalogs and dump to temporary file
        if (False):
            print "Writing matchcat files"
            dump_matchcat = numpy.zeros((nearby_matches.shape[0],4))
            dump_matchcat[:,0:2] = sources_ota[:,2:4]
            #dump_matchcat[:,2:4] = matched_ota[:,0:2]
            dump_matchcat[:,2:4] = matched_ota[:,2:4]
            x = open("matchcat_%02d.dat" % ota, "w")
            numpy.savetxt(x, dump_matchcat)
            x.close()
            del x

        if (numpy.sum(nearby_matches) < 5):
            continue

        matched_good = matched_ota[nearby_matches]
        sources_good = sources_ota[nearby_matches]

        n_matches = sources_good.shape[0]

        # Now we have to match up pixel coordinates with sky-coordinates
        xy = sources_good[:,2:4]
        radec = matched_good[:, 0:2]

        ### if (False):
        ###     pixel = 0.11 / 3600.
        ###     nrepeat = 0
        ###     mederr = 1e9
        ###     while (mederr > 1*pixel and nrepeat < 2):

        ###         xy_comp = wcs2odi(radec, wcs)
        ###         dxy = xy_comp - xy
        ###         print dxy

        ###         median_shift = numpy.median(dxy, axis=0)
        ###         print median_shift

        ###         wcs.header['CRPIX1'] -= median_shift[0]
        ###         wcs.header['CRPIX2'] -= median_shift[1]

        ###         nrepeat += 1

        ###     ota_list[ext].header['CRPIX1'] = wcs.header['CRPIX1']
        ###     ota_list[ext].header['CRPIX2'] = wcs.header['CRPIX2']


        wcs = astWCS.WCS(ota_list[ext].header, mode="pyfits")

        if (improve_wcs == "CRPIX"):
            errfunc = lambda p, xy, radec, wcs: wcserr_crpix(p, xy, radec, wcs)
            pinit = [wcs.header['CRPIX1'], wcs.header['CRPIX2']]

            out = scipy.optimize.leastsq(errfunc, pinit,
                                         args=(xy, radec, wcs), full_output=1)

            ota_list[ext].header['CRPIX1'] = out[0][0]
            ota_list[ext].header['CRPIX2'] = out[0][1]
        elif (improve_wcs == "CRVAL"):
            wcs = astWCS.WCS(ota_list[ext].header, mode="pyfits")
            errfunc = lambda p, xy, radec, wcs: wcserr_crval(p, xy, radec, wcs)
            pinit = [wcs.header['CRVAL1']+(2./3600.), wcs.header['CRVAL2']-(1./3600.)]

            out = scipy.optimize.leastsq(errfunc, pinit,
                                         args=(xy, radec, wcs), full_output=1)

            d_crval1, d_crval2 = out[0][0]-pinit[0], out[0][1]-pinit[1]
            ota_list[ext].header.update("D_CRVAL1", d_crval1)
            ota_list[ext].header.update("D_CRVAL2", d_crval2)

            print "BF_CRVAL:",out[0][0], wcs.header['CRVAL1']

            ota_list[ext].header['CRVAL1'] = out[0][0]
            ota_list[ext].header['CRVAL2'] = out[0][1]
        elif (improve_wcs == "shift"):

            ref = matched_good[:, 0:2]
            odi = matched_good[:, 2:4]
            errfunc = lambda p, ref, odi: wcserr_shift(p, ref, odi)

            pinit = [0, 0] 

            out = scipy.optimize.leastsq(errfunc, pinit,
                                         args=(ref, odi), full_output=1)

            print "Final shift", out[0], out[0]*3600.
            ota_list[ext].header['CRVAL1'] -= out[0][0]
            ota_list[ext].header['CRVAL2'] -= out[0][1]

        else:
            return

        radec_computed = odi2wcs(xy, wcs)
        d_radec = radec - radec_computed

        d_radec += numpy.random.randn(d_radec.shape[0], d_radec.shape[1]) * (0.5/3600.)
    
        d2 = matched_good[:, 0:2] - matched_good[:, 3:5]

        extname = ota_list[ext].header['EXTNAME']

        if (False):
            fig = matplotlib.pyplot.figure()
            matplotlib.pyplot.plot(d_radec[:,0]*3600., d_radec[:,1]*3600., 'rs')
            matplotlib.pyplot.plot(d2[:,0]*3600., d2[:,1]*3600., 'bx')
            matplotlib.pyplot.xlim(-2,2)
            matplotlib.pyplot.ylim(-2,2)
    #        matplotlib.pyplot.xlim(-2e-3,2e-3)
    #        matplotlib.pyplot.ylim(-2e-3,2e-3)

    #        print radec
    #        print radec_computed
    #        print d_radec
    #        matplotlib.pyplot.plot(radec[:,0], radec[:,1], 'o', radec_computed[:,0], radec_computed[:,1], 'x')
            #matplotlib.pyplot.scatter(radec_computed[:,0], radec_computed[:,1])
            matplotlib.pyplot.grid(True)
            matplotlib.pyplot.xlabel("d_RA [arcsec]")
            matplotlib.pyplot.ylabel("d_DEC [arcsec]")
            #matplotlib.pyplot.show()
            fig.savefig(extname+".png")
            matplotlib.pyplot.close()

            fsdump_file = "source_%s.dmp" % (extname)
            wcsx = astWCS.WCS(ota_list[ext].header, mode="pyfits")
            sources_ota[:,0:2] = odi2wcs(sources_ota[:,2:4], wcsx)
            x = open(fsdump_file, "w")
            numpy.savetxt(x, sources_ota)
            x.close()

            fsdump_file = "source_%s.dmp" % (extname)
            wcsx = astWCS.WCS(ota_list[ext].header, mode="pyfits")
            sources_ota[:,0:2] = odi2wcs(sources_ota[:,2:4], wcsx)
            x = open(fsdump_file, "w")
            numpy.savetxt(x, sources_ota)
            x.close()

        
        #if (n_matches>10):
        #    print "Found enough matches to refine CDx_y"
        #    pinit.append(wcs.header['CD1_1'])
        #    pinit.append(wcs.header['CD1_2'])
        #    pinit.append(wcs.header['CD2_1'])
        #    pinit.append(wcs.header['CD2_2'])


        #print "Before fit", pinit
        #print "after fit", out[0]

        #if (n_matches>10):
        #    print "Found enough matches to refine CDx_y"
        #    ota_list[ext].header['CD1_1'] = out[0][2]
        #    ota_list[ext].header['CD1_2'] = out[0][3]
        #    ota_list[ext].header['CD2_1'] = out[0][4]
        #    ota_list[ext].header['CD2_2'] = out[0][5]

    #print matched_cat.shape, fixwcs_odi_sourcecat.shape


