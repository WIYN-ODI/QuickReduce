#!/usr/bin/env python

import os, sys
import numpy
import pyfits
import scipy
import scipy.spatial
import scipy.interpolate
import bottleneck
import math
import multiprocessing
import Queue
import ctypes

import astLib.astWCS as astWCS

from podi_commandline import *
import podi_associations
import podi_logging
from podi_definitions import *

class PhotFlatFrame(object):

    def __init__(self, filename):

        self.logger = logging.getLogger("PhotFlatFrame")
        self.nmax = 250

        # load the FITS and extract the photcalib table
        self.filename = filename

        hdulist = pyfits.open(self.filename)

        # load and store all WCS structures
        self.is_valid = False
        self.wcs = {}
        self.ref_header = None
        self.ota_list = []
        self.extname_list = []
        for ext in hdulist:
            if (not is_image_extension(ext)):
                continue
            ota = ext.header['OTA']
            self.logger.debug("reading wcs for %s in %s" % (ext.name, self.filename))
            #self.wcs[ext.name] = astWCS.WCS(ext.header, mode='pyfits')
            self.wcs[ota] = astWCS.WCS(ext.header, mode='pyfits')
            if (self.ref_header is None):
                self.ref_header = ext.header
            print "Reading:", self.filename, ext.header['OTA'], ext.name
            self.ota_list.append(ext.header['OTA'])
            self.extname_list.append(ext.name)

        #
        # now read the photcalib table
        #
        self.field_names = {}
        try:
            photcalib_tbhdu = hdulist['CAT.PHOTCALIB']
        except:
            self.logger.warning("No PHOTCALIB table found in %s" % (self.filename))
            photcalib_tbhu = None
            raise Exception("no PHOTCALIB extension")

        #
        # Now extract all sources for each OTA
        #
        tfields = photcalib_tbhdu.header['TFIELDS']
        nstars = photcalib_tbhdu.header['NAXIS2']
        catalog = numpy.empty((nstars, tfields))

        for field in range(tfields):
            catalog[:, field] = photcalib_tbhdu.data.field(field)
            ttype = 'TTYPE%d' % (field + 1)
            if (ttype in photcalib_tbhdu.header):
                fieldname = photcalib_tbhdu.header[ttype]
                self.field_names[fieldname] = field

        self.catalog = catalog
        self.logger.info("Read %d matched sources from %s" % (
            catalog.shape[0], self.filename))

        _, bn = os.path.split(self.filename)
        dbg_catfile = "cat_%s.cat" % (bn[:-5])
        logger.debug("Saving debug catalog to %s" % (dbg_catfile))
        numpy.savetxt(dbg_catfile, self.catalog)

        # now compute all zeropoints
        mag0size = hdulist[0].header['MAG0SIZE']
        photfilt = hdulist[0].header['PHOTFILT'].upper()
        self.logger.debug("PHOTOMETRIC FILTER: %s" % (photfilt))
        if (photfilt.lower() not in ['u', 'g', 'r', 'i', 'z']):
            self.logger.critical("Unable to identify PHOTFILT value: %s" % (photfilt))
            return None

        odimag = self.catalog[:, self.field_names['ODI_MAG_D%d' % (int(mag0size*10))]]
        refmag = self.catalog[:, self.field_names['REF_%s' % (photfilt)]]
        odierr = self.catalog[:, self.field_names['ODI_ERR_D%d' % (int(mag0size*10))]]
        referr = self.catalog[:, self.field_names['REF_ERR_%s' % (photfilt)]]
        self.zeropoint = refmag - odimag
        self.zeropoint_error = numpy.hypot(odierr, referr)
        self.logger.info("computed %d zeropoints" % (odimag.shape[0]))

        # also store some info about the zeropoints already computed
        self.ref_ra = self.ref_header['CRVAL1']
        self.ref_dec = self.ref_header['CRVAL2']
        self.cos_dec = numpy.cos(numpy.radians(self.ref_dec))
        print("cos-dec=", self.cos_dec)

        self.deprojected_coords = numpy.array(self.catalog[:, 0:2])
        self.deprojected_coords[:,0] *= self.cos_dec
        self.coord_tree = scipy.spatial.cKDTree(self.deprojected_coords)
        self.is_valid = True
        #print self.field_names

    def valid(self):
        return self.is_valid

    def get_source_indices(self, ra, dec, radius=1, nmax=None, relative_coords=False):
        _ra = ra * self.cos_dec
        k = self.nmax if nmax is None else nmax

        if (relative_coords):
            _ra = self.ref_ra*self.cos_dec - ra/60
            dec = self.ref_dec - dec/60

        self.logger.debug("Searching for %d sources within %f arcmin of %f, %f" % (k, radius, _ra, dec))
        # query which sources match
        d, i = self.coord_tree.query(
            [_ra,dec],
            distance_upper_bound=radius/60.,
            k=k,
            p=2) #, n_jobs=-1)
        valid_indices = numpy.isfinite(d)
        idx = i[valid_indices]

        return idx

    def get_zeropoints(self, ra, dec, radius=1, nmax=None, relative_coords=False, max_error=None):
        """

        :param ra: right ascension of search cone
        :param dec: declination of search cone
        :param radius: radius of search cone, in arcmin
        :return: list of zeropoints as numpy.ndarray

        """

        idx = self.get_source_indices(ra, dec, radius, nmax, relative_coords)
        print idx
        zp = self.zeropoint[idx]
        zperr = self.zeropoint_error[idx]

        if (max_error is not None):
            good_error = zperr < max_error
            zp = zp[good_error]
            zperr = zperr[good_error]

        return zp

    def get_ota_list(self):
        return self.ota_list

    def get_extname_list(self):
        return self.extname_list


    def get_ota_zeropoints(self, ota=None, x=2048, y=2048, radius=2048, strict_ota=False, return_error=True):
        if (ota is None):
            return None

        #print ota

        #select_ota = (self.catalog[:, self.field_names['OTA']] == ota)
        #ota_cat = self.catalog[select_ota]

        # convert pixel position to ra/dec
        radec = self.wcs[ota].pix2wcs(x,y)
        #print radec
        #return radec

        indices = self.get_source_indices(radec[0], radec[1],
                                          radius=radius*0.11/60)
        return_cat = self.catalog[indices]
        zp = self.zeropoint[indices]
        zperr = self.zeropoint_error[indices]

        if (strict_ota):
            # only return sources from the selected ota
            ota_num = ota #int(ota[3:5])
            this_ota = (return_cat[:, self.field_names['OTA']] == ota_num)
            return_cat = return_cat[this_ota]
            zp = zp[this_ota]
            zperr = zperr[this_ota]

        #print return_cat.shape,
        if (return_error):
            return zp, zperr

        return zp





# class PhotFrame(object):
#     def __init__(self, filename):
#
#         self.field_names = {}
#         self.catalog = None
#         self.wcs = {}
#         self.file_loaded = False
#         self.zeropoints_computed = False
#         self.apertures = numpy.array([20, 30, 40, 50, 60, 80, 100, 120])
#
#         self.filename = filename
#         _, bn = os.path.split(filename)
#         self.filebase = bn
#
#         self.logger = logging.getLogger("PF(%s)" % (bn))
#
#         self.neighbor_count = 250
#         self.neighbor_radius = 0
#
#         self.scaling_factor = 0
#
#         self.read_frame()
#
#     def read_frame(self):
#
#         if (self.file_loaded):
#             return self.file_loaded
#
#         self.logger.debug("Reading %s" % (self.filename))
#         hdulist = pyfits.open(self.filename)
#
#         #
#         # Read general relavant properties
#         #
#
#         self.phot_reference = hdulist[0].header['PHOTMCAT']
#         # self.magzero = hdulist[0].header['MAGZERO']
#         self.magzero = hdulist[0].header['PHOTZP_X']
#         self.ref_filter = hdulist[0].header['PHOTFILT']
#
#         #
#         # Read the photometric catalog(s)
#         #
#
#         try:
#             photcalib_tbhdu = hdulist['CAT.PHOTCALIB']
#         except:
#             self.logger.warning("No PHOTCALIB table found in %s" % (self.filename))
#             photcalib_tbhu = None
#             return False
#
#         #
#         # Now extract all sources for each OTA
#         #
#         tfields = photcalib_tbhdu.header['TFIELDS']
#         nstars = photcalib_tbhdu.header['NAXIS2']
#         catalog = numpy.empty((nstars, tfields))
#
#         for field in range(tfields):
#             catalog[:, field] = photcalib_tbhdu.data.field(field)
#             ttype = 'TTYPE%d' % (field + 1)
#             if (ttype in photcalib_tbhdu.header):
#                 fieldname = photcalib_tbhdu.header[ttype]
#                 self.field_names[fieldname] = field
#
#         self.catalog = catalog
#         self.logger.info("Read %d matched sources from %s" % (
#             catalog.shape[0], self.filename))
#
#         #
#         # Read all WCS headers
#         #
#         for ext in hdulist:
#             if (not is_image_extension(ext)):
#                 continue
#             wcs = astWCS.WCS(ext.header, mode='pyfits')
#             extname = ext.name
#
#             self.wcs[extname] = wcs
#
#         #
#         # compute a simplified WCS solution, ignoring distortion and the fact
#         # there are multiple OTAs
#         #
#         self.setup_quick_wcs()
#
#         #
#         # Update the neighbor radius
#         #
#         self.neighbor_radius = 4. * 60. / self.pixelscale
#
#         hdulist.close()
#         self.file_loaded = True
#         self.logger.debug("File read completed!")
#         return True
#
#     def setup_quick_wcs(self):
#         ota44 = self.wcs['OTA44.SCI'].header
#         simple = pyfits.ImageHDU()
#         for key in ['CD1_1', 'CD2_2', 'CD1_2', 'CD2_1',
#                     'CRVAL1', 'CRVAL2',
#                     'CRPIX1', 'CRPIX2',
#                     'NAXIS', 'NAXIS1', 'NAXIS2']:
#             simple.header[key] = ota44[key]
#         self.simple_wcs = self.wcs['OTA44.SCI']  # astWCS.WCS(simple.header, mode='pyfits')
#
#         self.simple_coords = numpy.array(self.simple_wcs.wcs2pix(self.catalog[:, 0], self.catalog[:, 1]))
#
#         self.simple_tree = scipy.spatial.cKDTree(self.simple_coords)
#         self.pixelscale = self.simple_wcs.getPixelSizeDeg() * 3600.
#
#         # numpy.savetxt("dummy_%s.cat" % (self.filebase),
#         #               numpy.append(self.simple_coords, self.catalog, axis=1))
#         # numpy.savetxt("dummy2_%s.cat" % (self.filebase), self.simple_coords)
#
#     def compute_zeropoints(self):
#
#         self.logger.debug("computing zeropoints")
#
#         ref_mag = "SDSS_MAG_%s" % (self.ref_filter.upper())
#         ref_err = "SDSS_ERR_%s" % (self.ref_filter.upper())
#
#         self.zeropoints = numpy.empty((self.catalog.shape[0], self.apertures.shape[0]))
#         self.zeropoints[:, :] = numpy.NaN
#
#         for idx, ap in enumerate(self.apertures):
#
#             odi_mag = "ODI_MAG_D%d" % (ap)
#             odi_err = "ODI_ERR_D%d" % (ap)
#
#             if (odi_mag in self.field_names and
#                         ref_mag in self.field_names):
#                 zp = self.catalog[:, self.field_names[ref_mag]] - \
#                      self.catalog[:, self.field_names[odi_mag]]
#                 self.zeropoints[:, idx] = zp
#
#         self.zeropoints_computed = True
#
#     def convert_to_relative(self):
#
#         # # compute position in the bottom left corner of OTA44
#         # ra_dec = self.wcs['OTA44.SCI'].pix2wcs(0,0)
#         # print ra_dec
#
#         #
#         # select all objects within 4 arcmin of the assumed center.
#         # define the median zeropoint as reference zeropoint
#         #
#         d, i = self.simple_tree.query(
#             [0, 0],
#             p=2,
#             k=1000,  # use 1000 sources at most
#             distance_upper_bound=(4 * 60 / self.pixelscale),
#         )
#         good_match = numpy.isfinite(d) & (i < self.catalog.shape[0])
#
#         near_center = self.zeropoints[i[good_match]]
#         # numpy.savetxt("center.cat", near_center)
#
#         center_zp = bottleneck.nanmedian(near_center, axis=0)
#
#         # self.zeropoints_relative = self.zeropoints - center_zp
#         self.zeropoints_relative = self.zeropoints - self.magzero
#
#         # print center_zp.shape
#         # print center_zp
#
#     def prep(self):
#         self.read_frame()
#         self.compute_zeropoints()
#         self.convert_to_relative()
#
#     def get_correction(self, ra_dec):
#
#         # make sure Ra/Dec has the right dimensions
#         if (ra_dec.ndim == 1):
#             ra_dec = ra_dec.reshape((1, -1))
#
#         #
#         # convert Ra/Dec to X/Y in the simple projected image
#         #
#         xy = numpy.array(self.simple_wcs.wcs2pix(ra_dec[:, 0], ra_dec[:, 1]))
#
#         #
#         # Query all stars around the given coordinates
#         #
#         d, i = self.simple_tree.query(xy, p=2,
#                                       k=self.neighbor_count,
#                                       distance_upper_bound=self.neighbor_radius,
#                                       )
#         # print d.shape
#
#         valid = numpy.isfinite(d) & (i < self.catalog.shape[0])
#
#         # prepare the result buffer
#         all_corrections = numpy.empty((ra_dec.shape[0], self.zeropoints_relative.shape[1]))
#         all_corrections[:, :] = numpy.NaN
#
#         for idx in range(ra_dec.shape[0]):
#             # print idx, ra_dec[idx], xy[idx]
#
#             nearby_sources = i[idx, :][valid[idx, :]]
#             # this_valid = valid[idx, :]
#             # print this_valid.shape, numpy.sum(this_valid), nearby_sources.shape
#
#             rel_zp = self.zeropoints_relative[nearby_sources]
#             # print rel_zp.shape
#
#             all_corrections[idx, :] = numpy.median(rel_zp, axis=0)
#             # print correction.shape
#
#         return all_corrections
#
#     def get_weight(self, targetzp=25.):
#         scaling_factor = math.pow(10, 0.4 * (targetzp - self.magzero))
#         return 1. / scaling_factor
#
#
#     def get_zeropoint(self, ra, dec):
#         pass
#



class PhotFlatHandler(object):

    def __init__(self, filelist):
        self.logger = logging.getLogger("PhotFlat")

        self.filelist = filelist

        self.phot_frames = {}

        self.ota_from_extname = {}
        self.extname_from_ota = {}
        pass

    def read_catalogs(self):

        self.logger.info("Reading PHOTCALIB catalogs")
        for idx, fn in enumerate(self.filelist):

            if (not os.path.isfile(fn)):
                self.logger.warning("File %s does not exist" % (fn))
                continue

            new_frame = PhotFlatFrame(fn)
            if (new_frame.is_valid):
                self.phot_frames[fn] = new_frame

                for i in range(len(new_frame.get_ota_list())):#
                    ota = new_frame.get_ota_list()[i]
                    extname = new_frame.get_extname_list()[i]
                    # ,extname in new_frame.get_ota_list(),new_frame.get_extname_list():
                    self.ota_from_extname[extname] = ota
                    self.extname_from_ota[ota] = extname

    def extname2ota(self, extname):
        return self.ota_from_extname[extname]

    def ota2extname(self, ota):
        return self.extname_from_ota[ota]





def expand_to_fullres_worker(job_queue, photflat, blocksize, shmem_out, shmem_shape, memlock):

    # Prepare the relative coordinates
    x, y = numpy.indices((blocksize, blocksize), dtype=numpy.float)
    y /= blocksize
    x /= blocksize
    # print x[:5,:5], x[-5:,-5:]

    omx = 1. - x
    omy = 1. - y

    omx_omy = omx * omy
    omx_y = omx * y
    x_omy = x * omy
    x_y = x * y

    # out_buffer = numpy.zeros((4096,4096))
    # out_buffer[:,:] = numpy.NaN

    # print ix,iy

    out_buffer = shmem_as_ndarray(shmem_out).reshape(shmem_shape)
    out = numpy.empty((blocksize, blocksize))

    while (True):

        try:
            job = job_queue.get_nowait()
        except Queue.Empty:
            # print "done!"
            break

        (ix, iy) = job
        #print ix,iy

        #
        # Follow the algorithm outlined in
        # https://en.wikipedia.org/wiki/Bilinear_interpolation#Unit_Square
        #
        try:
            f_00 = photflat[ix, iy]
            f_01 = photflat[ix, iy+1]
            f_10 = photflat[ix+1, iy]
            f_11 = photflat[ix+1, iy+1]
            # f_00 = photflat[iy, ix]
            # f_01 = photflat[iy, ix + 1]
            # f_10 = photflat[iy + 1, ix]
            # f_11 = photflat[iy + 1, ix + 1]
        except IndexError:
            print "Index error accessing", ix, iy
            continue

        if (ix>=15 or iy>=15):
            print ix,iy, f_00, f_01, f_10, f_11

        out = f_00 * omx_omy + f_10 * omx_y + f_01 * x_omy + f_11 * x_y
        # out = f_00 * omx_omy + f_10 * x_omy + f_01 * omx_y + f_11 * x_y
        #out[:,:] = f_00


        memlock.acquire()
        #print "writing line",iy
        out_buffer[iy * blocksize:(iy + 1) * blocksize,
        ix * blocksize:(ix + 1) * blocksize] = out
        memlock.release()




def expand_to_fullres(photflat, blocksize, out_dimension=None):

    if (out_dimension is None):
        out_dimension = (4096, 4096)

    #
    # Prepare parallel interpolation upwards to full-resolution
    #
    _x,_y = out_dimension
    out_shmem = multiprocessing.RawArray(ctypes.c_float, _x*_y)
    out_buffer = shmem_as_ndarray(out_shmem).reshape(out_dimension)
    out_buffer[:, :] = numpy.NaN
    job_queue = multiprocessing.JoinableQueue()
    data_lock = multiprocessing.Lock()

    # prepare all jobs - each job interpolates one little block
    for ix, iy in itertools.product(range(photflat.shape[0]-1), repeat=2):
        job_queue.put((ix,iy))

    processes = []
    for i in range(7):
         p = multiprocessing.Process(target=expand_to_fullres_worker,
                                     kwargs={
                                         'job_queue': job_queue,
                                         'shmem_out': out_shmem,
                                         'shmem_shape': out_buffer.shape,
                                         'photflat': photflat, #correction_2d,
                                         'blocksize': blocksize,
                                         'memlock': data_lock,
                                     }
         )
         p.start()
         processes.append(p)

    for p in processes:
        p.join()


    # Now we have the photometric flat-field, convert it to scaling frame
    #print("writing results")
    photflat_2d = numpy.power(10, 0.4*out_buffer)

    return photflat_2d


def create_photometric_flatfield(
        filelist,
        strict_ota=False,
        resolution_arcsec=60.,
        debug=False,
):

    logger = logging.getLogger("PhotFlat")

    pf = PhotFlatHandler(filelist=filelist)
    logger.info("Input files:\n-- %s" % ("\n-- ".join(filelist)))

    pf.read_catalogs()

    reference_pos = [4., -4.]
    # that's in arc-min relative to reference point from CRVAL1/2

    reference_zp = {}
    list_of_otas = []
    list_of_extnames = []
    for framename in pf.phot_frames:
        logger.info("Adding photometric data from file %s" % (framename))
        frame = pf.phot_frames[framename]

        zps = frame.get_zeropoints(ra=reference_pos[0],
                                   dec=reference_pos[1],
                                   radius=3,
                                   relative_coords=True, max_error=0.05)
        print zps
        reference_zp[framename] = numpy.median(zps)

        # also collect a list of all available OTAs
        list_of_otas.extend(frame.get_ota_list())
        list_of_extnames.extend(frame.get_extname_list())

    print reference_zp

    unique_otas = set(list_of_otas)
    print list_of_otas

    unique_extnames = set(list_of_extnames)
    print unique_extnames

    #
    # Now extract the relative ZP differences for each of the sectors in each ota
    #
    pixel_resolution = 512
    enlarge = 2.
    n_samples = int(math.ceil(4096. / pixel_resolution))+1
    sample_pixels = int(math.floor((4096. / (n_samples-1))))

    logger.info("Using pixel-grid of %d^2 samples every %d pixels" % (n_samples, sample_pixels))
    ref_points = numpy.arange(n_samples) * sample_pixels
    print n_samples, sample_pixels, n_samples*sample_pixels

    iy,ix = numpy.indices((n_samples, n_samples)) * sample_pixels

    running_sum = 0
    otalist = [pyfits.PrimaryHDU()]

    for i, extname in enumerate(unique_extnames):

        logger.info("Computing photometric flat-field for OTA %s (%2d of %2d)" % (extname, i+1, len(unique_otas)))
        ota = pf.extname2ota(extname)

        photflat = numpy.empty((n_samples,n_samples))
        photflat_err = numpy.empty((n_samples,n_samples))
        dump = open("dump_%d" % (ota), "w")
        for _x,_y in itertools.product(range(n_samples), repeat=2):
            #print x,y

            x = ref_points[_x]
            y = ref_points[_y]

            full_zp_list = None
            full_zperr_list = None
            for framename in pf.phot_frames:
                frame = pf.phot_frames[framename]
                #print framename, ota, x, y, frame
                zp_list, zperr_list = frame.get_ota_zeropoints(
                    ota=ota, x=x, y=y, radius=pixel_resolution*enlarge,
                    strict_ota=strict_ota,
                    return_error=True)
                # correct zp for the specific offset
                zp_list = zp_list - reference_zp[framename]
                if (full_zp_list is None):
                    full_zp_list = zp_list
                    full_zperr_list = zperr_list
                else:
                    full_zp_list = numpy.append(full_zp_list, zp_list, axis=0)
                    full_zperr_list = numpy.append(full_zperr_list, zperr_list, axis=0)

            # print ota, ota[3:5], x, y, full_zp_list.shape[0]
            running_sum += full_zp_list.shape[0]

            # ideally, do some outlier rejection here
            _, mask = three_sigma_clip(full_zp_list, return_mask=True)
            full_zp_list = full_zp_list[mask]
            full_zperr_list = full_zperr_list[mask]

            # compute the weighted mean correction factor
            small_error = full_zperr_list < 0.03
            full_zperr_list = full_zperr_list[small_error]
            full_zp_list = full_zp_list[small_error]
            photflat[_x,_y] = numpy.sum(full_zp_list/full_zperr_list) / numpy.sum(1./full_zperr_list)
            photflat_err[_x,_y] = numpy.std(full_zp_list)

            numpy.savetxt(dump,
                          numpy.array([full_zp_list, full_zperr_list, numpy.ones((full_zp_list.shape[0]))*photflat[_x,_y]]).T)
                          #numpy.append(full_zp_list.reshape((-1,1)), full_zperr_list.reshape((-1,1)), axis=1))
            print >>dump, "\n"*10

            # photflat[_x,_y] = numpy.median(full_zp_list)

        combined = numpy.empty((n_samples**2, 3))
        combined[:,0] = ix.ravel()
        combined[:,1] = iy.ravel()
        combined[:,2] = photflat.ravel()
        numpy.savetxt("photflat.%02d.combined" % (ota), combined)
        numpy.savetxt("photflat.%02d.photflat" % (ota), photflat)

        combined[:,2] = photflat_err.ravel()
        numpy.savetxt("photflat.%02d.err" % (ota), combined)

        fullres = expand_to_fullres(photflat, blocksize=pixel_resolution)
        imghdu = pyfits.ImageHDU(data=fullres, name=extname)
        # add some headers to allow for mosaic viewing
        otax, otay = int(math.floor(ota/10)), int(ota%10)
        s = 4096
        detsec = "[%d:%d,%d:%d]" % (otax * s, (otax + 1) * s, otay * s, (otay + 1) * s)
        imghdu.header["DETSEC"] = (detsec, "position of OTA in focal plane")
        otalist.append(imghdu)

    logger.info("Total sum of reference values: %d" % (running_sum))
    # break

    hdulist = pyfits.HDUList(otalist)
    return hdulist


    #
    # Put all frames on the same relative photometry grid, using the zeropoint
    # from the vicinity around the reference point for absolute calibration
    #





    # src_cat_info = get_clean_cmdline()[2]
    # src_cat_fn = src_cat_info.split(",")[0]
    # logger.info("Reading source catalog (%s)" % (src_cat_fn))
    # src_cat = numpy.loadtxt(src_cat_fn)
    #
    # logger.info("Finding corrections for %d sources/positions" % (src_cat.shape[0]))
    # #pf.get_correction(numpy.array([[49.400,41.286],[49.35,41.20]]))
    # correction = pf.get_correction(src_cat[:,0:2], debug=debug)
    #
    # out_fn = "%s.corr" % (src_cat_fn)
    # numpy.savetxt(out_fn, correction)



if __name__ == "__main__":

    options = {}
    podi_logging.setup_logging(options)

    logger = logging.getLogger("PhotFlat")

    output_filename = sys.argv[1]
    filelist = get_clean_cmdline()[2:]
    debug = cmdline_arg_isset("-debug")
    strict_ota = cmdline_arg_isset("-ota")
    resolution = float(cmdline_arg_set_or_default("-resolution", 60.))

    hdulist = create_photometric_flatfield(
        filelist=filelist,
        debug=debug,
        strict_ota=strict_ota,
        resolution_arcsec=resolution,
    )
    clobberfile(output_filename)
    hdulist.writeto(output_filename, clobber=True)

    podi_logging.shutdown_logging(options)


