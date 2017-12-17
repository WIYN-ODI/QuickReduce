import logging
import pyfits
import numpy
import podi_logging

def read_fits_catalog(fn, extension, flatten=True):
    logger = logging.getLogger("ReadFITScat")

    logger.debug("Opening FITS catalog from %s" % (fn))
    if (type(fn) is str):
        hdulist = pyfits.open(fn)
    else:
        hdulist = fn

    if (type(extension) is list):
        ext_list = extension
    else:
        ext_list = [extension]

    return_tables = []

    for ext_id in ext_list:

        if (not ext_id in hdulist):
            logger.warning("Extension %s not found" % (ext_id))
            return_tables.append(None)
            continue

        ext = hdulist[ext_id]

        n_fields = ext.header['TFIELDS']
        n_rows = ext.header['NAXIS2']
        # table = numpy.empty((n_rows, n_fields))
        logger.debug("Reading data for %d fields" % (n_fields))
        # print "Reading data for %d fields" % (n_fields)
        table = []
        for f in range(n_fields):
            fd = ext.data.field(f)
            if (fd.ndim > 2 and flatten):
                logger.warning("Unable to handle catalog field %s" % (
                    ext.header['TTYPE%d' % (f + 1)]))
                continue
            if (fd.ndim > 1 and flatten):
                for c2 in range(fd.shape[1]):
                    table.append(fd[:, c2])
            else:
                table.append(fd)

        try:
            if (flatten):
                table = numpy.array(table).T
                logger.debug("Table data: %s" % (str(table.shape)))
            else:
                logger.debug("Table data: %d columns" % (len(table)))
            # np_table = numpy.array(table)
            # print table
        except:
            podi_logging.log_exception()
            pass

        # print table.shape

        # print table[1:3]

        # print ext.header
        return_tables.append(table)

    if (type(extension) is not list):
        return return_tables[0]

    return return_tables

