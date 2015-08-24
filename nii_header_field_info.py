import binascii
import logging
import numpy as np
import struct


class NiiHdrFieldInfo:
    def __init__(self, name, offset, dtype, shape=1):
        self.name = name
        self.offset = offset
        self.dtype = dtype
        self.shape = shape
        self.val = None

    def read_value(self, header_bytes):
        logger = logging.getLogger('raw2nii')
        data = header_bytes[self.offset:]
        logger.debug("Reading value for: {0}[{1}]{{{2}*{3}}}".format(self.name,
            self.offset, self.dtype, self.shape))
        logger.debug("Data: {{{0}}}({1})".format(data, binascii.hexlify(data)))
        if isinstance(self.shape, tuple):
            count = reduce(lambda x, y: x * y, self.shape)
            arrdata = np.fromstring(data, dtype=self.dtype, count=count)
            logger.debug("Count: {0} Arrdata: {1}".format(count, arrdata))
            self.val = arrdata.reshape(self.shape)
            if len(self.shape) > 1:
                self.val = self.val.T
            logger.debug("Value: {0}".format(self.val))
        elif self.dtype == 's':
            self.val = ""
            for c in data:
                if ord(c) == 0:
                    #Leave off the last null byte of the C string
                    break
                self.val += c
            logger.debug("Value: \"{0}\"".format(self.val))
        else:
            #logger.debug("Data: {0}".format(data))
            self.val = struct.unpack_from(self.dtype, data)
            logger.debug("Value: {0}".format(self.val))
        return self.val

    def __repr__(self):
        return "Field {0}[{1}] (prec={2}, val='{3}')".format(self.name,
            self.offset, self.dtype, self.val)
