import binascii
import logging
import numpy as np
import struct


_VALUE_FORMATTERS = {
    's': lambda x: "'{0}'".format(x),
}

class NiiHdrFieldInfo:
    def __init__(self, name, offset, dtype, shape=1):
        self.name = name
        self.offset = offset
        self.dtype = dtype
        self.shape = shape
        self.val = None
        self._formatter = _VALUE_FORMATTERS.get(self.dtype)
        if not self._formatter:
            self._formatter = lambda x:x

    def get_value(self):
        return self._formatter(self.val)

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
            if len(self.shape) > 1:
                #Transpose matrices from Fortran order
                self.val = arrdata.reshape(self.shape[::-1]).T
            else:
                self.val = arrdata.reshape(self.shape)
            logger.debug("Value: {0}".format(self.val))
        elif self.dtype == 's':
            #self.shape is the length of the string
            self.val = ""
            for i, c in enumerate(data):
                if i == self.shape:
                    break
                self.val += c
            logger.debug("Value: \"{0}\"".format(self.val))
        else:
            #logger.debug("Data: {0}".format(data))
            self.val = struct.unpack_from(self.dtype, data)[0]
            logger.debug("Value: {0}".format(self.val))
        return self._formatter(self.val)

    def __repr__(self):
        return "Field {0}[{1}] (val={2}, prec={3})".format(self.name,
            self.offset, self._formatter(self.val), self.dtype)
