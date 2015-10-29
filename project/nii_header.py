import binascii
import logging
import numpy as np
import struct


class NiiHdr:
    FIELD_NAMES = ('HdrSz', 'Data_Type', 'db_name', 'extents', 'session_error',
        'regular', 'dim_info', 'dim', 'intent_p123', 'intent_code',
        'datatype', 'bitpix', 'slice_start', 'pixdim', 'vox_offset',
        'scl_slope', 'scl_inter', 'slice_end', 'slice_code', 'xyzt_units',
        'cal_maxmin', 'slice_duration', 'toffset', 'glmaxmin', 'descrip',
        'aux_file', 'qform_code', 'sform_code', 'quatern_bcd',
        'qoffset_xyz', 'srow_xyz', 'intent_name', 'magic')

    def __init__(self, **kwargs):
        for name in self.FIELD_NAMES:
            setattr(self, name, kwargs.get(name))

    def write_to_file(self, fd):
        for key in self.FIELD_NAMES:
            hdrfield = getattr(self, key)
            prec = hdrfield.prec
            val = hdrfield.val
            if not isinstance(val, np.ndarray):  # If writing simple value
                if prec == 's':  # Writing a string
                    packed = val
                    #logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                    #    "packed string='{4}')".format(key, fd.tell(), prec,
                    #    val, packed))
                    fd.write(packed)
                else:
                    packed = struct.pack(prec, val)
                    #logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                    #    "packed bytes='{4}')".format(key, fd.tell(), prec,
                    #    val, binascii.hexlify(packed)))
                    fd.write(packed)
            else: # If writing array
                if len(val.shape) == 1:  # If 1D array
                    packed = struct.pack(prec * len(val), *list(val))
                    #logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                    #    "packed array='{4}')".format(key, fd.tell(), prec,
                    #    val, binascii.hexlify(packed)))
                    fd.write(packed)
                else:
                    #Transpose matrix before writing to get Fortran order
                    packed = val.astype(prec).T
                    #logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                    #    "packed matrix='{4}')".format(key, fd.tell(), prec,
                    #    val, binascii.hexlify(packed.tobytes())))
                    packed.tofile(fd)

    def __repr__(self):
        s = []
        for key in self.fieldnames():
            val = getattr(self, key)
            s.append('{0}="{1}"'.format(key, val))
        return "<NiiHdr {0}>".format(" ".join(s))


class NiiHdrField:
    """ Arbitrary structure """
    def __init__(self, val, prec):
        self.val = val
        self.prec = prec

    def __repr__(self):
        return "<struct val={0} prec={1}>".format(self.val, self.prec)
