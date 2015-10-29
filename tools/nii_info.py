import logging
import numpy as np

from nii_header_field_info import NiiHdrFieldInfo


HEADER = [
    NiiHdrFieldInfo('HdrSz', 0, 'i'),
    NiiHdrFieldInfo('DataType', 4, 's', 10),
    NiiHdrFieldInfo('db_name', 14, 's', 18),
    NiiHdrFieldInfo('extents', 32, 'i'),
    NiiHdrFieldInfo('session_error', 36, 'h'),
    NiiHdrFieldInfo('regular', 38, 'B'),
    NiiHdrFieldInfo('dim_info', 39, 'B'),
    NiiHdrFieldInfo('dim', 40, 'h', (8,)),
    NiiHdrFieldInfo('intent_p123', 56, 'f', (3,)),
    NiiHdrFieldInfo('intent_code', 68, 'h'),
    NiiHdrFieldInfo('datatype', 70, 'h'),
    NiiHdrFieldInfo('bitpix', 72, 'h'),
    NiiHdrFieldInfo('slice_start', 74, 'h'),
    NiiHdrFieldInfo('pixdim', 76, 'f', (8,)),
    NiiHdrFieldInfo('vox_offset', 108, 'f'),
    NiiHdrFieldInfo('scl_slope', 112, 'f'),
    NiiHdrFieldInfo('scl_inter', 116, 'f'),
    NiiHdrFieldInfo('slice_end', 120, 'h'),
    NiiHdrFieldInfo('slice_code', 122, 'B'),
    NiiHdrFieldInfo('xyzt_units', 123, 'B'),
    NiiHdrFieldInfo('cal_maxmin', 124, 'f', (2,)),
    NiiHdrFieldInfo('slice_duration', 132, 'f'),
    NiiHdrFieldInfo('toffset', 136, 'f'),
    NiiHdrFieldInfo('glmaxmin', 140, 'i', (2,)),
    NiiHdrFieldInfo('descrip', 148, 's', 80),
    NiiHdrFieldInfo('aux_file', 228, 's', 24),
    NiiHdrFieldInfo('qform_code', 252, 'h'),
    NiiHdrFieldInfo('sform_code', 254, 'h'),
    NiiHdrFieldInfo('quatern_bcd', 256, 'f', (3,)),
    NiiHdrFieldInfo('qoffset_xyz', 268, 'f', (3,)),
    NiiHdrFieldInfo('srow_xyz', 280, 'f', (4, 3)),
    NiiHdrFieldInfo('intent_name', 328, 's', 16),
    NiiHdrFieldInfo('magic', 344, 'i'),
]

def read_nii_header(filename):
    logger = logging.getLogger('raw2nii')
    with open(filename, 'rb') as f:
        header = _read_nii_header(f, logger)
    return header

def read_nii_body(filename):
    logger = logging.getLogger('raw2nii')
    with open(filename, 'rb') as f:
        header = _read_nii_header(f, logger)
        body = _read_nii_body(f, header, logger)
    return header, body

def read_nii(filename):
    logger = logging.getLogger('raw2nii')
    with open(filename, 'rb') as f:
        header = _read_nii_header(f, logger)
        body = _read_nii_body(f, header, logger)
    return header, body

def _read_nii_header(f, logger=None):
    header = {}
    header_bytes = f.read(348)
    if len(header_bytes) == 348:
        magic = header_bytes[-4:]
    else:
        magic = None
    if magic not in (bytes('n+1\0'), bytes('ni1\0')):
        if logger:
            logger.info('magic number: {0}'.format(magic))
        raise RuntimeError('This is not a NIFTI file!')
    for info in HEADER:
        header[info.name] = info.read_value(header_bytes)
    return header 

def _read_nii_body(f, header, logger):
    dim = header['dim'][1:5]
    bitpix = header['bitpix']
    datatype = header['datatype']
    unsigned_prefix = ('', 'u')[datatype ==
        {8: 2, 16: 512, 32: 768, 64: 1280}[bitpix]]
    dtype = '{0}int{1}'.format(unsigned_prefix, bitpix)
    data_shape = (dim[0], dim[1], dim[2] * dim[3])
    arr = np.fromstring(f.read(), dtype=dtype, count=np.product(data_shape))
    data = arr.reshape(data_shape[::-1]).T
    body = {'data': data, 'dtype': dtype}
    return body

