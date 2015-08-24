import argparse
import logging

from nii_header_field_info import NiiHdrFieldInfo


HEADER = [
    NiiHdrFieldInfo('HdrSz', 0, 'i'),
    NiiHdrFieldInfo('DataType', 4, 's'),
    NiiHdrFieldInfo('db_name', 14, 's'),
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
    NiiHdrFieldInfo('descrip', 148, 's'),
    NiiHdrFieldInfo('aux_file', 228, 's'),
    NiiHdrFieldInfo('qform_code', 252, 'h'),
    NiiHdrFieldInfo('sform_code', 254, 'h'),
    NiiHdrFieldInfo('quatern_bcd', 256, 'f', (3,)),
    NiiHdrFieldInfo('qoffset_xyz', 268, 'f', (3,)),
    NiiHdrFieldInfo('srow_xyz', 280, 'f', (4, 3)),
    NiiHdrFieldInfo('intent_name', 328, 's'),
    NiiHdrFieldInfo('magic', 344, 'i'),
]

def print_nii_header(filename=None):
    logger = logging.getLogger('raw2nii')

    with open(filename, 'rb') as f:
        header_bytes = f.read(348)
        for info in HEADER:
            info.read_value(header_bytes)
            logger.info(info)


if __name__ == "__main__":
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.INFO)
    _formatter = logging.Formatter("%(levelname)s %(filename)s: %(message)s")
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--debug", "-d", action="store_true")
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    print_nii_header(options.filename)
