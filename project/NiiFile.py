#Reference for NIFTI header values can be found at:
#http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h
HEADER_FIELD_NAMES = ('HdrSz', 'Data_Type', 'db_name', 'extents',
    'session_error', 'regular', 'dim_info', 'dim', 'intent_p123', 'intent_code',
    'datatype', 'bitpix', 'slice_start', 'pixdim', 'vox_offset', 'scl_slope',
    'scl_inter', 'slice_end', 'slice_code', 'xyzt_units', 'cal_maxmin',
    'slice_duration', 'toffset', 'glmaxmin', 'descrip', 'aux_file',
    'qform_code', 'sform_code', 'quatern_bcd', 'qoffset_xyz', 'srow_xyz',
    'intent_name', 'magic')

class NiiHdr:
    def __repr__(self):
        s = []
        for key in self.fieldnames():
            val = getattr(self, key)
            s.append('{0}="{1}"'.format(key, val))
        return '<NiiHdr {0}>'.format(' '.join(s))

class NiiHdrField:
    """ Arbitrary structure """
    def __init__(self, val, prec):
        self.val = val
        self.prec = prec

    def __repr__(self):
        return '<struct val={0} prec={1}>'.format(self.val, self.prec)
