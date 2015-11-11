import logging
import numpy as np

import nifti_defines
import raw2nii_version
from calc_angulation import calc_angulation
from nifti_mat44_to_quatern import nifti_mat44_to_quatern
from nii_header import NiiHdr, NiiHdrField


#Reference for NIFTI header values can be found at:
#http://nifti.nimh.nih.gov/pub/dist/src/niftilib/nifti1.h

_HEADER_SIZE = 348
_FILLER_CHAR = ' '  # Used to pad strings
_DATATYPE_TABLE = {
    8: nifti_defines.kDT_UNSIGNED_CHAR,
    16: nifti_defines.kDT_SIGNED_SHORT,
    32: nifti_defines.kDT_FLOAT,
    64: nifti_defines.kDT_DOUBLE,
}

def CreateNiiHdr(par, angulation=None, dim=None):
    """
        create Nifti header from parameters as read from PAR file
    """
    logger = logging.getLogger('raw2nii')
    M, realvoxsize = calc_angulation(par, angulation)
    qoffset_xyz, quatern_bcd, qfac = nifti_mat44_to_quatern(M)
    NHdr = NiiHdr()
    NHdr.HdrSz = NiiHdrField(_HEADER_SIZE, 'i')
    NHdr.Data_Type = NiiHdrField(_FILLER_CHAR * 10, 's')
    NHdr.db_name = NiiHdrField(_FILLER_CHAR * 18, 's')
    NHdr.extents = NiiHdrField(0, 'i')
    NHdr.session_error = NiiHdrField(0, 'h')
    NHdr.regular = NiiHdrField(ord('r'), 'B')
    NHdr.dim_info = NiiHdrField(0, 'B')
    # Note on pixdim: use real voxel dimensions as calculated from
    # FOV/matrixsize in approp direction (CHECK!). Because for older Philips
    # releases, voxel dimensions in PAR file slice lines are rounded to 0.1!
    if np.nan is par.RT:
        par.RT = 1
    filedim = 4
    realdim = np.array([1, 1, 1])
    NHdr.dim = NiiHdrField(np.concatenate(([filedim], par.dim, dim,
        realdim)), 'h')
    NHdr.pixdim = NiiHdrField(np.concatenate(([qfac], realvoxsize,
        [par.RT, 1, 1, 1])), 'f')
    NHdr.intent_p123 = NiiHdrField(np.array([0, 0, 0]), 'f')
    NHdr.intent_code = NiiHdrField(0, 'h')
    if not par.issue:
        NHdr.datatype = NiiHdrField(_DATATYPE_TABLE[par.bit], 'h')
        NHdr.bitpix = NiiHdrField(par.bit, 'h')
    else:
        NHdr.datatype = NiiHdrField(nifti_defines.kDT_FLOAT, 'h')
        NHdr.bitpix = NiiHdrField(32, 'h')
    NHdr.slice_start = NiiHdrField(0, 'h')
    #vox_offset=352.0 means that the data starts immediately after the NIFTI-1
    #header
    NHdr.vox_offset = NiiHdrField(352, 'f')
    #Using a form of ternary statements here
    rs = {False: par.rescale_slope, True: 1}[par.issue]
    ri = {False: par.rescale_interc, True: 0}[par.issue]
    NHdr.scl_slope = NiiHdrField(rs, 'f')
    NHdr.scl_inter = NiiHdrField(ri, 'f')
    NHdr.slice_end = NiiHdrField(0, 'h')
    NHdr.slice_code = NiiHdrField(0, 'B')
    NHdr.xyzt_units = NiiHdrField(nifti_defines.kNIFTI_UNITS_MM
        | nifti_defines.kNIFTI_UNITS_SEC, 'B')
    NHdr.cal_maxmin = NiiHdrField(np.array([0, 0]), 'f')
    NHdr.slice_duration = NiiHdrField(0, 'f')
    NHdr.toffset = NiiHdrField(0, 'f')
    NHdr.glmaxmin = NiiHdrField(np.array([255, 0]), 'i')
    #Max 80 characters, left justify with spaces
    descrip = '{0}; converted by raw2nii {1}'.format(par.gen_info.protocol_name,
        raw2nii_version.VERSION)[:80].ljust(80)
    NHdr.descrip = NiiHdrField(descrip, 's')
    NHdr.aux_file = NiiHdrField(_FILLER_CHAR * 24, 's')
    #OVERRIDE VANDERBILT :
    # WARNING : the calcul for the qform is wrong ( calc quaternion not
    # working) only using sform, set qform to 0
    NHdr.qform_code = NiiHdrField(0, 'h')
    #NHdr.qform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    NHdr.sform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    NHdr.quatern_bcd = NiiHdrField(quatern_bcd, 'f')
    NHdr.qoffset_xyz = NiiHdrField(qoffset_xyz, 'f')
    NHdr.srow_xyz = NiiHdrField(M[0:3,:].T, 'f')
    NHdr.intent_name = NiiHdrField(_FILLER_CHAR * 16, 's')
    NHdr.magic = NiiHdrField(nifti_defines.kNIFTI_MAGIC_EMBEDDED_HDR, 'i')
    return NHdr
