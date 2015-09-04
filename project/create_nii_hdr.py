import logging
import numpy as np

import nifti_defines
import raw2nii_version
from calc_angulation import calc_angulation
from nifti_mat44_to_quatern import nifti_mat44_to_quatern
from nii_header import NiiHdr, NiiHdrField


FILLER_CHAR = '\0'  # Used to pad strings

#           Struct precision values reference
#                                          Standard size
# format  C type              Python type  (in bytes)
#---------------------------------------------------------
# x       pad byte	          no value
# c	      char	              string of    1
#                               length 1
# b	      signed char	      integer	   1
# B       unsigned char	      integer	   1
# ?	      _Bool	              bool	       1
# h       short	              integer	   2
# H       unsigned short      integer	   2
# i       int	              integer	   4
# I       unsigned int	      integer	   4
# l       long	              integer	   4
# L	      unsigned long	      integer	   4
# q       long long	          integer	   8
# Q       unsigned long long  integer	   8
# f       float	              float	       4
# d       double	          float	       8
# s       char[]	          string
# p       char[]	          string
# P       void *	          integer

def CreateNiiHdr(par, angulation=None, rescale=None, dim=None):
    """
        create Nifti header from parameters as read from PAR file
    """
    logger = logging.getLogger('raw2nii')
    M, realvoxsize = calc_angulation(par, angulation)
    qoffset_xyz, quatern_bcd, qfac = nifti_mat44_to_quatern(M)
    NHdr = NiiHdr()
    NHdr.HdrSz = NiiHdrField(348, 'i')
    NHdr.Data_Type = NiiHdrField(FILLER_CHAR * 10, 's')
    NHdr.db_name = NiiHdrField(FILLER_CHAR * 18, 's')
    NHdr.extents = NiiHdrField(0, 'i')
    NHdr.session_error = NiiHdrField(0, 'h')
    NHdr.regular = NiiHdrField(114, 'B')
    NHdr.dim_info = NiiHdrField(0, 'B')
    # Note on pixdim: use real voxel dimensions as calculated from
    # FOV/matrixsize in approp direction (CHECK!). Because for older Philips
    # releases, voxel dimensions in PAR file slice lines are rounded to 0.1!
    if 3 == dim[0]:
        NHdr.dim = NiiHdrField(np.concatenate(([dim[0]], par.dim,
            [1, 1, 1, 1])), 'h')
        if np.nan is par.RT:
            par.RT = 1
        NHdr.pixdim = NiiHdrField(np.concatenate(([qfac], realvoxsize,
            [par.RT, 0, 0, 0])), 'f')
    elif 4 == dim[0]:
        NHdr.dim = NiiHdrField(np.concatenate(([dim[0]], par.dim,
            [dim[1], 1, 1, 1])), 'h')
        NHdr.pixdim = NiiHdrField(np.concatenate(([qfac], realvoxsize,
            [par.RT, 0, 0, 0])), 'f')
    elif 5 == dim[0]:
        #Find real dimension
        if dim[2] == 1:
            filedim = 4
            realdim = np.array([dim[2], 1, 1, 1])
        else:
            filedim = 5
            realdim = np.array([dim[2], 1, 1])
        NHdr.dim = NiiHdrField(np.concatenate(([filedim], par.dim, [dim[1]],
            realdim)), 'h')
        NHdr.pixdim = NiiHdrField(np.concatenate(([qfac], realvoxsize,
            [par.RT, 0, 0, 0])), 'f')
    NHdr.intent_p123 = NiiHdrField(np.array([0, 0, 0]), 'f')
    NHdr.intent_code = NiiHdrField(0, 'h')
    if 8 == par.bit:
        dt = 2
    elif 16 == par.bit:
        dt = 4
    elif 32 == par.bit:
        dt = 16
    elif 64 == par.bit:
        dt = 64
    NHdr.datatype = NiiHdrField(dt, 'h')
    NHdr.bitpix = NiiHdrField(par.bit, 'h')
    NHdr.slice_start = NiiHdrField(0, 'h')
    NHdr.vox_offset = NiiHdrField(352, 'f')
    if rescale == 1:
        rs = par.rescale_slope
        ri = par.rescale_interc
    else:
        rs = 1
        ri = 0
    NHdr.scl_slope = NiiHdrField(rs, 'f')
    NHdr.scl_inter = NiiHdrField(ri, 'f')
    NHdr.slice_end = NiiHdrField(0, 'h')
    NHdr.slice_code = NiiHdrField(0, 'B')
    NHdr.xyzt_units = NiiHdrField(10, 'B')
    NHdr.cal_maxmin = NiiHdrField(np.array([0, 0]), 'f')
    NHdr.slice_duration = NiiHdrField(0, 'f')
    NHdr.toffset = NiiHdrField(0, 'f')
    NHdr.glmaxmin = NiiHdrField(np.array([255, 0]), 'i')
    #Max 80 characters, left justify with spaces
    descrip = "{0}; converted by raw2nii {1}".format(par.name,
        raw2nii_version.VERSION)[:80].ljust(80)
    NHdr.descrip = NiiHdrField(descrip, 's')
    NHdr.aux_file = NiiHdrField(FILLER_CHAR * 24, 's')
    NHdr.qform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    NHdr.sform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    NHdr.quatern_bcd = NiiHdrField(quatern_bcd, 'f')
    NHdr.qoffset_xyz = NiiHdrField(qoffset_xyz, 'f')
    NHdr.srow_xyz = NiiHdrField(M[0:3,:].T, 'f')
    NHdr.intent_name = NiiHdrField(FILLER_CHAR * 16, 's')
    NHdr.magic = NiiHdrField(nifti_defines.kNIFTI_MAGIC_EMBEDDED_HDR, 'i')
    return NHdr

if __name__ == "__main__":
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.DEBUG)
    CreateNiiHdr()
