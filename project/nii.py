from __future__ import division
import array
import binascii
import itertools
import logging
import numpy as np
import os
import struct

import nifti_defines
import raw2nii_version
from NiiFile import NiiHdr, NiiHdrField, HEADER_FIELD_NAMES


__all__ = ['NiiHdr', 'NiiHdrField', 'write_nii_from_par']

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
#Maps bit width -> unsigned datatype
_UDATATYPE_TABLE = {
    8: nifti_defines.kDT_UNSIGNED_CHAR,
    16: nifti_defines.kDT_UINT16,
    32: nifti_defines.kDT_UINT32,
    64: nifti_defines.kDT_UINT64,
}

def _create_nii_header(par):
    """ create Nifti header from parameters as read from PAR file """
    M, realvoxsize = _calc_angulation(par, True)
    qoffset_xyz, quatern_bcd, qfac = _nifti_mat44_to_quatern(M)
    hdr = NiiHdr()
    hdr.HdrSz = NiiHdrField(_HEADER_SIZE, 'i')
    hdr.Data_Type = NiiHdrField(_FILLER_CHAR * 10, 's')
    hdr.db_name = NiiHdrField(_FILLER_CHAR * 18, 's')
    hdr.extents = NiiHdrField(0, 'i')
    hdr.session_error = NiiHdrField(0, 'h')
    hdr.regular = NiiHdrField(ord('r'), 'B')
    hdr.dim_info = NiiHdrField(0, 'B')
    # Note on pixdim: use real voxel dimensions as calculated from
    # FOV/matrixsize in approp direction (CHECK!). Because for older Philips
    # releases, voxel dimensions in PAR file slice lines are rounded to 0.1!
    if np.nan is par.RT:
        par.RT = 1
    filedim = 4
    realdim = np.array([1, 1, 1])
    dim = np.array([par.nr_dyn * par.nr_diffgrads * par.nr_echos
        * max(par.nr_mrtypes, par.nr_realmrtypes)])
    hdr.dim = NiiHdrField(np.concatenate(([filedim], par.dim, dim, realdim)),
        'h')
    hdr.pixdim = NiiHdrField(np.concatenate(([qfac], realvoxsize,
        [par.RT, 1, 1, 1])), 'f')
    hdr.intent_p123 = NiiHdrField(np.array([0, 0, 0]), 'f')
    hdr.intent_code = NiiHdrField(0, 'h')
    hdr.multi_scaling_factors = par.multi_scaling_factors
    if not hdr.multi_scaling_factors:
        hdr.datatype = NiiHdrField(_DATATYPE_TABLE[par.bit], 'h')
        hdr.bitpix = NiiHdrField(par.bit, 'h')
    else:
        hdr.datatype = NiiHdrField(nifti_defines.kDT_FLOAT, 'h')
        hdr.bitpix = NiiHdrField(32, 'h')
    #Determine the type that should actually be written to the nifti
    if hdr.bitpix.val in (8, 16, 32, 64):
        if hdr.multi_scaling_factors and hdr.bitpix.val == 32:
            hdr.bitpixstr = 'float32'
        else:
            #Using a form of ternary statements here
            u_prefix = {True: 'u', False: ''}[
                hdr.datatype.val == _UDATATYPE_TABLE[hdr.bitpix.val]]
            hdr.bitpixstr = '{0}int{1}'.format(u_prefix, hdr.bitpix.val)
    hdr.slice_start = NiiHdrField(0, 'h')
    #vox_offset=352.0 means that the data starts immediately after the NIFTI-1
    #header
    hdr.vox_offset = NiiHdrField(352, 'f')
    #Using a form of ternary statements here
    rs = {False: par.rescale_slope, True: 1}[hdr.multi_scaling_factors]
    ri = {False: par.rescale_interc, True: 0}[hdr.multi_scaling_factors]
    hdr.scl_slope = NiiHdrField(rs, 'f')
    hdr.scl_inter = NiiHdrField(ri, 'f')
    hdr.slice_end = NiiHdrField(0, 'h')
    hdr.slice_code = NiiHdrField(0, 'B')
    hdr.xyzt_units = NiiHdrField(nifti_defines.kNIFTI_UNITS_MM
        | nifti_defines.kNIFTI_UNITS_SEC, 'B')
    hdr.cal_maxmin = NiiHdrField(np.array([0, 0]), 'f')
    hdr.slice_duration = NiiHdrField(0, 'f')
    hdr.toffset = NiiHdrField(0, 'f')
    hdr.glmaxmin = NiiHdrField(np.array([255, 0]), 'i')
    #Max 80 characters, left justify with spaces
    descrip = '{0}; converted by raw2nii {1}'.format(par.gen_info.protocol_name,
        raw2nii_version.VERSION)[:80].ljust(80)
    hdr.descrip = NiiHdrField(descrip, 's')
    hdr.aux_file = NiiHdrField(_FILLER_CHAR * 24, 's')
    #OVERRIDE VANDERBILT :
    # WARNING : the calcul for the qform is wrong ( calc quaternion not
    # working) only using sform, set qform to 0
    hdr.qform_code = NiiHdrField(0, 'h')
    #hdr.qform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    hdr.sform_code = NiiHdrField(nifti_defines.kNIFTI_XFORM_SCANNER_ANAT, 'h')
    hdr.quatern_bcd = NiiHdrField(quatern_bcd, 'f')
    hdr.qoffset_xyz = NiiHdrField(qoffset_xyz, 'f')
    hdr.srow_xyz = NiiHdrField(M[0:3,:].T, 'f')
    hdr.intent_name = NiiHdrField(_FILLER_CHAR * 16, 's')
    hdr.magic = NiiHdrField(nifti_defines.kNIFTI_MAGIC_EMBEDDED_HDR, 'i')
    return hdr

def _calc_angulation(par, angulation):
    if angulation:
        # trying to incorporate AP FH RL rotation angles: determined using some
        # common sense, Chris Rordon's help + source code and trial and error,
        # this is considered EXPERIMENTAL!
        rads = np.deg2rad(par.angRL)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r1 = np.array([[1, 0, 0], [0, cosrads, -sinrads], [0, sinrads, cosrads]])
        rads = np.deg2rad(par.angAP)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r2 = np.array([[cosrads, 0, sinrads], [0, 1, 0], [-sinrads, 0, cosrads]])
        rads = np.deg2rad(par.angFH)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r3 = np.array([[cosrads, -sinrads, 0], [sinrads, cosrads, 0], [0, 0, 1]])
        col = np.array([0, 0, 0, 1])[np.newaxis].T
        R_tot = np.concatenate((np.concatenate((r1.dot(r2).dot(r3),
            [np.zeros(3)])), col), axis=1)
    else:
        R_tot = np.eye(4)
    if 1 == par.sliceorient:  # Traversal
        lmm = np.eye(4)  # Do not rotate
        lXmm = par.fov_apfhrl[2] / par.dim[0]
        lYmm = par.fov_apfhrl[0] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[1] / par.dim[2]
    elif 2 == par.sliceorient:  # Sagittal
        lmm = np.array([[0, 0, -1, 0], [1, 0, 0, 0], [0, -1, 0, 0],
            [0, 0, 0, 1]])
        #Vanderbilt override
        lXmm = par.fov_apfhrl[0] / par.dim[0]
        lYmm = par.fov_apfhrl[1] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[2] / par.dim[2]
    elif 3 == par.sliceorient:  # Coronal
        lmm = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0],
            [0, 0, 0, 1]])  # Rotate 90 degrees
        #Vanderbilt override
        lXmm = par.fov_apfhrl[1] / par.dim[0]
        lYmm = par.fov_apfhrl[2] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[0] / par.dim[2]
    Zm = np.array([[lXmm, 0, 0, 0], [0, lYmm, 0, 0], [0, 0, lZmm, 0],
        [0, 0, 0, 1]])
    #realvoxsize is used to fill in pixdim nifti header info
    realvoxsize = np.array([lXmm, lYmm, lZmm])
    patient_to_tal = np.diag([-1, -1, 1, 1])
    analyze_to_dicom = np.diag([1, -1, 1, 1])
    A_tot = patient_to_tal.dot(R_tot).dot(Zm).dot(lmm).dot(analyze_to_dicom)
    p_orig = np.array([(par.dim[0] - 1) / 2, (par.dim[1] - 2) / 2,
        (par.dim[2] - 1) / 2, 1])
    offsetA = A_tot.dot(p_orig.T)
    if angulation:
        # trying to incorporate AP FH RL translation: determined using some
        # common sense, Chris Rordon's help + source code and trial and error,
        # this is considered EXPERIMENTAL!
        A_tot[0:3,3] = [-offsetA[0] - par.offRL,
            -offsetA[1] - par.offAP, -offsetA[2] + par.offFH]
    else:
        A_tot[0:3,3] = [-offsetA[0], -offsetA[1], -offsetA[2]]
    HdrMat = A_tot
    return HdrMat, realvoxsize

def _set_larger(A, B):
    C = A
    D = B
    if A > B:
        D = C
    else:
        C = D
    return C, D

def _nifti_mat44_to_quatern(A):
    A = A.copy()
    #offset outputs are read write out of input matrix
    qoffset_xyz = A[0:3,3]
    #compute lengths of each column; these determine grid spacings
    d1 = np.sqrt(np.sum(A[0:3,0:3]*A[0:3, 0:3], axis=0))
    #if a column length is zero, patch the trouble
    if d1[0] == 0:
        A[:,0] = [[1], [0], [0]]
        d1[0] = 1
    if d1[1] == 0:
        A[:,1] = [[0], [1], [0]]
        d1[1] = 1
    if d1[2] == 0:
        A[:,2] = [[0], [0], [1]]
        d1[2] = 1
    #normalize the columns
    A[0:3,0] = A[0:3,0] / d1[0]
    A[0:3,1] = A[0:3,1] / d1[1]
    A[0:3,2] = A[0:3,2] / d1[2]
    # At this point, the matrix has normal columns, but we have to allow
    # for the fact that the hideous user may not have given us a matrix
    # with orthogonal columns.
    # So, now find the orthogonal matrix closest to the current matrix.
    # One reason for using the polar decomposition to get this
    # orthogonal matrix, rather than just directly orthogonalizing
    # the columns, is so that inputting the inverse matrix to R
    # will result in the inverse orthogonal matrix at this point.
    # If we just orthogonalized the columns, this wouldn't necessarily hold.
    Q = A[0:3,0:3]
    U, S, V = np.linalg.svd(Q)
    #numpy's svd function transposes the V matrix compared to Matlab
    P = U.dot(V.T)
    #                            [ r11 r12 r13 ]
    # at this point, the matrix  [ r21 r22 r23 ] is orthogonal
    #                            [ r31 r32 r33 ]
    #
    # compute the determinant to determine if it is proper
    zd = np.linalg.det(P)
    if zd > 0:
        qfac = 1.0
    else: # improper ==> flip 3rd column
        qfac = -1.0
        P[:,2] = -P[:,2]
    #now, compute quaternion parameters
    a = np.trace(P) + 1
    if a > 0.5:  # simplest case
        a = 0.5 * np.sqrt(a)
        b = 0.25 * (P[2,1] - P[1,2]) / a
        c = 0.25 * (P[0,2] - P[2,0]) / a
        d = 0.25 * (P[1,0] - P[0,1]) / a
    else:  # trickier case
        xd = 1.0 + P[0,0] - (P[1,1] + P[2,2])
        yd = 1.0 + P[1,1] - (P[0,0] + P[2,2])
        zd = 1.0 + P[2,2] - (P[0,0] + P[1,1])
        if xd > 1.0:
            b = 0.5 * np.sqrt(xd)
            c = 0.25 * (P[0,1] + P[1,0]) / b
            d = 0.25 * (P[0,2] + P[2,0]) / b
            a = 0.25 * (P[2,1] - P[1,2]) / b
        elif yd > 1.0:
            c = 0.5 * np.sqrt(yd)
            b = 0.25 * (P[0,1] + P[1,0]) / c
            d = 0.25 * (P[1,2] + P[2,1]) / c
            a = 0.25 * (P[0,2] - P[2,0]) / c
        else:
            d = 0.5 * np.sqrt(zd)
            b = 0.25 * (P[0,2] + P[2,0]) / d
            c = 0.25 * (P[1,2] + P[2,1]) / d
            a = 0.25 * (P[1,0] - P[0,1]) / d
        if a < 0.0:
            b = -b
            c = -c
            d = -d
            a = -a
    quatern_bcd = np.array([b, c, d])
    return qoffset_xyz, quatern_bcd, qfac

def write_nii(fname, hdr, Data3D):
    """ Write the nifti to a file """
    try:
        with open(fname, 'wb') as fd:
            _write_nii_header(hdr, fd)  # write header to nii binary
            if hdr.bitpix.val in (8, 16, 32, 64):
                if hdr.multi_scaling_factors and hdr.bitpix.val == 32:
                    bitpixstr = 'float32'
                else:
                    #Using a form of ternary statements here
                    u_prefix = {True: 'u', False: ''}[
                        hdr.datatype.val == _UDATATYPE_TABLE[hdr.bitpix.val]]
                    bitpixstr = '{0}int{1}'.format(u_prefix, hdr.bitpix.val)
            #now add 4 extra bytes in space between header and offset for data
            #indicating that single .nii file ("n+1\0") rather than separate
            #img/hdr files were written. see http://nifti.nimh.nih.gov
            fd.write(bytearray(binascii.unhexlify('6e2b3100')))
            #add remaining 0s (probably not required)
            fd.write(bytearray([0] * (hdr.vox_offset.val - hdr.HdrSz.val
                - 4)))
            for s3 in range(Data3D.shape[2]):
                #flip data in order to write out
                Data3D[:,:,s3] = np.fliplr(Data3D[:,:,s3])
            #Transpose matrix before writing to get Fortran order
            Data3D.astype(bitpixstr).T.tofile(fd)
    except IOError as e:
        logger = logging.getLogger('raw2nii')
        logger.error('Write failed: {0}'.format(e))

def write_nii_from_par(nii_fname, par):
    """ Write the nifti to a file """
    logger = logging.getLogger('raw2nii')
    hdr = _create_nii_header(par)
    if par.dti_revertb0:
        _write_dynamics_files(nii_fname, par)
    try:
        logger.info('Writing file: {0}...'.format(nii_fname))
        with open(nii_fname, 'wb') as fd:
            _write_nii_header(hdr, fd)  # write header to nii binary
            #now add 4 extra bytes in space between header and offset for data
            #indicating that single .nii file ("n+1\0") rather than separate
            #img/hdr files were written. see http://nifti.nimh.nih.gov
            fd.write(bytearray(binascii.unhexlify('6e2b3100')))
            #add remaining 0s (probably not required)
            fd.write(bytearray([0] * (hdr.vox_offset.val - hdr.HdrSz.val - 4)))
            #Get the datatype to actually write the slices with
            if hdr.bitpix.val in (8, 16, 32, 64):
                if hdr.multi_scaling_factors and hdr.bitpix.val == 32:
                    bitpixstr = 'float32'
                else:
                    #Using a form of ternary statements here to prepend 'u' to
                    #bitpixstr if datatype is present in the _UDATATYPE_TABLE
                    u_prefix = {True: 'u', False: ''}[
                        hdr.datatype.val == _UDATATYPE_TABLE[hdr.bitpix.val]]
                    bitpixstr = '{0}int{1}'.format(u_prefix, hdr.bitpix.val)
            par_dt = {8: 'b', 16: 'h', 32: 'i'}[par.bit]
            #Read the REC file slice by slice and write to the nii right away.
            #This significantly reduces the memory required to process the REC.
            with open(par.rec_fname, 'rb') as rec:
                for slicenr, sl in enumerate(par.slices_sorted):
                    rec.seek(sl.index_in_rec_file * sl.recon_resolution_x *
                        sl.recon_resolution_y * 2)
                    sl_arr = array.array(par_dt)
                    sl_arr.fromfile(rec, sl.recon_resolution_x *
                        sl.recon_resolution_y)
                    sl_data = np.reshape(sl_arr, (sl.recon_resolution_x,
                        sl.recon_resolution_y)).T
                    if par.multi_scaling_factors:
                        sl_data = ((sl_data * sl.rescale_slope +
                            sl.rescale_intercept) / (sl.scale_slope *
                            sl.rescale_slope))
                    #Flip data left-to-right for radiological order then
                    #transpose matrix before writing to get Fortran order
                    np.fliplr(sl_data).astype(bitpixstr).T.tofile(fd)
        logger.info('  ...done')
    except IOError as e:
        logger = logging.getLogger('raw2nii')
        logger.error('Write failed: {0}'.format(e))
    return fd

def _write_nii_header(hdr, fd):
    logger = logging.getLogger('raw2nii')
    logger.debug('Writing NHdr...')
    for key in HEADER_FIELD_NAMES:
        hdrfield = getattr(hdr, key)
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

def _write_dynamics_files(nii_fname, par):
    logger = logging.getLogger('raw2nii')
    nslice = par.dim[2]
    if par.nr_diffgrads != par.NumberOfVolumes:
        logger.warning('VANDERBILT hack, the number of diffusion gradients is '
            'not coherent, taking the info from PAR header.')
        #Problem with slices_sorted that point the b0 at
        #Sorted the list to keep the b0 at the end
        #Get the first index for the b0
        r = par.slices_sorted.diffusion_b_value_number
        index = np.extract(r == 1, np.arange(r.shape[0]))
        #Put the info at the end
        B0 = np.copy(par.slices_sorted[index])
        par.slices_sorted[index] = par.slices_sorted[-nslice:]
        par.slices_sorted[-nslice:] = B0
        #par.slices_sorted = np.roll(par.slices_sorted, 1)  #This may work too
        #Keep only the number of volumes from par header
        par.nr_diffgrads = par.NumberOfVolumes
    logger.warning('Doing a very dirty hack to DTI data: putting b0 data as '
        'first in ND nii file.')
    #The last shall be first and the first shall be last
    par.slices_sorted = np.concatenate((par.slices_sorted[-nslice:],
        par.slices_sorted[:-nslice]))
    #Write the bval after sorting the slices. Put the b0 first if needed from
    #the par.slices_sorted in the bval and bvec
    Img_size = par.slices.shape[0]
    numberofslices_from_header = par.gen_info.max_number_of_slices_locations
    numberofslices = Img_size // par.NumberOfVolumes
    if not np.allclose(numberofslices_from_header, numberofslices):
        logger.warning('DTI incomplete. Number of slices different from header '
            'and reality.')
        sl_i = np.arange(0, par.NumberOfVolumes * numberofslices_from_header,
            numberofslices_from_header)
    else:  # No error with the number of slices
        sl_i = np.arange(0, par.NumberOfVolumes * numberofslices,
            numberofslices)
    #b0moved = (par.slices_sorted[-1].index_in_rec_file + 1 == Img_size)
    logger.warning('VANDERBILT hack -> putting last value in front for '
        'bval/bvec.')
    sl_i = np.roll(sl_i, 1)
    b_slices = par.slices[sl_i]
    name, ext = os.path.splitext(nii_fname)
    bval_filename = name + '-x-bval.txt'
    bvec_filename = name + '-x-bvec.txt'
    try:
        with open(bval_filename, 'wb') as bval_file:
            bval_file.writelines('{0:.6f} '.format(x)
                for x in b_slices.diffusion_b_factor)
    except OSError as e:
        logger.error('Failed to write bval text file "{0}": {1}'.format(
            bval_filename, e))
    try:
        with open(bvec_filename, 'wb') as bvec_file:
            bvec_file.writelines('{0:.6f} '.format(x)
                for x in b_slices.diffusion_rl)
            bvec_file.write('\n')
            #After checking the dtiqa process, need to flip Y data (so minus)
            bvec_file.writelines('{0:.6f} '.format(-y)
                for y in b_slices.diffusion_ap)
            bvec_file.write('\n')
            bvec_file.writelines('{0:.6f} '.format(z)
                for z in b_slices.diffusion_fh)
    except OSError as e:
        logger.error('Failed to write bvec text file "{0}": {1}'.format(
            bvec_filename, e))
