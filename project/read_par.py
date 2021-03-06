""" Philips PAR file interpreter. Reads several important fields from PAR
files, and returns them in a structure. Reads version number of Philips
research tools and interpretes accordingly. Research tools are used to
extract data from database; dataformats differ considerably between
versions. We now handle V3 and V4

function read_par
    par_fname: string with complete par-file name (with path)
    rec_fname: string with complete rec-file name (with path)
returns:
    par: A PARFile instance
"""
from __future__ import division
import logging
import numpy as np
import re

import par_defines
from PARFile import PARFile


__all__ = ['read_par']

#Maps image def values that have multiple fields -> names of those fields
_SUBVAR_NAMES = {
    'recon_resolution': ('x', 'y'),
    'image_angulation': ('ap', 'fh', 'rl'),
    'image_offcentre': ('ap', 'fh', 'rl'),
    'pixel_spacing': ('x', 'y'),
    'diffusion': ('ap', 'fh', 'rl'),
}

_IMAGE_INFORMATION_LINE = ('# === IMAGE INFORMATION ==========================='
    '===============================')

def read_par(par_fname, rec_fname):
    logger = logging.getLogger('raw2nii')
    par = PARFile()
    par.par_fname = par_fname
    par.rec_fname = rec_fname
    try:
        with open(par_fname, 'rb') as parfile:
            _skip_lines(parfile, 7)  # Skip first 7 lines
            par.version = parfile.readline().split()[-1]
            logger.debug('PAR version: {0}'.format(par.version))
            if 'V3' == par.version:
                raise NotImplementedError
            elif par.version in ('V4', 'V4.1', 'V4.2'):
                _skip_lines(parfile, 5)
                gen_info = _parse_general_info_V4X(par, parfile)
                logger.debug('Parameters name: {0}'.format(par.gen_info))
                _parse_definition_V4X(par, parfile)
                _skip_comment_lines(parfile)
                slices = _parse_slices_V4X(par, parfile)
    except OSError as e:
        par.problem_reading = True
        logger.error('Failed to read par file "{0}": {1}'.format(par_fname, e))
        return par
    if par.version in ('V4', 'V4.1', 'V4.2'):
        first_row = slices[0]
        last_row = slices[-1]
        #If there is more than one slice, check the order of the
        #slice numbers. 1 = ascending order, 2 = descending order
        if slices.shape[0] > 1:
            par.are_slices_sorted = (2, 1)[slices.slice_number[0] >
                slices.slice_number[1]]
        else:
            par.are_slices_sorted = True
        if gen_info.max_number_of_dynamics > 1:
            #estimate scan-duration from dtime PAR file row
            par.RT = (last_row.dyn_scan_begin_time -
                first_row.dyn_scan_begin_time) / (
                gen_info.max_number_of_dynamics - 1)
        else:
            par.RT = np.nan
        par.sliceorient = first_row.slice_orientation
        x = first_row.recon_resolution_x
        y = first_row.recon_resolution_y
        z = gen_info.max_number_of_slices_locations
        par.dim = np.array([x, y, z])
        par.multi_scaling_factors = (
            np.product(np.unique(par.slices.scale_slope).shape) != 1
            or np.product(np.unique(par.slices.rescale_intercept).shape) != 1
            or np.product(np.unique(par.slices.rescale_slope).shape) != 1)
        if par.multi_scaling_factors:
            logger.warning('Multiple scaling factors detected. Switching to '
                'float 32 nifti and rescaling')
            par.rescale_slope = slices.rescale_slope
            par.rescale_interc = slices.rescale_intercept
            par.scale_slope = slices.scale_slope
        else:
            par.rescale_slope = 1 / first_row.scale_slope
            par.rescale_interc = first_row.rescale_intercept
        par.bit = first_row.image_pixel_size
        par.slth = first_row.slice_thickness
        par.gap = first_row.slice_gap
        voxx = first_row.pixel_spacing_x
        voxy = first_row.pixel_spacing_y
        voxz = par.slth + par.gap
        par.vox = np.array([voxx, voxy, voxz])
        fovz, fovx, fovy = gen_info.fov
        par.fov = np.array([fovx, fovy, fovz])
        par.fov_apfhrl = np.array([fovz, fovx, fovy])
        par.angAP, par.angFH, par.angRL = gen_info.angulation_midslice
        par.offAP, par.offFH, par.offRL = gen_info.off_centre_midslice
    _check_number_of_volumes(par)
    _check_slice_orientation(par)
    _check_slice_order(par)
    _check_dti(par)
    logger.debug('PARFile {0}'.format(par))
    return par

def _parse_general_info_V4X(par, parfile):
    """ Reads the GENERAL INFORMATION section from the PAR file """
    line = None
    while line != '':
        pos = parfile.tell()
        line = parfile.readline()
        #Parse the useful parts of the general info entry on the left and
        #right of the colon: key and value
        m = re.search(r'\. ([^<>\(\)\[\]]*[a-zA-Z]).*: *(.*)', line)
        if not m:
            parfile.seek(pos)
            break
        key, val = m.group(1, 2)
        key = _sanitize_to_identifer(key).lower()
        #Try to guess the type of the field by conversion
        _val_split = val.split()
        if len(_val_split) > 1:
            try:
                val = np.array(tuple(float(x) for x in _val_split))
            except:
                pass
        else:
            try:
                val = int(val)
            except ValueError:
                pass
        #logger.debug("Key = '{0}' Val = '{1}'".format(key, val))
        setattr(par.gen_info, key, val)
    return par.gen_info

def _parse_definition_V4X(par, parfile):
    """ Reads the IMAGE INFORMATION DEFINITION section from the PAR file """
    line = None
    while line != '':
        pos = parfile.tell()
        line = parfile.readline().strip()
        #Parse the useful parts of the definition entry:
        #the identifier-valid name, the number of columns, and the type
        m = re.search(r'# ([^<>\(\)\[\]]*[a-zA-Z]).*\((\d+)?[\*]?(\w+)\)', line)
        if not m:
            if not par.fields:
                continue
            else:
                parfile.seek(pos)
                break
        var_descrip, type_len, type_descrip = m.group(1, 2, 3)
        var_name = _sanitize_to_identifer(var_descrip).lower()
        if type_len:
            type_len = int(type_len)
        else:
            type_len = 1
        #'string' should be interpreted as integer regardless
        if type_descrip == 'integer' or type_descrip == 'string':
            type_code = np.int64
        elif type_descrip == 'float':
            type_code = np.float64  # Same as MATLAB double
        else:
            raise ValueError(descrip)
        #Sub variables exist for variables that have size > 1
        #We add an underscore plus the name of the sub variable
        #i.e. image_angulation_x, image_angulation_y, image_angulation_z
        if type_len > 1:
            par.fields.extend(tuple((var_name + '_' + s, type_code)
                for s in _SUBVAR_NAMES[var_name]))
        else:
            par.fields.append((var_name, type_code))
    par.field_len = len(par.fields)
    return par.fields

def _parse_slices_V4X(par, parfile):
    """ Reads each slice line from the PAR file and calculates some metrics """
    logger = logging.getLogger('raw2nii')
    par.slices = np.loadtxt(parfile, dtype=par.fields, ndmin=1).view(
        np.recarray)
    if not par.slices.shape:
        logger.warning('par.slices has wrong shape: {0}, reshaping...'.format(
            par.slices.shape))
        par.slices = np.reshape(par.slices, (1,))
    if len(par.slices[0]) != par.field_len:
        raise ValueError('Slice tag format does not match the number of '
            'entries')
    #Determine number of interleaved image sequences (was:types,
    #name kept for historic reasons) (e.g. angio)
    par.nr_mrtypes = np.unique(par.slices.scanning_sequence).shape[0]
    #Determine number of interleaved echos
    par.nr_echos = np.unique(par.slices.echo_number).shape[0]
    #Determine number of interleaved image types (e.g. angio)
    par.nr_realmrtypes = np.unique(par.slices.image_type_mr).shape[0]
    #Determine number of diffusion gradients (e.g. DTI)
    par.nr_diffgrads = np.unique(par.slices.gradient_orientation_number
        ).shape[0]
    #Determine number of dynamics(directly from slice lines in
    #PAR file instead of PAR file header info!)
    par.nr_dyn = np.unique(par.slices.dynamic_scan_number).shape[0]
    if par.nr_dyn != par.gen_info.max_number_of_dynamics:
        logger.warning('Number of dynamics in header of PAR file does not '
            'match number of dynamics in the body')
    par.nr_bvalues = np.unique(par.slices.diffusion_b_value_number).shape[0]
    #Check if multishell
    par.is_multishell = par.nr_bvalues > 2
    #Sort the slices
    sort_order = (par.slices.slice_number, par.slices.dynamic_scan_number,
        par.slices.diffusion_b_value_number,
        par.slices.gradient_orientation_number, par.slices.echo_number,
        par.slices.image_type_mr, par.slices.scanning_sequence)
    if par.is_multishell:
        sort_order = sort_order[:2] + (sort_order[3], sort_order[2]) + (
            sort_order[4:])  # Swap diffusion b value and gradient orientation
    else:
        pass  # B0 and B1 diffusion weighting
    par.slices_sorted = par.slices[np.lexsort(sort_order)]
    return par.slices

def _check_number_of_volumes(par):
    logger = logging.getLogger('raw2nii')
    slices = par.slices
    par.NumberOfVolumes = np.unique(slices.dynamic_scan_number).shape[0]
    NoV = (slices.dynamic_scan_number.shape[0] //
        np.unique(slices.slice_number).shape[0])
    if NoV != par.NumberOfVolumes:
        logger.warning('Dynamic Scan Number does not match number of slices.'
            'Assuming slices are ordered.')
        #This code, transliterated from original override, is not quite right.
        #It also seems nonessential
        #cnt = np.zeros((max(slices.slice_number),))
        #for s in slices:
        #    cnt[s.slice_number - 1] += 1
        #    s.dynamic_scan_number = cnt[s.slice_number - 1]
        #par.NumberOfVolumes = np.unique(slices.dynamic_scan_number).shape[0]
        par.NumberOfVolumes = NoV

def _check_slice_orientation(par):
    logger = logging.getLogger('raw2nii')
    fov = par.gen_info.fov
    if par.sliceorient == par_defines.ORIENT_TRA:
        if not np.allclose(fov[0], fov[2]):
            logger.warning('AXIAL (TRA): par.gen_info.fov[0] != '
                'par.gen_info.fov[2]. Setting to max')
            fov[[0, 2]] = np.max(fov[[0, 2]])
        fov_div_slices = (fov[1] / par.gen_info.max_number_of_slices_locations)
    elif par.sliceorient == par_defines.ORIENT_COR:
        if not np.allclose(fov[1], fov[2]):
            logger.warning('AXIAL (COR): par.gen_info.fov[1] != '
                'par.gen_info.fov[2]. Setting to max')
            fov[[1, 2]] = np.max(fov[[1, 2]])
        fov_div_slices = (fov[0] / par.gen_info.max_number_of_slices_locations)
    elif par.sliceorient == par_defines.ORIENT_SAG:
        if not np.allclose(fov[1], fov[0]):
            logger.warning('AXIAL (SAG): par.gen_info.fov[1] != '
                'par.gen_info.fov[0]. Setting to max')
            fov[[0, 1]] = np.max(fov[[0, 1]])
        fov_div_slices = (fov[2] / par.gen_info.max_number_of_slices_locations)
    if not np.allclose(par.slth, fov_div_slices):
        logger.warning('Slice Thickness does not match fov/num_slices. '
            'ADJUSTING!!!')
        par.slth = fov_div_slices
    if not np.allclose(par.gap, 0):
        logger.warning('Non-zero slice gap: adjusting slice thickness')
        #Read the header variables
        par.slth += par.gap
        par.gap = 0

def _check_slice_order(par):
    """ Determine File Volume/Slice Order (inputVolumeSliceOrder):
        a) volume - all slices are listed (in order) for each volume before
        the next volume
        b) slices - the same slice is listed for all volumes before the next
        slice of the first volume (volumes are ordered)
        c) other - some other ordering (any ordering of volumes/slices is
        supported in the PAR file format)
        Procedure: Build a matrix with each row having: [VOLUME SLICE IDX] """
    logger = logging.getLogger('raw2nii')
    slices = par.slices
    N = slices.shape[0]
    SRT = np.concatenate((slices.dynamic_scan_number[np.newaxis].T,
        slices.slice_number[np.newaxis].T, np.arange(N)[np.newaxis].T),
        axis=1).reshape(N, 3)
    SRT = SRT[np.lexsort((SRT[:,0], SRT[:,1], SRT[:,2]))]
    is_in_volume_order = np.all(SRT[:,2] == np.arange(N))
    if is_in_volume_order:
        par.inputVolumeSliceOrder = 'volume'
    else:
        SRT = SRT[np.lexsort((SRT[:,1], SRT[:,0], SRT[:,2]))]
        is_in_slice_order = np.all(SRT[:,2] == np.arange(N))
        if is_in_slice_order:
            par.inputVolumeSliceOrder = 'slice'
        else:
            par.inputVolumeSliceOrder = 'unknown'
            logger.warning('Slice ordering is not a predefined type.')
            logger.info('This toolbox is compatible with arbitrary '
                'ordering of slices in PAR/REC files.')
            logger.info('However, other toolboxes or REC readers may '
                'assume a specific ordering.')

def _check_dti(par):
    par.dti_revertb0 = np.allclose(par.gen_info.diffusion, 1)  # True if dti

def _sanitize_to_identifer(name):
    """ Changes a string with various characters into an identifier that is
        valid for python. i.e. 'minimum RR-interval' -> 'minimum_RR_interval'
        """
    n = name.strip()
    n = re.sub('/', ' ', n)
    n = re.sub('-', ' ', n)
    n = re.sub(' +', '_', n)
    n = re.sub('[\W]+', '', n)
    return n

def _skip_lines(parfile, n):
    for i in range(n):
        parfile.readline()

def _skip_comment_lines(parfile):
    line = None
    while line != '':
        pos = parfile.tell()
        line = parfile.readline()
        if not line.startswith('#') and not line.strip() == '':
            parfile.seek(pos)
            break
