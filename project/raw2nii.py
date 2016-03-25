#!/bin/env python
from __future__ import division
import argparse
import logging
import numpy as np
import os
import re
import sys

from nii import write_nii_from_par
from write_parrec_from_dicom import write_parrec_from_dicom
from read_dicom import read_dicom
from read_par import read_par
from spm_type import spm_type


def raw_convert(input_file, output_file, **options):
    logger = logging.getLogger('raw2nii')
    #Convert from DCM to PAR
    if re.search('.dcm$', input_file, re.I) and re.search('.par$', output_file,
            re.I):
        return convert_dcm2par(input_file, output_file, **options)
    #Convert from PAR to NII
    if re.search('.par$', input_file, re.I) and re.search('.nii$', output_file,
            re.I):
        return convert_par2nii(input_file, output_file, **options)
    #Error
    logger.error('Conversion not supported')

def _get_rec_fname(par_fname):
    rec_fname, ext = os.path.splitext(par_fname)
    if '.par' == ext:
        rec_fname += '.rec'
    elif '.PAR' == ext:
        rec_fname += '.REC'
    return rec_fname

def convert_par2nii(par_fname, nii_fname, no_angulation, no_rescale,
        dti_revertb0):
    """
        no_angulation   : when True: do NOT include affine transformation as defined in PAR
                       file in hdr part of Nifti file (nifti only, EXPERIMENTAL!)
        no_rescale      : when True: do NOT store intensity scale as found in PAR
                       file (assumed equall for all slices). do NOT yield DV values.
        dti_revertb0 : when False (default), philips ordering is used for DTI data
                       (eg b0 image last). When True, b0 is saved as first image
                       in 3D or 4D data
    """
    logger = logging.getLogger('raw2nii')
    rec_fname = _get_rec_fname(par_fname)
    #extract the bval and bvec from the PAR file
    par = read_par(par_fname, rec_fname)
    if par.problem_reading:
        logger.warning('Skipping volume {0} because of reading errors.'
            .format(par_fname))
        return 1
    if 'V3' == par.version:
        raise NotImplementedError
    elif par.version in ('V4', 'V4.1', 'V4.2'):
        #new: loop slices (as in slice_index) and open and close
        #files along the way (according to info in index on dynamic
        #and mr_type)
        if not no_angulation:
            logger.warning('Assuming angulation parameters are identical '
                'for all scans in (4D) volume!')
        if not no_rescale:
            logger.warning('Assuming rescaling parameters (see PAR-file) '
                'are identical for all slices in volume and all scans in '
                '(4D) volume!')
        write_nii_from_par(nii_fname, par)
    else:
        logger.warning('Sorry, but data format extracted using Philips '
            'Research File format {0} was not known at the time the '
            'raw2nii software was developed'.format(par.version))
    return 0

def convert_dcm2par(dcm_fname, par_fname, **options):
    logger = logging.getLogger('raw2nii')
    dcm = read_dicom(dcm_fname)
    rec_fname = _get_rec_fname(par_fname)
    write_parrec_from_dicom(par_fname, rec_fname, dcm)
    return 0

def _generate_filename(sl):
    if nr_dyn > 1:
        dyn_suffix = '-{0:04d}'.format(
            sl.dynamic_scan_number)
        dyn_ndsuffix = '-d{0:04d}'.format(nr_dyn)
    else:
        dyn_suffix = '-{0:04d}'.format(1)
        dyn_ndsuffix = ''
    if nr_mrtypes > 1:
        mrtype_suffix = '-s{0:3d}'.format(
            sl.scanning_sequence)
        mrtype_ndsuffix = '-s{0:03d}'.format(nr_mrtypes)
    else:
        mrtype_suffix = ''
        mrtype_ndsuffix = ''
    if nr_realmrtypes > 1:
        realmrtype_suffix = '-t{0:03d}'.format(
            sl.image_type_mr)
        realmrtype_ndsuffix = '-t{0:03d}'.format(
            nr_realmrtypes)
    else:
        realmrtype_suffix = ''
        realmrtype_ndsuffix = ''
    if nr_echos > 1:
        echo_suffix = '-e{0:03d}'.format(
            sl.echo_number)
        echo_ndsuffix = '-e{0:03d}'.format(nr_echos)
    else:
        echo_suffix = ''
        echo_ndsuffix = ''
    if nr_diffgrads > 1:
        diffgrad_suffix = '-g{0:03d}'.format(
            sl.gradient_orientation_number)
        diffgrad_ndsuffix = '-g{0:03d}'.format(
            nr_diffgrads)
    else:
        diffgrad_suffix = ''
        diffgrad_ndsuffix = ''
    if nr_bvalues > 1:
        bval_suffix = '-b{0:03d}'.format(
            sl.diffusion_b_value_number)
        bval_ndsuffix = '-b{0:03d}'.format(nr_bvalues)
    else:
        bval_suffix = ''
        bval_ndsuffix = ''

def main():
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.INFO)
    _formatter = logging.Formatter('%(levelname)s %(asctime)s %(filename)s: '
        '%(message)s')
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('--no-rescale', action='store_false')
    parser.add_argument('--no-angulation', action='store_false')
    parser.add_argument('--dti_revertb0', action='store_true')
    parser.add_argument('input_file', type=str)
    parser.add_argument('output_file', type=str)
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    options = vars(options)
    options.pop('debug', None)
    sys.exit(raw_convert(**options))

if __name__ == '__main__':
    main()
