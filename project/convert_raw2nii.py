from __future__ import division
import argparse
import logging
import numpy as np
import os

from nii import write_nii_from_par
from read_par import read_par
from spm_type import spm_type


def convert_raw2nii(filelist, prefix, suffix, pathpar, outfolder, outputformat,
        angulation, rescale, dti_revertb0):
    """
        filelist: list containing PAR file names (without path info) to convert
            into analyze/nifti
        prefix       : characters to prepend to all output filenames.
                       Blank by default. The prefix, PAR filename (without the
                       extension), suffix, and the file number are used
                       respectively to form the output filename.
        suffix       : characters to append to all output filenames.
                       Blank by default. The prefix, PAR filename (without the
                       extension), suffix, and the file number are used
                       respectively to form the output filename.
        pathpar      : complete path containing PAR files (with trailing /)
        outfolder    : when set, files will be written to this directory,
                       including lowest level folder containing parfile
        outputformat : 1 for Nifty output format (spm5), 2 for Analyze (spm2)
        angulation   : when 1: include affine transformation as defined in PAR
                       file in hdr part of Nifti file (nifti only, EXPERIMENTAL!)
        rescale      : when 1: store intensity scale as found in PAR
                       file (assumed equall for all slices). Yields DV values.
        dti_revertb0 : when 0 (default), philips ordering is used for DTI data
                       (eg b0 image last). When 1, b0 is saved as first image
                       in 3D or 4D data
    """
    logger = logging.getLogger('raw2nii')
    pathpar = os.path.expanduser(pathpar)
    outfiles = []
    for nbf, par_fname in enumerate(filelist):
        full_par_fname = os.path.join(pathpar, par_fname)
        full_rec_fname, ext = os.path.splitext(full_par_fname)
        if '.par' == ext:
            full_rec_fname += '.rec'
        elif '.PAR' == ext:
            full_rec_fname += '.REC'
        #extract the bval and bvec from the PAR file
        par = read_par(full_par_fname, full_rec_fname)
        if par.problem_reading:
            logger.warning('Skipping volume {0} because of reading errors.'
                .format(full_par_fname))
            continue
        if outfolder:
            absDir = os.path.abspath(os.path.expanduser(outfolder))
        else:
            absDir = os.path.abspath('NIFTI')
        if not os.path.exists(absDir):
            logger.debug('Makedirs on: {0}'.format(absDir))
            os.makedirs(absDir)
        logger.debug('Check dir: {0}'.format(absDir))
        out_basename, _ = os.path.splitext(os.path.basename(par_fname))
        nii_fname = prefix + out_basename + suffix
        if 'V3' == par.version:
            raise NotImplementedError
        elif par.version in ('V4', 'V4.1', 'V4.2'):
            #new: loop slices (as in slice_index) and open and close
            #files along the way (according to info in index on dynamic
            #and mr_type)
            if angulation:
                logger.warning('Assuming angulation parameters are identical '
                    'for all scans in (4D) volume!')
            if rescale:
                logger.warning('Assuming rescaling parameters (see PAR-file) '
                    'are identical for all slices in volume and all scans in '
                    '(4D) volume!')
            volname = os.path.join(absDir, nii_fname) + '.nii'
            write_nii_from_par(volname, par)
        else:
            logger.warning('Sorry, but data format extracted using Philips '
                'Research File format {0} was not known at the time the '
                'raw2nii software was developed'.format(par.version))
    return outfiles

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
    parser.add_argument('filelist', nargs='+', type=str)
    parser.add_argument('--prefix', default='')
    parser.add_argument('--suffix', default='')
    parser.add_argument('--pathpar', default='')
    parser.add_argument('--outfolder', default=None)
    parser.add_argument('--outputformat', type=int, default=1)
    parser.add_argument('--angulation', type=int, default=1)
    parser.add_argument('--rescale', type=int, default=1)
    parser.add_argument('--dti_revertb0', type=int, default=0)
    parser.add_argument('--debug', '-d', action='store_true')
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    convert_raw2nii(options.filelist, options.prefix, options.suffix,
        options.pathpar, options.outfolder, options.outputformat,
        options.angulation, options.rescale, options.dti_revertb0)

if __name__ == '__main__':
    main()
