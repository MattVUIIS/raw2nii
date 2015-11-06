import argparse
import glob
import logging
import numpy as np
import os
import pprint
import subprocess
import tempfile

from project import convert_raw2nii
from tools.nii_info import read_nii


DEFAULT_TEST_DATA = [
    ('canonicalData', 'philips/philips_3T_achieva/fMRI/'),
    ('canonicalData', 'philips/philips_3T_achieva/dti/'),
    #('canonicalData', 'philips/philips_3T_achieva/b0map/'),
    #('canonicalData', 'philips/philips_1_5T_intera/volume_2/'),
]


def run_test(canon_data_folder, trial_folder, rtol, atol):
    """ Runs convert_raw2nii on a folder containing PARREC test data and
        compares the NIFTI output """
    logger = logging.getLogger('raw2nii_test')
    canon_data_folder = os.path.abspath(canon_data_folder)
    path_par = os.path.join(canon_data_folder, trial_folder, 'PARREC')
    input_par = glob.glob('{0}/*.PAR'.format(path_par))[0]
    test_folder = os.path.join('testOutput', trial_folder, 'NIFTI')
    #Run NIFTI conversion function and output to testOutput directory
    logger.info('Running test: python convert_raw2nii.py {0} --pathpar={1} '
        '--outfolder={2}'.format(input_par, path_par, test_folder))
    output = convert_raw2nii(filelist=[input_par], prefix='', suffix='',
        pathpar=path_par, outfolder=test_folder, outputformat=1, angulation=1,
        rescale=1, dti_revertb0=0)
    #Compare to the 4D data since we have the Vanderbilt override behavior
    canon_nifti_folder = os.path.join(canon_data_folder, trial_folder, 'NIFTI',
        '4D')
    #Produce listing of each file in NIFTI folders in lexographical order and
    #prepare each pair for comparison
    canon_niftis = sorted(glob.glob('{0}/*.nii'.format(canon_nifti_folder)))
    test_niftis = sorted(glob.glob('{0}/*.nii'.format(test_folder)))
    if len(canon_niftis) != len(test_niftis):
        logger.error('Different number of output files! Test niftis: {0} '
            'Output niftis: {1}'.format(pprint.pformat(canon_niftis),
            pprint.pformat(test_niftis)))
    file_pairs = zip(canon_niftis, test_niftis)
    for canon_nifti, test_nifti in file_pairs:
        logger.info('File pair: {0} -> {1}'.format(canon_nifti, test_nifti))
        canon_header, canon_body = read_nii(canon_nifti)
        test_header, test_body = read_nii(test_nifti)
        #First compare the headers of the two niftis
        for canon_key, canon_val in canon_header.items():
            test_val = test_header[canon_key]
            #logger.debug('Compare header.{0} : {1} -> {2}'.format(canon_key,
            #   canon_val, test_val))
            if isinstance(canon_val, np.ndarray):
                is_different = not np.allclose(canon_val, test_val, rtol, atol)
            elif isinstance(canon_val, float):
                is_different = abs(canon_val - test_val) > (atol + rtol *
                    abs(test_val))
            else:
                is_different = (canon_val != test_val)
            if is_different:
                logger.warning('header.{0} are different: {1} -> {2}'.format(
                    canon_key, canon_val, test_val))
        #Then compare the bodies of the niftis
        if canon_body['dtype'] != test_body['dtype']:
            logger.warning('Body dtypes are different: {0} -> {1}'.format(
                canon_body['dtype'], test_body['dtype']))
        elif canon_body['data'].shape != test_body['data'].shape:
            logger.warning('Body shapes are different: {0} -> {1}'.format(
                canon_body['dtype'].shape, test_body['dtype'].shape))
        elif not np.allclose(canon_body['data'], test_body['data'], rtol, atol):
            diff_output, err_output = subprocess.Popen([
                './tools/nii_body_diff.sh', canon_nifti, test_nifti],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            logger.warning('Body data is different:\n{0}{1}'.format(
                diff_output, err_output))
    return 0

if __name__ == '__main__':
    logger = logging.getLogger('raw2nii_test')
    logger.setLevel(logging.INFO)
    _formatter = logging.Formatter('%(levelname)s: %(message)s')
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', '-d', action='store_true')
    parser.add_argument('canon_data_folder', nargs='?')
    parser.add_argument('trial_folder', nargs='?')
    parser.add_argument('--rel-tolerance', '-r', type=float, default=1e-05)
    parser.add_argument('--abs-tolerance', '-a', type=float, default=1e-08)
    options = parser.parse_args()
    if options.canon_data_folder and options.trial_folder:
        test_args = [(options.canon_data_folder, options.trial_folder)]
    else:
        test_args = DEFAULT_TEST_DATA
    if options.debug:
        logger.setLevel(logging.DEBUG)
    #Set up convert_raw2nii logger
    raw2nii_logger = logging.getLogger('raw2nii')
    #raw2nii_logger.setLevel(logging.DEBUG)
    raw2nii_logger.addHandler(_stream_handler)
    for canon_data_folder, trial_folder in test_args:
        run_test(canon_data_folder, trial_folder, options.rel_tolerance,
            options.abs_tolerance)
