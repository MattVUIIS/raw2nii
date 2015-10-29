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


def main(test_data_folder, test_folder, rtol, atol):
    """
        Runs convert_raw2nii on a folder containing PARREC test data and
        compares the NIFTI output
    """
    logger = logging.getLogger('raw2nii_test')
    test_data_folder = os.path.abspath(test_data_folder)
    path_par = os.path.join(test_data_folder, test_folder, "PARREC")
    input_par = glob.glob("{0}/*.PAR".format(path_par))[0]
    out_folder = os.path.join("testOutput", test_folder, "NIFTI")
    #Set up convert_raw2nii logger
    raw2nii_logger = logging.getLogger('raw2nii')
    #raw2nii_logger.setLevel(logging.DEBUG)
    _formatter = logging.Formatter('test_data: %(levelname)s: %(message)s')
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    raw2nii_logger.addHandler(_stream_handler)
    #Run NIFTI conversion function and output to testOutput directory
    output = convert_raw2nii(filelist=[input_par], prefix="", suffix="",
        pathpar=path_par, outfolder=out_folder, outputformat=1, angulation=1,
        rescale=1, dti_revertb0=0)
    #Compare to the 4D data since we have the Vanderbilt override behavior
    test_nifti_folder = os.path.join(test_data_folder, test_folder, "NIFTI",
        "4D")
    #Produce listing of each file in NIFTI folders in lexographical order and
    #prepare each pair for comparison
    test_niftis = sorted(glob.glob("{0}/*.nii".format(test_nifti_folder)))
    output_niftis = sorted(glob.glob("{0}/*.nii".format(out_folder)))
    if len(test_niftis) != len(output_niftis):
        logger.error('Different number of output files! Test niftis: {0} '
            'Output niftis: {1}'.format(pprint.pformat(test_niftis),
            pprint.pformat(output_niftis)))
    output_pairs = zip(test_niftis, output_niftis)
    for test_nifti, output_nifti in output_pairs:
        logger.info("Output pair: {0} -> {1}".format(test_nifti, output_nifti))
        test_header, test_body = read_nii(test_nifti)
        output_header, output_body = read_nii(output_nifti)
        #First compare the headers of the two niftis
        for test_key, test_val in test_header.items():
            output_val = output_header[test_key]
            #logger.debug("Compare header.{0} : {1} -> {2}".format(
            #    test_key, test_val, output_val))
            if isinstance(test_val, np.ndarray):
                is_different = not np.allclose(test_val, output_val, rtol, atol)
            elif isinstance(test_val, float):
                is_different = abs(test_val - output_val) > (atol + rtol *
                    abs(output_val))
            else:
                is_different = (test_val != output_val)
            if is_different:
                logger.warning("header.{0} are different: {1} -> {2}".format(
                    test_key, test_val, output_val))
        #Then compare the bodies of the niftis
        if test_body['dtype'] != output_body['dtype']:
            logger.warning("Body dtypes are different: {0} -> {1}".format(
                test_body['dtype'], output_body['dtype']))
        elif test_body['data'].shape != output_body['data'].shape:
            logger.warning("Body shapes are different: {0} -> {1}".format(
                test_body['dtype'].shape, output_body['dtype'].shape))
        elif not np.allclose(test_body['data'], output_body['data'], rtol,
                atol):
            diff_output, err_output = body_diff_script(logger, test_body,
                output_body)
            logger.warning("Body data is different: <{0}> <{1}>".format(
                diff_output, err_output))
    return 0

def body_diff_script(logger, test_body, output_body):
    #Write out test binary data to a tempfile
    test_bin_fd, test_bin_name = tempfile.mkstemp()
    with os.fdopen(test_bin_fd, 'wb') as f:
        test_body['data'].astype(test_body['dtype']).T.tofile(f)
    #Write out canonical binary data to a tempfile
    output_bin_fd, output_bin_name = tempfile.mkstemp()
    with os.fdopen(output_bin_fd, 'wb') as f:
        output_body['data'].astype(output_body['dtype']).T.tofile(f)
    #Compare two tempfiles
    return subprocess.Popen(['./tools/nii_body_diff.sh', test_bin_name,
        output_bin_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).communicate()

if __name__ == "__main__":
    logger = logging.getLogger('raw2nii_test')
    logger.setLevel(logging.INFO)
    _formatter = logging.Formatter("%(levelname)s: %(message)s")
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", "-d", action="store_true")
    parser.add_argument("test_data_folder")
    parser.add_argument("test_folder")
    parser.add_argument("--rel-tolerance", "-r", type=float, default=1e-05)
    parser.add_argument("--abs-tolerance", "-a", type=float, default=1e-08)
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    main(options.test_data_folder, options.test_folder, options.rel_tolerance,
        options.abs_tolerance)
