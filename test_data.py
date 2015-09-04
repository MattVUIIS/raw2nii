import argparse
import bsdiff4
import glob
import logging
import numpy as np
import os
import pprint
import subprocess

from tools.nii_info import read_nii


def main(test_data_folder, test_folder, rtol, atol):
    """
        Runs convert_raw2nii on a folder containing PARREC test data and
        compares the NIFTI output
    """
    logger = logging.getLogger('raw2nii_test')
    path_par = os.path.join(test_data_folder, test_folder, "PARREC")
    input_par = glob.glob("{0}/*.PAR".format(path_par))[0]
    out_folder = os.path.join("testOutput", test_folder, "NIFTI")
    #Run NIFTI conversion script and output to testOutput directory
    cmd = ["python", "convert_raw2nii.py", input_par, "--pathpar", path_par,
        "--outfolder", out_folder]
    logger.debug("Running command: {0}".format(cmd))
    os.chdir("project")
    status = subprocess.call(cmd)
    if status != 0:
        logger.error("Error occurred")
        return 1
    test_nifti_folder = os.path.join(test_data_folder, test_folder, "NIFTI",
        "3D")
    #Produce listing of each file in NIFTI folders in lexographical order and
    #prepare each pair for comparison
    test_niftis = sorted(glob.glob("{0}/*.nii".format(test_nifti_folder)))
    output_niftis = sorted(glob.glob("{0}/*.nii".format(out_folder)))
    output_pairs = zip(test_niftis, output_niftis)
    #print("Test niftis: {0}".format(pprint.pformat(test_niftis)))
    #print("Output niftis: {0}".format(pprint.pformat(output_niftis)))
    #print("Output pairs: {0}".format(pprint.pformat(output_pairs)))
    for test_nifti, output_nifti in output_pairs:
        logger.info("Output pair: {0} -> {1}".format(test_nifti, output_nifti))
        test_header, test_body = read_nii(test_nifti)
        output_header, output_body = read_nii(output_nifti)
        #logger.debug("Test header: {0}".format(pprint.pformat(test_header)))
        #logger.debug("Test body: {0}".format(len(test_body)))
        #logger.debug("Output header: {0}".format(pprint.pformat(output_header)))
        #logger.debug("Output body: {0}".format(len(output_body)))
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
        if test_body != output_body:
            diff_output = bsdiff4.diff(test_body, output_body)
            logger.warning("Bodies are different: {0}".format(diff_output))
    return 0

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
