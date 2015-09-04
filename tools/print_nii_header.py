import argparse
import logging

from nii_info import read_nii_header


if __name__ == "__main__":
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.INFO)
    _formatter = logging.Formatter("%(levelname)s %(filename)s: %(message)s")
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--debug", "-d", action="store_true")
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    read_nii_header(options.filename)
