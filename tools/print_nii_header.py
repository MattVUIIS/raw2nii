#!/usr/bin/env python
from __future__ import print_function
import argparse
import logging
import pprint

from nii_info import read_nii_header


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    options = parser.parse_args()
    is_many = len(options.filenames) > 1
    for filename in options.filenames:
        if is_many:
            print('{0}\n-----------------------------'.format(filename))
        header = read_nii_header(filename)
        pprint.pprint(header)
        if is_many:
            print('\n')
