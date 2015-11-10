#!/usr/bin/env python
from __future__ import print_function
import argparse
import logging
import numpy as np
import pprint

from nii_info import read_nii


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    options = parser.parse_args()
    is_many = len(options.filenames) > 1
    for filename in options.filenames:
        if is_many:
            print('{0}\n-----------------------------'.format(filename))
        header, body = read_nii(filename)
        pprint.pformat('dtype: {0}'.format(body['dtype']))
        #np.set_printoptions(threshold=np.nan)
        pprint.pprint(body['data'])
        if is_many:
            print('\n')
