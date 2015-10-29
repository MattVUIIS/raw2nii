import argparse
import logging
import numpy as np
import pprint

from nii_info import read_nii


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    options = parser.parse_args()
    header, body = read_nii(options.filename)
    pprint.pformat('dtype: {0}'.format(body['dtype']))
    #np.set_printoptions(threshold=np.nan)
    pprint.pprint(body['data'])
