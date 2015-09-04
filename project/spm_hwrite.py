import logging
import numpy as np
import os
import re
import struct
import sys

from spm_type import spm_type


def spm_hwrite(P, DIM, VOX, SCALE, TYPE, OFFSET, ORIGIN=None,
        DESCRIP='spm compatible'):
    """ writes a header
     (function copied from spm99, so spm99 does not have to be present)
     FORMAT [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);

     P       - filename 	     (e.g 'spm' or 'spm.img')
     DIM     - image size       [i j k [l]] (voxels)
     VOX     - voxel size       [x y z [t]] (mm [sec])
     SCALE   - scale factor
     TYPE    - datatype (integer - see spm_type)
     OFFSET  - offset (bytes)
     ORIGIN  - [i j k] of origin  (default = [0 0 0])
     DESCRIP - description string (default = 'spm compatible')

     s       - number of elements successfully written (should be 348)
    __________________________________________________________________________

     spm_hwrite writes variables from working memory into a SPM/ANALYZE
     compatible header file.  The 'originator' field of the ANALYZE format has
     been changed to ORIGIN in the SPM version of the header. funused1
     of the ANALYZE format is used for SCALE

     see also dbh.h (ANALYZE) spm_hread.m and spm_type.m

    __________________________________________________________________________
     @(#)spm_hwrite.m	2.2 99/10/29


     ensure correct suffix {.hdr} and open header file
    ---------------------------------------------------------------------------
    """
    logger = logging.getLogger('raw2nii')
    if ORIGIN is None:
        origin = np.array([0, 0, 0, 0, 0])
    else:
        origin = np.array([ORIGIN[:].T, 0, 0])
    #Remove spaces, file extension, and add '.hdr'
    P, _ = os.path.splitext(re.sub(' ', '', P))
    P += '.hdr'
    #For byte swapped data-types, also swap the bytes around in the headers.
    mach = 'native'
    if spm_type(TYPE, 'swapped'):
        if sys.byteorder == 'big':
            mach = 'ieee-le'
        else:
            mach = 'ieee-be'
        TYPE = spm_type(spm_type(TYPE))
    try:
        with open(P, 'wb') as fid:
            data_type = 'dsr      ' + chr(0)
            P += '                  '
            db_name = P[0:17] + chr(0)
            #set header variables
            DIM = DIM[:].T
            if DIM.shape[1] < 4:
                DIM = np.array([DIM, 1])
            VOX = VOX[:].T
            if VOX.shape[1] < 4:
                VOX = np.array([VOX, 0])
            dim = np.array([4, DIM[0:4], 0, 0, 0])
            pixdim = np.array([0, VOX[0:4], 0, 0, 0])
            vox_offset = OFFSET
            funused1 = SCALE
            glmax = 1
            glmin = 0
            bitpix = 0
            descrip = np.zeros((1, 80))
            aux_file = 'none                   ' + chr(0)
            origin = np.array([0, 0, 0, 0, 0])
            if TYPE == 1:
                bitpix = 1
                glmax = 1
                glmin = 0
            elif TYPE == 2:
                bitpix = 8
                glmax = 255
                glmin = 0
            elif TYPE == 4:
                bitpix = 16
                glmax = 32767
                glmin = 0
            elif TYPE == 8:
                bitpix = 32
                glmax = 2 ** 31 - 1
                glmin = 0
            elif TYPE == 16:
                bitpix = 32
                glmax = 1
                glmin = 0
            elif TYPE == 64:
                bitpix = 64
                glmax = 1
                glmin = 0
            d = range(0, min(len(DESCRIP), 79))
            descrip[d] = DESCRIP[d]
            #write (struct) header_key
            fid.seek(0)
            fid.write(struct.pack('i', 348))
            fid.write(data_type)
            fid.write(db_name)
            fid.write(struct.pack('ihcc', 0, 0, 'r', '0'))
            #write (struct) image_dimension
            fid.seek(40)
            fid.write(struct.pack('hccBB', dim, 'm', 'm', 0, 0))
            fid.write(np.zeros((1, 8)))
            fid.write(struct.pack('hhhhfffffffiiii', 0, TYPE, bitpix, 0,
                pixdim, vox_offset, funused1, 0, 0, 0, 0, 0, 0, glmax, glmin))
            #write (struct) image_dimension
            fid.write(descrip)
            fid.write(aux_file)
            fid.write(struct.pack('Bh', 0, origin))
            fid.write(np.zeros((1, 85)))
    except IOError as e:
        logger.error("Error opening/writing {0}: {1}".format(P, e))
    return s
