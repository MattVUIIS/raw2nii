from __future__ import division
import binascii
import logging
import numpy as np
import struct


def WriteNii(fname, NHdr, Data3D):
    logger = logging.getLogger('raw2nii')
    try:
        logger.debug("Writing NHdr...")
        with open(fname, 'wb') as fid:
            #write header to nii binary
            fn = NHdr.fieldnames()
            for key in fn:
                hdrfield = getattr(NHdr, key)
                prec = hdrfield.prec
                val = hdrfield.val
                if not isinstance(val, np.ndarray):  # If writing simple value
                    if prec == 's':  # Writing a string
                        packed = val
                        logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                            "packed string='{4}')".format(key, fid.tell(), prec,
                            val, packed))
                        fid.write(packed)
                    else:
                        packed = struct.pack(prec, val)
                        logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                            "packed bytes='{4}')".format(key, fid.tell(), prec,
                            val, binascii.hexlify(packed)))
                        fid.write(packed)
                else: # If writing array
                    if len(val.shape) == 1:  # If 1D array
                        packed = struct.pack(prec * len(val), *list(val))
                        logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                            "packed array='{4}')".format(key, fid.tell(), prec,
                            val, binascii.hexlify(packed)))
                        fid.write(packed)
                    else:
                        #Transpose matrix before writing to get Fortran order
                        packed = val.astype(prec).T
                        logger.debug("Field {0}[{1}] (prec={2}, val='{3}', "
                            "packed matrix='{4}')".format(key, fid.tell(), prec,
                            val, binascii.hexlify(packed.tobytes())))
                        packed.tofile(fid)
                fid.flush()
            if NHdr.bitpix.val in (8, 16, 32, 64):
                bitpixstr = 'int{0}'.format(NHdr.bitpix.val)
            #now add 4 extra bytes in space between header and offset for data
            #indicating that single .nii file ("n+1\0") rather than separate
            #img/hdr files were written. see http://nifti.nimh.nih.gov
            niftitag = bytearray([int('6e', 16), int('2b', 16), int('31', 16),
                int('00', 16)])
            fid.write(niftitag)
            #add remaining 0s (probably not required)
            if NHdr.vox_offset.val - NHdr.HdrSz.val - 4 > 0:
                arr_len = NHdr.vox_offset.val - NHdr.HdrSz.val - 4
                fid.write(bytearray([0] * arr_len))
            for s3 in range(0, Data3D.shape[2]):
                #flip data in order to write out
                Data3D[:,:,s3] = np.fliplr(Data3D[:,:,s3])
            #Transpose matrix before writing to get Fortran order
            Data3D.astype(bitpixstr).T.tofile(fid)
    except IOError as e:
        logger.error("Write failed: {0}".format(e))
