from __future__ import division
import logging
import numpy as np
import struct


def WriteNii(fname, NHdr, Data3D):
    logger = logging.getLogger('raw2nii')
    try:
        with open(fname, 'wb') as fid:
            #write header to nii binary
            fn = NHdr.fieldnames()
            for key in fn:
                hdrfield = getattr(NHdr, key)
                prec = hdrfield.prec
                val = hdrfield.val
                if not isinstance(val, np.ndarray):  # If writing simple value
                    if prec == 's':  # Writing a string
                        fid.write(val)
                    else:
                        fid.write(struct.pack(prec, val))
                else: # If writing array
                    if len(val.shape) == 1:  # If 1D array
                        fid.write(struct.pack(prec * len(val),
                            *list(val)))
                    else:
                        #Transpose matrix before writing to get Fortran order
                        val.astype(prec).T.tofile(fid)
                fid.flush()
            if 8 == NHdr.bitpix.val:
                bitpixstr = 'int8'
            elif 16 == NHdr.bitpix.val:
                bitpixstr = 'int16'
            elif 32 == NHdr.bitpix.val:
                bitpixstr = 'int32'
            elif 64 == NHdr.bitpix.val:
                bitpixstr = 'int64'
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
