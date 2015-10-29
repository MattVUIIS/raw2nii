import binascii
import logging
import numpy as np

import nifti_defines


#Maps bit width -> unsigned datatype
_UDATATYPE_TABLE = {
    8: nifti_defines.kDT_UNSIGNED_CHAR,
    16: nifti_defines.kDT_UINT16,
    32: nifti_defines.kDT_UINT32,
    64: nifti_defines.kDT_UINT64,
}

def WriteNii(fname, NHdr, Data3D, issue):
    logger = logging.getLogger('raw2nii')
    try:
        logger.debug("Writing NHdr...")
        with open(fname, 'wb') as fid:
            NHdr.write_to_file(fid)  # write header to nii binary
            fid.flush()
            if NHdr.bitpix.val in (8, 16, 32, 64):
                if NHdr.bitpix.val == 32 and issue:
                    bitpixstr = 'float'
                else:
                    #Using a form of ternary statements here
                    u_prefix = {True: 'u', False: ''}[
                        NHdr.datatype.val == _UDATATYPE_TABLE[NHdr.bitpix.val]]
                    bitpixstr = '{0}int{1}'.format(u_prefix, NHdr.bitpix.val)
            #now add 4 extra bytes in space between header and offset for data
            #indicating that single .nii file ("n+1\0") rather than separate
            #img/hdr files were written. see http://nifti.nimh.nih.gov
            niftitag = bytearray(binascii.unhexlify('6e2b3100'))
            fid.write(niftitag)
            #add remaining 0s (probably not required)
            if NHdr.vox_offset.val - NHdr.HdrSz.val - 4 > 0:
                arr_len = NHdr.vox_offset.val - NHdr.HdrSz.val - 4
                fid.write(bytearray([0] * arr_len))
            for s3 in range(Data3D.shape[2]):
                #flip data in order to write out
                Data3D[:,:,s3] = np.fliplr(Data3D[:,:,s3])
            #Transpose matrix before writing to get Fortran order
            Data3D.astype(bitpixstr).T.tofile(fid)
    except IOError as e:
        logger.error("Write failed: {0}".format(e))
