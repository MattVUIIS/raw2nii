import logging
import numpy as np


_TYPES = np.array((2, 4, 8, 16, 64, 130, 132, 136, 512, 1024, 2048, 4096,
    16384, 33280, 33792, 34816))
_PREC = ('uint8', 'int16', 'int32', 'float', 'double', 'int8', 'uint16',
    'uint32', 'uint8', 'int16', 'int32', 'float', 'double', 'int8',
    'uint16', 'uint32')
_SWAPPED = np.array((0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1))
_MAXVAL = np.array((2**8 - 1, 2**15 - 1, 2**31 - 1, np.Inf, np.Inf, 2**7 - 1,
    2**16 - 1, 2**32 - 1, 2**8 - 1, 2**15 - 1, 2**31 - 1, np.Inf, np.Inf,
    2**8 - 1, 2**16 - 1, 2**32 - 1))
_MINVAL = np.array((0, -2**15, -2**31, -np.Inf, -np.Inf, -2**7, 0, 0, 0, -2**15,
    -2**31, -np.Inf, -np.Inf, -2**7, 0, 0))
_NANREP = np.array((0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0))
_BITS = np.array((8, 16, 32, 32, 64, 8, 16, 32, 8, 16, 32, 32, 64, 8, 16,
    32))
_INTT = np.array((1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1))


def spm_type(x=None, arg=None):
    """
     Translates data type specifiers between SPM & Matlab representations.
     Returns the type
     x   : specifier
     arg : optional string argument, can be
    	 'swapped' - if type is byteswapped return 1.
    	 'maxval'  - return maximum allowed value.
    	 'minval'  - return minimum allowed value.
    	 'nanrep'  - return 1 if there is a NaN representation.
    	 'bits'    - return the number of bits per voxel.
    	 'intt'    - return 1 if values rounded to nearest integer.
    _______________________________________________________________________

     Original format specifiers are based on ANALYZE.  If the input is
     a number then the corresponding matlab string is returned by default.
     If the input is a string then the appropriate TYPE is returned.
     However, if the optional arg argument is supplied then other
     information will be returned instead.

     With no arguments, a list of data types is returned.

     Additional support was added for signed bytes, unsigned short and
     unsigned int (by adding 128 to the format specifiers for unsigned bytes
     signed short and signed int).  Byte swapped datatypes have the same
     identifiers as the non-byte-swapped versions, multiplied by a factor of
     256.
    """
    if x is None and arg is None:
        return _TYPES

    if isinstance(x, str):
        x = x.rstrip()
        sel = None
        msk = np.where(_SWAPPED == 0)[0]
        for i in msk:
            if _PREC[i] == x:
                sel = i
                break
    else:
        sel = np.where(_TYPES == x)[0]

    if arg is None:
        if isinstance(x, str):
            if sel is None:
                T = np.nan
            else:
                T = _TYPES[sel]
        elif sel is None:
            T = 'unknown'
        else:
            T = _PREC[sel]
    else:
        arg = arg.lower()
        if sel is None:
            T = np.nan
        elif 'swapped' == arg:
            T = _SWAPPED[sel]
        elif 'maxval' == arg:
            T = _MAXVAL[sel]
        elif 'minval' == arg:
            T = _MINVAL[sel]
        elif 'nanrep' == arg:
            T = _NANREP[sel]
        elif 'bits' == arg:
            T = _BITS[sel]
        elif 'intt' == arg:
            T = _INTT[sel]
        else:
            T = np.nan
    return T
