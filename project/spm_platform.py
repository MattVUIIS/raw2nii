import logging
import numpy as np


def spm_platform(*args, **kwargs):
    """
     Platform specific configuration parameters for SPM

     FORMAT ans = spm_platform(arg)
     arg  - optional string argument, can be
            - 'bigend'  - return whether this architecture is bigendian
                          - Inf - is not IEEE floating point
                          - 0   - is little end
                          - 1   - big end
            - 'filesys' - type of filesystem
                          - 'unx' - UNIX
                          - 'win' - DOS
                          - 'mac' - Macintosh
                          - 'vms' - VMS
            - 'sepchar' - returns directory separator
            - 'rootlen' - returns number of chars in root directory name
            - 'user'    - returns username
            - 'tempdir' - returns name of temp directory

     FORMAT PlatFontNames = spm_platform('fonts')
     Returns structure with fields named after the generic (UNIX) fonts, the
     field containing the name of the platform specific font.

     FORMAT PlatFontName = spm_platform('font',GenFontName)
     Maps generic (UNIX) FontNames to platform specific FontNames

     FORMAT SPM_PLATFORM = spm_platform('init',comp)
     Initialises platform specific parameters in global SPM_PLATFORM
     (External gateway to init_platform(comp) subfunction)
     comp         - computer to use [defaults to MatLab's `computer`]
     SPM_PLATFORM - copy of global SPM_PLATFORM

     FORMAT spm_platform
     Initialises platform specific parameters in global SPM_PLATFORM
     (External gateway to init_platform(computer) subfunction)

     FORMAT spm_platform('clear')
     Clears global SPM_PLATFORM containing platform specific parameters

                               ----------------
     SUBFUNCTIONS:

     FORMAT init_platform(comp)
     Initialise platform specific parameters in global SPM_PLATFORM
     comp         - computer to use [defaults to MatLab's `computer`]

    -----------------------------------------------------------------------

     Since calls to spm_platform will be made frequently, most platform
     specific parameters are stored as a structure in the global variable
     SPM_PLATFORM. Subsequent calls use the information from this global
     variable, if it exists.

     Platform specific difinitions are contained in the data structures at
     the beginning of the init_platform subfunction at the end of this
     file.
    _______________________________________________________________________
    """
    logger = logging.getLogger('raw2nii')
    nargin = len(args)
    varargout = []
    global SPM_PLATFORM
    if isempty(SPM_PLATFORM):
        init_platform
    if nargin == 0:
        return varargout
    if 'init' == lower(args[0]):  # (re)initialise
        init_platform(args[1:])
        varargout=[SPM_PLATFORM]
    elif 'clear' == lower(args[0]):  # Clear SPM_PLATFORM
        clear('global','SPM_PLATFORM')
    elif 'bigend' == lower(args[0]):  # Return endian for this architecture
        varargout=[SPM_PLATFORM.bigend]
        if not finite(SPM_PLATFORM.bigend):
            if isnan(SPM_PLATFORM.bigend):
                error(cat(char("I don't know if ""),computer,char('" is big-endian.')))
            else:
                error(cat(char("I don't think that ""),computer,char('" uses IEEE floating point ops.')))
    elif 'filesys' == lower(args[0]):  # Return file system
        varargout=[SPM_PLATFORM.filesys]
    elif 'sepchar' == lower(args[0]):  # Return file separator character
        warning(char('use filesep instead (supported by MathWorks)'))
        varargout=[SPM_PLATFORM.sepchar]
    elif 'rootlen' == lower(args[0]):  # Return length in chars of root directory name
        varargout=[SPM_PLATFORM.rootlen]
    elif 'user' == lower(args[0]):  # Return user string
        varargout=[SPM_PLATFORM.user]
    elif 'tempdir' == lower(args[0]):  # Return temporary directory
        twd=getenv('SPMTMP')
        if isempty(twd):
            twd=copy(tempdir)
        varargout=[twd]
    elif ['font','fonts'] == lower(args[0]):  # Map default font names to platform font names
        if nargin < 2:
            varargout=[SPM_PLATFORM.font]
            return varargout
        if 'times' == lower(args[0]):
            varargout=[SPM_PLATFORM.font.times]
        elif 'courier' == lower(args[0]):
            varargout=[SPM_PLATFORM.font.courier]
        elif 'helvetica' == lower(args[0]):
            varargout=[SPM_PLATFORM.font.helvetica]
        elif 'symbol' == lower(args[0]):
            varargout=[SPM_PLATFORM.font.symbol]
        else:
            logger.warning("Unknown font {0}, using default".format(args[0]))
            varargout = [SPM_PLATFORM.font.helvetica]
    else:  # Unknown Action string
        error(char('Unknown Action string'))
    return varargout
