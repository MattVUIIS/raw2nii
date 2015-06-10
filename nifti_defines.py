from __future__ import division


kDT_BINARY = 1
kDT_UNSIGNED_CHAR = 2
kDT_SIGNED_SHORT = 4
kDT_SIGNED_INT = 8
kDT_FLOAT = 16
kDT_COMPLEX = 32
kDT_DOUBLE = 64
kDT_RGB = 128
kDT_INT8 = 256
kDT_UINT16 = 512
kDT_UINT32 = 768
kDT_INT64 = 1024
kDT_UINT64 = 1280
kDT_FLOAT128 = 1536
kDT_COMPLEX128 = 1792
kDT_COMPLEX256 = 2048
kNIFTI_SLICE_SEQ_UNKNOWN = 0
kNIFTI_SLICE_SEQ_INC = 1
kNIFTI_SLICE_SEQ_DEC = 2
kNIFTI_SLICE_ALT_INC = 3
kNIFTI_SLICE_ALT_DEC = 4
kNIFTI_UNITS_UNKNOWN = 0
kNIFTI_UNITS_METER = 1
kNIFTI_UNITS_MM = 2
kNIFTI_UNITS_MICRON = 3
kNIFTI_UNITS_SEC = 8
kNIFTI_UNITS_MSEC = 16
kNIFTI_UNITS_USEC = 24
kNIFTI_UNITS_HZ = 32
kNIFTI_UNITS_PPM = 40
kNIFTI_XFORM_UNKNOWN = 0
kNIFTI_XFORM_SCANNER_ANAT = 1
kNIFTI_XFORM_ALIGNED_ANAT = 2
kNIFTI_XFORM_TALAIRACH = 3
kNIFTI_XFORM_MNI_152 = 4
kNIFTI_MAGIC_SEPARATE_HDR = int('0031696E', 16)
kNIFTI_MAGIC_EMBEDDED_HDR = int('00312B6E', 16)
kswapNIFTI_MAGIC_SEPARATE_HDR = int('6E693100', 16)
kswapNIFTI_MAGIC_EMBEDDED_HDR = int('6E2B3100', 16)
kNIFTI_INTENT_NONE = 0
kNIFTI_INTENT_CORREL = 2
kNIFTI_INTENT_TTEST = 3
kNIFTI_INTENT_FTEST = 4
kNIFTI_INTENT_ZSCORE = 5
kNIFTI_INTENT_CHISQ = 6
kNIFTI_INTENT_BETA = 7
kNIFTI_INTENT_BINOM = 8
kNIFTI_INTENT_GAMMA = 9
kNIFTI_INTENT_POISSON = 10
kNIFTI_INTENT_NORMAL = 11
kNIFTI_INTENT_FTEST_NONC = 12
kNIFTI_INTENT_CHISQ_NONC = 13
kNIFTI_INTENT_LOGISTIC = 14
kNIFTI_INTENT_LAPLACE = 15
kNIFTI_INTENT_UNIFORM = 16
kNIFTI_INTENT_TTEST_NONC = 17
kNIFTI_INTENT_WEIBULL = 18
kNIFTI_INTENT_CHI = 19
kNIFTI_INTENT_INVGAUSS = 20
kNIFTI_INTENT_EXTVAL = 21
kNIFTI_INTENT_PVAL = 22
NIFTI_INTENT_LOGPVAL = 23
NIFTI_INTENT_LOG10PVAL = 24
kNIFTI_LAST_STATCODE = 24
kNIFTI_INTENT_ESTIMATE = 1001
kNIFTI_FIRST_NONSTATCODE = kNIFTI_INTENT_ESTIMATE
kNIFTI_INTENT_LABEL = 1002
kNIFTI_INTENT_NEURONAME = 1003
kNIFTI_INTENT_GENMATRIX = 1004
kNIFTI_INTENT_SYMMATRIX = 1005
kNIFTI_INTENT_DISPVECT = 1006
kNIFTI_INTENT_VECTOR = 1007
kNIFTI_INTENT_POINTSET = 1008
kNIFTI_INTENT_TRIANGLE = 1009
kNIFTI_INTENT_QUATERNION = 1010