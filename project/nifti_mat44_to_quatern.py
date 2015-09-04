from __future__ import division
import logging
import numpy as np


def nifti_mat44_to_quatern(A):
    logger = logging.getLogger('raw2nii')
    A = A.copy()
    #offset outputs are read write out of input matrix
    qoffset_xyz = A[0:3,3]
    #compute lengths of each column; these determine grid spacings
    d1 = np.sqrt(np.sum(A[0:3,0:3]*A[0:3, 0:3], axis=0))
    #if a column length is zero, patch the trouble
    if d1[0] == 0:
        A[:,0] = [[1], [0], [0]]
        d1[0] = 1
    if d1[1] == 0:
        A[:,1] = [[0], [1], [0]]
        d1[1] = 1
    if d1[2] == 0:
        A[:,2] = [[0], [0], [1]]
        d1[2] = 1
    #normalize the columns
    A[0:3,0] = A[0:3,0] / d1[0]
    A[0:3,1] = A[0:3,1] / d1[1]
    A[0:3,2] = A[0:3,2] / d1[2]
    # At this point, the matrix has normal columns, but we have to allow
    # for the fact that the hideous user may not have given us a matrix
    # with orthogonal columns.
    # So, now find the orthogonal matrix closest to the current matrix.
    # One reason for using the polar decomposition to get this
    # orthogonal matrix, rather than just directly orthogonalizing
    # the columns, is so that inputting the inverse matrix to R
    # will result in the inverse orthogonal matrix at this point.
    # If we just orthogonalized the columns, this wouldn't necessarily hold.
    Q = A[0:3,0:3]
    U, S ,V = np.linalg.svd(Q)
    P = U.dot(V)
    #                            [ r11 r12 r13 ]
    # at this point, the matrix  [ r21 r22 r23 ] is orthogonal
    #                            [ r31 r32 r33 ]
    #
    # compute the determinant to determine if it is proper
    zd = np.linalg.det(P)
    if zd > 0:
        qfac = 1.0
    else: # improper ==> flip 3rd column
        qfac = -1.0
        P[:,2] = -P[:,2]
    #now, compute quaternion parameters
    a = np.trace(P) + 1
    if a > 0.5:  # simplest case
        a = 0.5 * np.sqrt(a)
        b = 0.25 * (P[2,1] - P[1,2]) / a
        c = 0.25 * (P[0,2] - P[2,0]) / a
        d = 0.25 * (P[1,0] - P[0,1]) / a
    else:  # trickier case
        xd = 1.0 + P[0,0] - (P[1,1] + P[2,2])
        yd = 1.0 + P[1,1] - (P[0,0] + P[2,2])
        zd = 1.0 + P[2,2] - (P[0,0] + P[1,1])
        if xd > 1.0:
            b = 0.5 * np.sqrt(xd)
            c = 0.25 * (P[0,1] + P[1,0]) / b
            d = 0.25 * (P[0,2] + P[2,0]) / b
            a = 0.25 * (P[2,1] - P[1,2]) / b
        elif yd > 1.0:
            c = 0.5 * np.sqrt(yd)
            b = 0.25 * (P[0,1] + P[1,0]) / c
            d = 0.25 * (P[1,2] + P[2,1]) / c
            a = 0.25 * (P[0,2] - P[2,0]) / c
        else:
            d = 0.5 * np.sqrt(zd)
            b = 0.25 * (P[0,2] + P[2,0]) / d
            c = 0.25 * (P[1,2] + P[2,1]) / d
            a = 0.25 * (P[1,0] - P[0,1]) / d
        if a < 0.0:
            b = -b
            c = -c
            d = -d
            a = -a
    quatern_bcd = np.array([b, c, d])
    return qoffset_xyz, quatern_bcd, qfac
