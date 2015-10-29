from __future__ import division
import logging
import numpy as np


def _set_larger(A, B):
    C = A
    D = B
    if A > B:
        D = C
    else:
        C = D
    return C, D

def calc_angulation(par, angulation):
    logger = logging.getLogger('raw2nii')
    if angulation:
        # trying to incorporate AP FH RL rotation angles: determined using some
        # common sense, Chris Rordon's help + source code and trial and error,
        # this is considered EXPERIMENTAL!
        rads = np.deg2rad(par.angRL)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r1 = np.array([[1, 0, 0], [0, cosrads, -sinrads], [0, sinrads, cosrads]])
        rads = np.deg2rad(par.angAP)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r2 = np.array([[cosrads, 0, sinrads], [0, 1, 0], [-sinrads, 0, cosrads]])
        rads = np.deg2rad(par.angFH)
        cosrads, sinrads = np.cos(rads), np.sin(rads)
        r3 = np.array([[cosrads, -sinrads, 0], [sinrads, cosrads, 0], [0, 0, 1]])
        col = np.array([0, 0, 0, 1])[np.newaxis].T
        R_tot = np.concatenate((np.concatenate((r1.dot(r2).dot(r3),
            [np.zeros(3)])), col), axis=1)
    else:
        R_tot = np.eye(4)
    if 1 == par.sliceorient:  # Traversal
        lmm = np.eye(4)  # Do not rotate
        lXmm = par.fov_apfhrl[2] / par.dim[0]
        lYmm = par.fov_apfhrl[0] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[1] / par.dim[2]
    elif 2 == par.sliceorient:  # Sagittal
        lmm = np.array([[0, 0, -1, 0], [1, 0, 0, 0], [0, -1, 0, 0],
            [0, 0, 0, 1]])
        #Vanderbilt override
        lXmm = par.fov_apfhrl[0] / par.dim[0]
        lYmm = par.fov_apfhrl[1] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[2] / par.dim[2]
    elif 3 == par.sliceorient:  # Coronal
        lmm = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0],
            [0, 0, 0, 1]])  # Rotate 90 degrees
        #Vanderbilt override
        lXmm = par.fov_apfhrl[1] / par.dim[0]
        lYmm = par.fov_apfhrl[2] / par.dim[1]
        #Use smallest in plane resolution...
        lXmm, lYmm = _set_larger(lXmm, lYmm)
        lZmm = par.fov_apfhrl[0] / par.dim[2]
    Zm = np.array([[lXmm, 0, 0, 0], [0, lYmm, 0, 0], [0, 0, lZmm, 0],
        [0, 0, 0, 1]])
    #realvoxsize is used to fill in pixdim nifti header info
    realvoxsize = np.array([lXmm, lYmm, lZmm])
    patient_to_tal = np.diag([-1, -1, 1, 1])
    analyze_to_dicom = np.diag([1, -1, 1, 1])
    A_tot = patient_to_tal.dot(R_tot).dot(Zm).dot(lmm).dot(analyze_to_dicom)
    p_orig = np.array([(par.dim[0] - 1) / 2, (par.dim[1] - 2) / 2,
        (par.dim[2] - 1) / 2, 1])
    offsetA = A_tot.dot(p_orig.T)
    if angulation:
        # trying to incorporate AP FH RL translation: determined using some
        # common sense, Chris Rordon's help + source code and trial and error,
        # this is considered EXPERIMENTAL!
        A_tot[0:3,3] = (np.array([-offsetA[0], -offsetA[1], -offsetA[2]])
            - np.array([par.offRL, par.offAP, -par.offFH]))
    else:
        A_tot[0:3,3] = np.array([-offsetA[0], -offsetA[1], -offsetA[2]])
    HdrMat = A_tot
    return HdrMat, realvoxsize
