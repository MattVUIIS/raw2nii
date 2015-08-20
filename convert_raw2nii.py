from __future__ import division
import argparse
import array
import logging
import numpy as np
import os

from create_nii_hdr import CreateNiiHdr
from read_par import read_par
from spm_type import spm_type
from write_nii import WriteNii


def convert_raw2nii(filelist, prefix="", suffix="", pathpar="",
        outfolder=None, outputformat=None, angulation=None,
        rescale=None, dim=3, dti_revertb0=0):
    """
        filelist: list containing PAR file names (without path info) to convert
            into analyze/nifti
        prefix       : characters to prepend to all output filenames.
                       Blank by default. The prefix, PAR filename (without the
                       extension), suffix, and the file number are used
                       respectively to form the output filename.
        suffix       : characters to append to all output filenames.
                       Blank by default. The prefix, PAR filename (without the
                       extension), suffix, and the file number are used
                       respectively to form the output filename.
        pathpar      : complete path containing PAR files (with trailing /)
        outfolder    : when set, files will be written to this directory,
                       including lowest level folder containing parfile
        outputformat : 1 for Nifty output format (spm5), 2 for Analyze (spm2)
        angulation   : when 1: include affine transformation as defined in PAR
                       file in hdr part of Nifti file (nifti only, EXPERIMENTAL!)
        rescale      : when 1: store intensity scale as found in PAR
                       file (assumed equall for all slices). Yields DV values.
        dim          : when 3, single 3D nii files will be produced, when
                       4, one 4D nii file will be produced, for example
                       for time series or dti data
        dti_revertb0 : when 0 (default), philips ordering is used for DTI data
                       (eg b0 image last). When 1, b0 is saved as first image
                       in 3D or 4D data
    """
    logger = logging.getLogger('raw2nii')
    pathpar = os.path.expanduser(pathpar)
    outfiles = []
    for parfilename in filelist:
        full_parfilename = os.path.join(pathpar, parfilename)
        Parameters = read_par(full_parfilename)
        if Parameters.problemreading == 1:
            logger.warning('Skipping volume {0} because of reading errors.'
                .format(full_parfilename))
        else:
            if outfolder:
                absDir = os.path.abspath(os.path.expanduser(outfolder))
            else:
                absDir = os.path.abspath('NIFTI')
            if not os.path.exists(absDir):
                logger.debug("Makedirs on: {0}".format(absDir))
                os.makedirs(absDir)
            logger.debug("Check dir: {0}".format(absDir))
            Vox = Parameters.vox
            Recfile = parfilename
            basename, ext = os.path.splitext(os.path.basename(Recfile))
            if '.par' == ext:
                Recfile = basename + '.rec'
            elif '.PAR' == ext:
                Recfile = basename + '.REC'
            outBaseName, _ = os.path.splitext(os.path.basename(parfilename))
            outFileName = prefix + outBaseName + suffix
            Precision = 'int' + str(Parameters.bit)
            Size = Parameters.dim[0] * Parameters.dim[1] * Parameters.dim[2]
            chunk_size = Size * (Parameters.bit // 8)
            SizeSlice = Parameters.dim[0] * Parameters.dim[1]
            slice_size = SizeSlice * (Parameters.bit // 8)
            with open(os.path.join(pathpar, Recfile), 'rb') as ID1:
                Dim = Parameters.dim
                cDim = Dim
                VolData = np.zeros((Dim[0], Dim[1], Dim[2]))
                Type = spm_type(Precision)
                Offset = 0
                Descrip = Parameters.name
                Scale = 1
                Orign = Parameters.fov / Vox / 2
                BytesPerValue = Parameters.bit // 8
                if 'V3' == Parameters.ResToolsVersion:
                    nii_hdr_dim = np.array([dim, 1, 1, 1, 1, 1, 1])
                    NHdr = CreateNiiHdr(Parameters, angulation, rescale,
                        nii_hdr_dim)
                    logger.info(' Start to convert scan: {0}'.format(Recfile))
                    for j in range(1, Parameters.dyn):
                        if 1 == Parameters.slicessorted:
                            Data = ID1.read(chunk_size)
                            Inputvolume = np.zeros((Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                            Inputvolume = np.reshape(Data, (Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                        elif 2 == Parameters.slicessorted:
                            InputVolume = None #clear('InputVolume')
                            Inputvolume = np.zeros((Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                            for slice_ in range(1, Parameters.dim(2)):
                                ID1.seek((j - 1 + Parameters.dyn * (slice_ - 1))
                                    * SizeSlice * BytesPerValue)
                                Inputvolume[:,:,slice_] = ID1.read(slice_size)
                            Inputvolume = np.reshape(Inputvolume,
                                (Parameters.dim(1), Parameters.dim(0),
                                Parameters.dim(2)))
                        VolNameSinExt = os.path.join(absDir,
                            outFileName + '-{0:03.0f}'.format(j))
                        if 1 == outputformat:
                            VolName = VolNameSinExt + '.nii'
                            logger.info('Writing file: {0}...'.format(VolName))
                            WriteNii(VolName, NHdr, Inputvolume)
                            logger.info('  ...done')
                            outfiles.append(VolName)
                        elif 2 == outputformat:
                            VolName = VolNameSinExt + '.img'
                            logger.info('Writing file: {0}...'.format(VolName))
                            with open(VolName, 'wb') as ID2:
                                ID2.write(Inputvolume[:,:,:])
                            P = VolName
                            spm_hwrite(P, Dim, Vox, Scale, Type, Offset,
                                round(Orign),Descrip)
                            logger.info('  ...done')
                        VolData = np.array([])
                        outfiles.append(VolName)
                        outfcount += 1
                        logger.info('Write file: {0}-{1}'.format(outFileName,
                            '{0:03.0f}'.format(j)))
                elif Parameters.ResToolsVersion in ('V4', 'V4.1', 'V4.2'):
                    if angulation == 1:
                        logger.warning("Assuming angulation parameters are "
                            "identical for all scans in (4D) volume!")
                    if rescale == 1:
                        logger.warning("Assuming rescaling parameters "
                            "(see PAR-file) are identical for all slices in "
                            "volume and all scans in (4D) volume!")
                    iSlice = Parameters.slice_index
                    #add column containing size of slice in bytes
                    iSlice = np.concatenate((iSlice, (iSlice[:,7] *
                        iSlice[:,9] * iSlice[:,10] // 8)[np.newaxis].T), axis=1)
                    #get order in which slices are stored in REC file.
                    order_slices = iSlice[:,6]
                    i = np.argsort(order_slices)
                    order_slices = np.sort(order_slices)  # sort them
                    #sort bytes per slice info acc. to index in REC file
                    bytespslice_sorted = iSlice[i,12]
                    fileposSlice_sorted = np.cumsum(bytespslice_sorted)
                    fileposSlice_sorted = np.concatenate(([0],
                        fileposSlice_sorted[:-1]))
                    index_orig = np.array(range(0, len(order_slices)))
                    fileposSlice = fileposSlice_sorted[index_orig[i]]
                    iSlice = np.concatenate((iSlice,
                        fileposSlice[np.newaxis].T), axis=1)
                    #Now sort entire slice_index according to dynamics, diff.
                    #gradient(fastest varying) and mr_type parameters (slowest
                    #varying). To do this we must convert to recarray for
                    #np.argsort function
                    rectypes = [('slice_number', int), ('echo_number', int),
                        ('dynamic_scan_number', int),
                        ('cardiac_phase_number', int),
                        ('image_type_mr', int), ('scanning_sequence', int),
                        ('index_in_REC_file', int), ('image_pixel_size', int),
                        ('scan_percentage', int), ('recon_resolution_x', int),
                        ('recon_resolution_y', int),
                        ('diffusion_b_value_number', int)]
                    #par.slice_index_as_rec = np.rec.fromrecords(rows, dtype=rectypes)
                    iSlices_sorted = iSlice[
                        np.lexsort((iSlice[:,0], iSlice[:,2],
                        iSlice[:,11], iSlice[:,1], iSlice[:,4], iSlice[:,5]))]
                    nLine = 0
                    nr_mrtypes = len(np.unique(iSlices_sorted[:,5]))
                    nr_echos = len(np.unique(iSlices_sorted[:,1]))
                    nr_realmrtypes = len(np.unique(iSlices_sorted[:,4]))
                    nr_diffgrads = len(np.unique(iSlices_sorted[:,11]))
                    nr_dyn = len(np.unique(iSlices_sorted[:,2]))
                    if 3 == dim:
                        #Compare dynamic, mr_type, realmr_type, echo, diff grad
                        slice_bool = (
                            np.not_equal(0, np.diff(iSlices_sorted[:,2])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,5])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,1])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,4])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,11]))
                            )[np.newaxis].T
                        NewFile = np.concatenate((slice_bool, [[1]]))
                        VolData = np.zeros((Dim[0], Dim[1], Dim[2]))
                    elif 4 == dim:
                        slice_bool = (
                            np.not_equal(0, np.diff(iSlices_sorted[:,5])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,1])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,4])) |
                            np.not_equal(0, np.diff(iSlices_sorted[:,11]))
                            )[np.newaxis].T
                        NewFile = np.concatenate((slice_bool, [[1]]))
                        expnrslices = Dim[2] * nr_dyn
                        bytes_required = expnrslices * Dim[0] * Dim[1] * 8
                        logger.warning("4D data conversion to nii may require "
                            "huge amounts of memory (RAM)! Your mileage: {0} "
                            "bytes required".format(bytes_required))
                        VolData = np.zeros((Dim[0], Dim[1], Dim[2] *
                            Parameters.dyn))
                    elif 5 == dim:
                        NewFile = np.zeros((iSlices_sorted.shape[0], 1))
                        NewFile[-1] = 1
                        expnrslices = (Dim[2] * nr_mrtypes * nr_echos *
                            nr_realmrtypes * nr_diffgrads * nr_dyn)
                        bytes_required = expnrslices * Dim[0] * Dim[1] * 8
                        logger.warning("ND data conversion to nii may require "
                            "huge amounts of memory (RAM)! Your mileage: {0} "
                            "bytes required".format(bytes_required))
                        VolData = np.zeros((Dim[0], Dim[1], expnrslices))
                    slicenr = 0
                    while nLine < iSlices_sorted.shape[0]:
                        ID1.seek(iSlices_sorted[nLine,13])
                        cDim = np.concatenate((iSlices_sorted[nLine,9:11],
                            [Dim[2]]))
                        type_map = {8: 'b', 16: 'h', 32: 'i'}
                        SliceData = array.array(type_map[Parameters.bit])
                        SliceData.fromfile(ID1, cDim[0] * cDim[1])
                        ImageSlice = np.reshape(SliceData, (cDim[0], cDim[1]),
                            order='F')
                        VolData[:,:,slicenr] = ImageSlice
                        slicenr += 1
                        if NewFile[nLine]:
                            if nr_dyn > 1:
                                dyn_suffix = "-{0:04d}".format(
                                    iSlices_sorted[nLine,2])
                                dyn_ndsuffix = "-d{0:04d}".format(nr_dyn)
                            else:
                                dyn_suffix = "-{0:04d}".format(1)
                                dyn_ndsuffix = ""
                            if nr_mrtypes > 1:
                                mrtype_suffix = "-s{0:3d}".format(
                                    iSlices_sorted[nLine,5])
                                mrtype_ndsuffix = "-s{0:03d}".format(nr_mrtypes)
                            else:
                                mrtype_suffix = ""
                                mrtype_ndsuffix = ""
                            if nr_realmrtypes > 1:
                                realmrtype_suffix = "-t{0:03d}".format(
                                    iSlices_sorted[nLine,4])
                                realmrtype_ndsuffix = "-t{0:03d}".format(
                                    nr_realmrtypes)
                            else:
                                realmrtype_suffix = ""
                                realmrtype_ndsuffix = ""
                            if nr_echos > 1:
                                echo_suffix = "-e{0:03d}".format(
                                    iSlices_sorted[nLine,1])
                                echo_ndsuffix = "-e{0:03d}".format(nr_echos)
                            else:
                                echo_suffix = ""
                                echo_ndsuffix = ""
                            if nr_diffgrads > 1:
                                diffgrad_suffix = "-g{0:03d}".format(
                                    iSlices_sorted[nLine,11])
                                diffgrad_ndsuffix = "-g{0:03d}".format(
                                    nr_diffgrads)
                            else:
                                diffgrad_suffix = ""
                                diffgrad_ndsuffix = ""
                            cDim = np.concatenate((iSlices_sorted[nLine,9:11],
                                [Dim[2]]))
                            nii_hdr_dim = np.array([dim, nr_dyn, nr_diffgrads,
                                nr_echos, nr_mrtypes, nr_realmrtypes])
                            NHdr = CreateNiiHdr(Parameters, angulation, rescale,
                                nii_hdr_dim)
                            if 3 == dim:
                                VolNameSinExt = os.path.join(absDir, outFileName +
                                    echo_suffix + realmrtype_suffix +
                                    mrtype_suffix + diffgrad_suffix +
                                    dyn_suffix)
                            elif 4 == dim:
                                VolNameSinExt = os.path.join(absDir, outFileName +
                                    echo_suffix + realmrtype_suffix +
                                    mrtype_suffix + diffgrad_suffix +
                                    dyn_ndsuffix)
                            elif 5 == dim:
                                if dti_revertb0:
                                    logger.warning("Doing a very dirty hack "
                                        "to DTI data, putting b0 data as first "
                                        "in ND nii file. Only use this "
                                        "(hidden) raw2nii option when working "
                                        "at the UMC Utrecht with Rene Mandl's "
                                        "DTI tools")
                                    nslice = Parameters.dim[2]
                                    B0 = VolData[:,:,-1 - nslice + 1:-1]
                                    VolData[:,:,nslice + 1: end()] = (
                                        VolData[:,:,0:-1 - nslice])
                                    VolData[:,:,0:nslice] = B0
                                    NHdr.dim.val[0] = 4
                                    d = NHdr.dim.val(4)
                                    NHdr.dim.val[4] = NHdr.dim.val(5)
                                    NHdr.dim.val[5] = d
                                    NHdr.pixdim.val[4] = 1
                                VolNameSinExt = os.path.join(absDir, outFileName +
                                    echo_ndsuffix + realmrtype_ndsuffix +
                                    mrtype_ndsuffix + diffgrad_ndsuffix +
                                    dyn_ndsuffix)
                            logger.info("VolNameSinExt: {0}".format(VolNameSinExt))
                            Slice = np.zeros((cDim[0], cDim[1]))
                            if 1 == outputformat:
                                VolName = VolNameSinExt + '.nii'
                                logger.info("Writing file: {0}...".format(
                                    VolName))
                                WriteNii(VolName, NHdr, VolData)
                                logger.info("  ...done")
                                outfiles.append(VolName)
                            elif 2 == outputformat:
                                VolName = VolNameSinExt + '.img'
                                logger.info("Writing file: {0}...".format(
                                    VolName))
                                with open(VolName, 'wb') as ID2:
                                    ID2.write(VolData[:,:,:])
                                P = VolName
                                spm_hwrite(P, Dim, Vox,
                                    Parameters.rescale_slope, Type,
                                    Parameters.rescale_interc, round(Orign),
                                    Descrip)
                                logger.info("  ...done")
                                outfiles.append(VolName)
                            slicenr = 0
                        nLine += 1
                else:
                    logger.warning("Sorry, but data format extracted using "
                        "Philips Research File format {0} was not known at "
                        "the time the raw2nii software was developed".format(
                        Parameters.ResToolsVersion))
    return outfiles

if __name__ == "__main__":
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.DEBUG)
    _formatter = logging.Formatter("%(levelname)s %(asctime)s %(filename)s: "
        "%(message)s")
    _stream_handler = logging.StreamHandler()
    _stream_handler.setFormatter(_formatter)
    logger.addHandler(_stream_handler)
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", nargs="+", type=str)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--suffix", default="")
    parser.add_argument("--pathpar", default="")
    parser.add_argument("--outfolder", default=None)
    parser.add_argument("--outputformat", type=int, default=1)
    parser.add_argument("--angulation", type=int, default=1)
    parser.add_argument("--rescale", type=int, default=1)
    parser.add_argument("--dim", type=int, default=3)
    parser.add_argument("--dti_revertb0", type=int, default=0)
    options = parser.parse_args()
    convert_raw2nii(**vars(options))
