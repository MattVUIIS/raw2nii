from __future__ import division
import argparse
import array
import logging
import numpy as np
import os
import re

from create_nii_hdr import CreateNiiHdr
from read_par import read_par
from spm_type import spm_type
from write_nii import WriteNii


def convert_raw2nii(filelist, prefix, suffix, pathpar, outfolder, outputformat,
        angulation, rescale, dti_revertb0):
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
        dti_revertb0 : when 0 (default), philips ordering is used for DTI data
                       (eg b0 image last). When 1, b0 is saved as first image
                       in 3D or 4D data
    """
    logger = logging.getLogger('raw2nii')
    dim = 4  # always use 4D
    pathpar = os.path.expanduser(pathpar)
    outfiles = []
    outfcount = 1
    files_error_size = []
    count_error = 1
    for nbf, parfilename in enumerate(filelist):
        full_parfilename = os.path.join(pathpar, parfilename)
        #extract the bval and bvec from the PAR file
        Parameters = read_par(full_parfilename)
        if np.allclose(Parameters.gen_info.diffusion, 1): #it's a dti file
            dti_revertb0 = 1
            NumberOfVolumes = Parameters.NumberOfVolumes
            Img_size = Parameters.slices.shape[0]
            numberofslices_from_header = (
                Parameters.gen_info.max_number_of_slices_locations)
            numberofslices = Img_size // NumberOfVolumes
            bval = np.zeros((NumberOfVolumes, 1))
            bvec = np.zeros((NumberOfVolumes, 3))
            if not np.allclose(numberofslices_from_header, numberofslices):
                logger.warning('DTI incomplete. Number of slices different '
                    'from header and reality.')
                sl_i = np.arange(0, NumberOfVolumes *
                    numberofslices_from_header, numberofslices_from_header)
            else:  # No error with the number of slices
                sl_i = np.arange(0, NumberOfVolumes * numberofslices,
                    numberofslices)
            for volumeIndex in range(NumberOfVolumes):
                sl = Parameters.slices[sl_i[volumeIndex]]
                bval[volumeIndex, 0] = sl.diffusion_b_factor
                # ap fh rl -> rl ap fh
                bvec[volumeIndex, 0] = sl.diffusion_rl
                bvec[volumeIndex, 1] = sl.diffusion_ap
                bvec[volumeIndex, 2] = sl.diffusion_fh
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
            Recfile, ext = os.path.splitext(os.path.join(pathpar, parfilename))
            if '.par' == ext:
                Recfile += '.rec'
            elif '.PAR' == ext:
                Recfile += '.REC'
            outBaseName, _ = os.path.splitext(os.path.basename(parfilename))
            outFileName = prefix + outBaseName + suffix
            Precision = 'int' + str(Parameters.bit)
            Size = Parameters.dim[0] * Parameters.dim[1] * Parameters.dim[2]
            chunk_size = Size * (Parameters.bit // 8)
            SizeSlice = Parameters.dim[0] * Parameters.dim[1]
            slice_size = SizeSlice * (Parameters.bit // 8)
            with open(Recfile, 'rb') as ID1:
                Dim = Parameters.dim
                cDim = Dim
                VolData = np.zeros((Dim[0], Dim[1], Dim[2]))
                Type = spm_type(Precision)
                Offset = 0
                Descrip = Parameters.gen_info.protocol_name
                Scale = 1
                Orign = Parameters.fov / Vox / 2
                BytesPerValue = Parameters.bit // 8
                if 'V3' == Parameters.ResToolsVersion:
                    nii_hdr_dim = np.array([dim, 1, 1, 1, 1, 1])
                    NHdr = CreateNiiHdr(Parameters, angulation, nii_hdr_dim)
                    logger.info('Start to convert scan: {0}'.format(Recfile))
                    for j in range(1, Parameters.gen_info.max_number_of_dynamics):
                        if 1 == Parameters.slices_sorted:  # Ascending order
                            Data = ID1.read(chunk_size)
                            Inputvolume = np.zeros((Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                            Inputvolume = np.reshape(Data, (Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                        elif 2 == Parameters.slicessorted:  # Descending order
                            InputVolume = None #clear('InputVolume')
                            Inputvolume = np.zeros((Parameters.dim(1),
                                Parameters.dim(0), Parameters.dim(2)))
                            for slice_ in range(1, Parameters.dim(2)):
                                ID1.seek((j - 1 +
                                    Parameters.gen_info.max_number_of_dynamics *
                                    (slice_ - 1)) * SizeSlice * BytesPerValue)
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
                    #new: loop slices (as in slice_index) and open and close
                    #files along the way (according to info in index on dynamic
                    #and mr_type)
                    if angulation:
                        logger.warning("Assuming angulation parameters are "
                            "identical for all scans in (4D) volume!")
                    if rescale:
                        logger.warning("Assuming rescaling parameters "
                            "(see PAR-file) are identical for all slices in "
                            "volume and all scans in (4D) volume!")
                    iSlice = Parameters.slices
                    #add column containing size of slice in bytes
                    #byte_size = np.array(iSlice.image_pixel_size *
                    #    iSlice.recon_resolution_x * iSlice.recon_resolution_y //
                    #    8)
                    #bytespslice = np.concatenate((
                    #    iSlice.slice_number[np.newaxis].T,
                    #    byte_size[np.newaxis].T), axis=1)
                    #bytespslice = np.concatenate((iSlice.slice_number[np.newaxis].T,
                    #    (iSlice.image_pixel_size *
                    #    iSlice.recon_resolution_x * iSlice.recon_resolution_y //
                    #    8)[np.newaxis].T), axis=1)
                    #get order in which slices are stored in REC file.
                    #order_slices = iSlice.index_in_rec_file
                    #i = np.argsort(order_slices)
                    #order_slices = np.sort(order_slices)  # sort them
                    #sort bytes per slice info acc. to index in REC file
                    #bytespslice_sorted = iSlice[i,13]
                    #fileposSlice_sorted = np.cumsum(bytespslice_sorted)
                    #fileposSlice_sorted = np.concatenate(([0],
                    #    fileposSlice_sorted[:-1]))
                    #index_orig = np.arange(order_slices.shape[0])
                    #fileposSlice = fileposSlice_sorted[index_orig[i]]
                    #Add column containing start position in bytes of this slice
                    #in the file
                    #iSlice = np.concatenate((iSlice,
                    #    fileposSlice[np.newaxis].T), axis=1)
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
                    iSlices_sorted = iSlice[np.lexsort((iSlice.slice_number,
                        iSlice.dynamic_scan_number,
                        iSlice.diffusion_b_value_number,
                        iSlice.gradient_orientation_number, iSlice.echo_number,
                        iSlice.image_type_mr, iSlice.scanning_sequence))]
                    nLine = 0
                    #Determine number of interleaved image sequences (was:types,
                    #name kept for historic reasons) (e.g. angio)
                    nr_mrtypes = np.unique(iSlices_sorted.scanning_sequence
                        ).shape[0]
                    #Determine number of interleaved echos
                    nr_echos = np.unique(iSlices_sorted.echo_number).shape[0]
                    #Determine number of interleaved image types (e.g. angio)
                    nr_realmrtypes = np.unique(iSlices_sorted.image_type_mr
                        ).shape[0]
                    #Determine number of diffusion gradients (e.g. DTI)
                    nr_diffgrads = np.unique(
                        iSlices_sorted.gradient_orientation_number).shape[0]
                    #Determine number of dynamics(directly from slice lines in
                    #PAR file instead of PAR file header info!)
                    nr_dyn = np.unique(iSlices_sorted.dynamic_scan_number
                        ).shape[0]
                    if nr_dyn != Parameters.gen_info.max_number_of_dynamics:
                        files_error_size.append(nbf)
                    if dti_revertb0:
                        if nr_diffgrads != NumberOfVolumes:
                            logger.warning('VANDERBILT hack, the number of '
                                'diffusion gradients is not coherent, taking '
                                'the info from PAR header.')
                            #Problem with iSlices_sorted that point the b0 at
                            #Sorted the list to keep the b0 at the end
                            nslice = Parameters.dim[2]
                            #Get the first index for the b0
                            r = iSlices_sorted.diffusion_b_value_number
                            index = np.extract(r == 1, np.arange(r.shape[0]))
                            #Put the info at the end
                            B0 = np.copy(sl)
                            iSlices_sorted[index] = iSlices_sorted[-nslice:]
                            iSlices_sorted[-nslice:] = B0
                            #Keep only the number of volumes from par header
                            nr_diffgrads = NumberOfVolumes
                        #Starting at zero, so +1
                        b0moved = (iSlices_sorted[-1].index_in_rec_file + 1 ==
                            Img_size)
                        #%%%%%%%%%%%%%%% WRITING BVAL/BVEC %%%%%%%%%%%%%%%%
                        #Write the bval after sorting the slices
                        #Put the b0 first if needed from the iSliced_Sorted in
                        #the bval and bvec
                        volumeRange = list(range(NumberOfVolumes))
                        if b0moved:  # Write the last value first
                            logger.warning('VANDERBILT hack -> putting last '
                                'value in front for bval/bvec.')
                            volumeRange = ([volumeRange[-1]]
                                + volumeRange[:-1])
                        pathstr, name = os.path.split(parfilename)
                        name, ext = os.path.splitext(name)
                        bval_filename = (os.path.join(absDir, name) +
                            '-x-bval.txt')
                        try:
                            with open(bval_filename, 'wb') as bval_file:
                                bval_file.writelines('{0:.6f}'.format(x) + ' '
                                    for x in bval[volumeRange,0])
                        except OSError:
                            logger.error('File name for bval.txt invalid')
                        #Write the bvec
                        bvec_filename = (os.path.join(absDir, name) +
                            '-x-bvec.txt')
                        try:
                            with open(bvec_filename, 'wb') as bvec_file:
                                bvec_file.writelines('{0:.6f}'.format(x) + ' '
                                    for x in bvec[volumeRange,0])
                                bvec_file.write('\n')
                                #After checking the dtiqa process, need to flip
                                #Y data (so minus)
                                bvec_file.writelines('{0:.6f}'.format(-x) + ' '
                                    for x in bvec[volumeRange,1])
                                bvec_file.write('\n')
                                bvec_file.writelines('{0:.6f}'.format(x) + ' '
                                    for x in bvec[volumeRange,2])
                        except OSError:
                            logger.error('File name for bvec.txt invalid')
                        #%%%%%%%%%%%%%%%%
                    nr_bvalues = np.unique(
                        iSlices_sorted.diffusion_b_value_number).shape[0]
                    #Always use 4D - Guessing here on "type"
                    expnrslices = (Dim[2] * nr_echos * max(nr_mrtypes,
                        nr_realmrtypes) * nr_diffgrads * nr_dyn)
                    bytes_required = expnrslices * Dim[0] * Dim[1] * 8
                    logger.warning("ND data conversion to nii may require "
                        "huge amounts of memory (RAM)! Your mileage: {0} "
                        "bytes required".format(bytes_required))
                    VolData = np.zeros((Dim[0], Dim[1], expnrslices))
                    slicenr = 0
                    nRowSlices = iSlices_sorted.shape[0]
                    type_map = {8: 'b', 16: 'h', 32: 'i'}
                    while nLine < nRowSlices:
                        sl = iSlices_sorted[nLine]
                        ID1.seek(sl.index_in_rec_file * cDim[0] * cDim[1] * 2)
                        cDim = np.array([sl.recon_resolution_x,
                            sl.recon_resolution_y, Dim[2]])
                        SliceData = array.array(type_map[Parameters.bit])
                        SliceData.fromfile(ID1, cDim[0] * cDim[1])
                        ImageSlice = np.reshape(SliceData, cDim[0:2]).T
                        if not Parameters.issue:
                            VolData[:,:,slicenr] = ImageSlice
                        else:
                            VolData[:,:,slicenr] = ((ImageSlice
                                * Parameters.rescale_slope[slicenr]
                                + Parameters.rescale_interc[slicenr])
                                / (Parameters.scale_slope[slicenr]
                                * Parameters.rescale_slope[slicenr]))
                        slicenr += 1
                        if (nLine + 1) % nRowSlices == 0:
                            if nr_dyn > 1:
                                dyn_suffix = "-{0:04d}".format(
                                    sl.dynamic_scan_number)
                                dyn_ndsuffix = "-d{0:04d}".format(nr_dyn)
                            else:
                                dyn_suffix = "-{0:04d}".format(1)
                                dyn_ndsuffix = ""
                            if nr_mrtypes > 1:
                                mrtype_suffix = "-s{0:3d}".format(
                                    sl.scanning_sequence)
                                mrtype_ndsuffix = "-s{0:03d}".format(nr_mrtypes)
                            else:
                                mrtype_suffix = ""
                                mrtype_ndsuffix = ""
                            if nr_realmrtypes > 1:
                                realmrtype_suffix = "-t{0:03d}".format(
                                    sl.image_type_mr)
                                realmrtype_ndsuffix = "-t{0:03d}".format(
                                    nr_realmrtypes)
                            else:
                                realmrtype_suffix = ""
                                realmrtype_ndsuffix = ""
                            if nr_echos > 1:
                                echo_suffix = "-e{0:03d}".format(
                                    sl.echo_number)
                                echo_ndsuffix = "-e{0:03d}".format(nr_echos)
                            else:
                                echo_suffix = ""
                                echo_ndsuffix = ""
                            if nr_diffgrads > 1:
                                diffgrad_suffix = "-g{0:03d}".format(
                                    sl.gradient_orientation_number)
                                diffgrad_ndsuffix = "-g{0:03d}".format(
                                    nr_diffgrads)
                            else:
                                diffgrad_suffix = ""
                                diffgrad_ndsuffix = ""
                            if nr_bvalues > 1:
                                bval_suffix = "-b{0:03d}".format(
                                    sl.diffusion_b_value_number)
                                bval_ndsuffix = "-b{0:03d}".format(nr_bvalues)
                            else:
                                bval_suffix = ""
                                bval_ndsuffix = ""
                            nii_hdr_dim = np.array([nr_dyn * nr_diffgrads *
                                nr_echos * max(nr_mrtypes, nr_realmrtypes)])
                            NHdr = CreateNiiHdr(Parameters, angulation,
                                nii_hdr_dim)
                            if dti_revertb0:
                                #Do a very dirty dti hack, putting b0 first
                                #(UMCU only)
                                logger.warning("Doing a very dirty hack "
                                    "to DTI data, putting b0 data as first "
                                    "in ND nii file. Only use this "
                                    "(hidden) raw2nii option when working "
                                    "at the UMC Utrecht with Rene Mandl's "
                                    "DTI tools")
                                nslice = Parameters.dim[2]
                                #Last shall be first and first shall be last
                                VolData = np.roll(VolData, nslice, 2)
                            VolNameSinExt = os.path.join(absDir, outFileName)
                            #Construct 4D transformation matrix from T
                            #(translation) and Zm (zoom, from voxel size), and
                            #rotation. Write data to img or nii file
                            VolName = VolNameSinExt + '.nii'
                            logger.info("Writing file: {0}...".format(VolName))
                            WriteNii(VolName, NHdr, VolData, Parameters.issue)
                            logger.info("  ...done")
                            outfiles.append((VolName, VolData))
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
    logger.setLevel(logging.INFO)
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
    parser.add_argument("--dti_revertb0", type=int, default=0)
    parser.add_argument("--debug", "-d", action="store_true")
    options = parser.parse_args()
    if options.debug:
        logger.setLevel(logging.DEBUG)
    convert_raw2nii(options.filelist, options.prefix, options.suffix,
        options.pathpar, options.outfolder, options.outputformat,
        options.angulation, options.rescale, options.dti_revertb0)
