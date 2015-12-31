from __future__ import division
import nibabel
import nibabel.dft
import numpy as np
import pprint
from decimal import Decimal

from DICOMFile import DICOMFile, ImageFrame, MRFrame

__all__ = ['read_dicom']


def _read_functional_group(dcm, fg, shared_fgs):
    """ Reads a functional group from the DICOM file and returns a slice/frame
        dcm: DICOMFile
        fg: Per-Frame Functional Group
        shared_fgs: Shared Functional Group Sequence
    """
    frame = ImageFrame()
    #Main frame data
    #Phillips-specific tag data is in tag 0x2005, 0x140f
    private_tag = fg[0x2005140f][0]
    frame.slice_num = np.int(private_tag[0x2001100a].value)
    frame.echo_num = np.int(private_tag[0x00180086].value)
    frame.dynamic = np.int(private_tag[0x00200100].value)
    frame.cardiac_phase = np.int(private_tag[0x20011008].value)
    frame.bvalue = np.int(private_tag[0x20051412].value)
    frame.grad_orient = np.int(private_tag[0x20051413].value)
    #If label type is not present, default to 1
    label_type = private_tag.get(0x20051429)
    if label_type:
        frame.label_type = label_type.value
    else:
        frame.label_type = '-'
    frame.img_type = private_tag[0x20051011].value
    frame.img_seq = private_tag[0x2005106e].value
    frame.px_size = dcm.px_size

    if hasattr(fg, 'MRFOVGeometrySequence'):
        frame.scan_prct = Decimal(fg.MRFOVGeometrySequence[0].PercentSampling)
    else:
        frame.scan_prct = Decimal(shared_fgs[0].MRFOVGeometrySequence[0].PercentSampling)

    frame.res_x = dcm.res_x
    frame.res_y = dcm.res_y
    pvts = fg.PixelValueTransformationSequence[0]
    frame.rescale_interc = Decimal(pvts.RescaleIntercept)
    frame.rescale_slope = Decimal(pvts.RescaleSlope)
    frame.scale_slope = Decimal(private_tag[0x2005100e].value)
    fvls = fg.FrameVOILUTSequence[0]
    frame.window_center = np.int(fvls.WindowCenter)
    frame.window_width = np.int(fvls.WindowWidth)

    pms = fg.PixelMeasuresSequence[0]
    if hasattr(pms, 'SliceThickness'):
        frame.sl_thickness = Decimal(pms.SliceThickness)
    else:
        frame.sl_thickness = Decimal(private_tag.SliceThickness)

    frame.sl_gap = dcm.space_bt_slices - frame.sl_thickness
    frame.dsply_orient = private_tag[0x20051004].value
    frame.fmri_stat_indication = np.int(private_tag[0x20051063].value)
    frame.img_type_ed_es = private_tag[0x20011007].value
    frame.px_spacing = np.array(pms.PixelSpacing, Decimal)
    frame.echo_t = Decimal(private_tag[0x00180081].value)
    frame.dyn_scan_begin_t = Decimal(private_tag[0x200510a0].value)
    #If trigger time is not present, default to 0.0
    trigger_t = private_tag.get(0x00181060)
    if trigger_t:
        frame.trigger_t = Decimal(trigger_t.value)
    else:
        frame.trigger_t = 0.0
    frame.diff_b_factor = Decimal(private_tag[0x20011003].value)
    frame.no_avgs = np.int(private_tag[0x00180083].value)
    frame.img_flip_ang = Decimal(private_tag[0x00181314].value)
    frame.cardiac_freq = np.int(private_tag[0x00181088].value)
    frame.min_rr_interval = np.int(private_tag[0x00181081].value)
    frame.max_rr_interval = np.int(private_tag[0x00181082].value)
    frame.turbo_factor = np.int(private_tag[0x00180091].value)
    frame.inversion_delay = Decimal(private_tag[0x00180082].value)
    frame.diff_anisotropy_type = private_tag[0x00189147].value
    frame.contrast_type = fg.MRImageFrameTypeSequence[0].AcquisitionContrast
    frame.diff_ap = Decimal(private_tag[0x200510b1].value)
    frame.diff_fh = Decimal(private_tag[0x200510b2].value)
    frame.diff_rl = Decimal(private_tag[0x200510b0].value)
    fcs = fg.FrameContentSequence[0]
    frame.mr_stack_id = np.int(fcs[0x00209056].value) - 1
    mr_frame = dcm.mr_stack[frame.mr_stack_id]
    frame.ang_ap = mr_frame.ang_ap
    frame.ang_fh = mr_frame.ang_fh
    frame.ang_rl = mr_frame.ang_rl
    img_pos = np.array(fg.PlanePositionSequence[0].ImagePositionPatient,
        Decimal)
    img_orient = np.array(fg.PlaneOrientationSequence[0]
        .ImageOrientationPatient, Decimal)
    frame.offcenter_ap = img_pos[1] + (
        (dcm.res_y * frame.px_spacing[0] * img_orient[1]) +
        (dcm.res_x * frame.px_spacing[1] * img_orient[4])) / 2
    frame.offcenter_fh = img_pos[2] + (
        (dcm.res_y * frame.px_spacing[0] * img_orient[2]) +
        (dcm.res_x * frame.px_spacing[1] * img_orient[5])) / 2
    frame.offcenter_rl = img_pos[0] + (
        (dcm.res_y * frame.px_spacing[0] * img_orient[0]) +
        (dcm.res_x * frame.px_spacing[1] * img_orient[3])) / 2
    frame.slice_orient = mr_frame.view_axis
    frame.temporal_pos_index = np.int(fcs[0x00209128].value) - 1
    dcm.stacks[frame.temporal_pos_index].append(frame)

"""
$series_info_special_hash{"MRSeriesNrOfStacks"}{"Tag"}           = "0x20011060";
$series_info_special_hash{"MFNumberOfFrames"}{"Tag"}             = "0x00280008";

my %stack_info_special_hash = ();
$stack_info_special_hash{"MRStackAngulationAP"}{"Tag"}          = "0x20051071";
$stack_info_special_hash{"MRStackAngulationFH"}{"Tag"}          = "0x20051072";
$stack_info_special_hash{"MRStackAngulationRL"}{"Tag"}          = "0x20051073";
$stack_info_special_hash{"MRStackViewAxis"}{"Tag"}              = "0x20051081";

my %frame_info_special_hash = ();
$frame_info_special_hash{"MRImageSpacingBetweenSlices"}{"Tag"}  = "0x00180088";
$frame_info_special_hash{"StackID"}{"Tag"}                      = "0x00209056";
$frame_info_special_hash{"ImagePlanePositionPatient"}{"Tag"}    = "0x00200032";
$frame_info_special_hash{"ImagePlaneOrientationPatient"}{"Tag"} = "0x00200037";
$frame_info_special_hash{"MFInstanceNumber"}{"Tag"}             = "0x00200013";
"""

def read_dicom(dcm_fname):
    dcm = DICOMFile()
    ds = nibabel.dft.dicom.read_file(dcm_fname)
    dcm.patient_name = ds.PatientName
    dcm.exam_name = ds.PerformedProcedureStepDescription
    dcm.protocol_name = ds.ProtocolName
    dcm.exam_date = ds.PerformedProcedureStepStartDate
    dcm.exam_time = ds.PerformedProcedureStepStartTime
    dcm.series_data_type = ds[0x20051035].value
    dcm.acquisition_nr = np.int(ds.AcquisitionNumber)
    dcm.recon_nr = np.int(ds[0x2001101d].value)
    dcm.acquisition_dur = Decimal(ds.AcquisitionDuration)
    dcm.max_n_phases_mr = np.int(ds[0x20011017].value)
    dcm.max_n_echoes = np.int(ds[0x20011014].value)
    dcm.max_n_slices = np.int(ds[0x20011018].value)
    dcm.max_n_dyn = np.int(ds[0x20011081].value)
    dcm.max_n_mixes = np.int(ds[0x20051021].value)
    dcm.max_n_bvalues = np.int(ds[0x20051414].value)
    dcm.max_n_grad_orients = np.int(ds[0x20051415].value)
    dcm.num_label_types = np.int(ds[0x20051428].value)
    dcm.patient_pos = ds.PatientPosition
    dcm.technique = ds[0x20011020].value
    dcm.scan_res_x = np.int(ds[0x2005101d].value)
    sfgs = ds.SharedFunctionalGroupsSequence[0]
    dcm.scan_res_y = np.int(sfgs[0x2005140e][0].NumberOfPhaseEncodingSteps)
    dcm.scan_mode = ds[0x2005106f].value
    dcm.rep_time = Decimal(ds[0x20051030].value)
    dcm.water_fat_shift = Decimal(ds[0x20011022].value)
    dcm.flow_compensation = ds[0x20051016].value
    dcm.presaturation = ds[0x2005102f].value
    dcm.phase_encoding_velocity = np.array(ds[0x2001101a].value, Decimal)
    dcm.mtc = ds[0x2005101c].value
    dcm.spir = ds[0x20011021].value
    dcm.epi_factor = np.int(ds[0x20011013].value)
    dcm.dyn_scan = ds[0x20011012].value
    dcm.diffusion = ds[0x20051014].value
    dcm.diff_echo_time = ds[0x20011011].value
    dcm.study_id = ds.StudyID
    dcm.patient_id = ds.PatientID
    dcm.study_inst_uid = ds.StudyInstanceUID
    dcm.series_inst_uid = ds.SeriesInstanceUID
    dcm.sop_class_uid = ds.SOPClassUID
    dcm.series_date = ds.SeriesDate
    dcm.series_time = ds.SeriesTime
    dcm.inst_creation_date = ds.InstanceCreationDate
    dcm.inst_creation_time = ds.InstanceCreationTime
    dcm.mr_series_development_mode = ds[0x20051013].value
    dcm.data_dict_contents_version = ds[0x2005143a].value
    #Shared data
    dcm.px_size = ds.BitsAllocated
    dcm.res_x = ds.Columns
    dcm.res_y = ds.Rows
    dcm.space_bt_slices = Decimal(ds[0x00180088].value)
    #Read MR stack
    dcm.mr_stack = []
    for item in ds[0x2001105f]:
        mr = MRFrame()
        mr.preparation_dir = item[0x2005107b].value
        mr.fov_ap = Decimal(item[0x20051074].value)
        mr.fov_fh = Decimal(item[0x20051075].value)
        mr.fov_rl = Decimal(item[0x20051076].value)
        mr.ang_ap = Decimal(item[0x20051071].value)
        mr.ang_fh = Decimal(item[0x20051072].value)
        mr.ang_rl = Decimal(item[0x20051073].value)
        mr.offcenter_ap = Decimal(item[0x20051078].value)
        mr.offcenter_fh = Decimal(item[0x20051079].value)
        mr.offcenter_rl = Decimal(item[0x2005107a].value)
        mr.view_axis = item[0x20051081].value
        mr.num_stack_slices = np.int(item[0x2001102d].value)
        dcm.mr_stack.append(mr)
    last_mr = dcm.mr_stack[-1]
    dcm.preparation_dir = last_mr.preparation_dir
    dcm.fov_ap = last_mr.fov_ap
    dcm.fov_fh = last_mr.fov_fh
    dcm.fov_rl = last_mr.fov_rl
    dcm.ang_ap = last_mr.ang_ap
    dcm.ang_fh = last_mr.ang_fh
    dcm.ang_rl = last_mr.ang_rl
    dcm.offcenter_ap = last_mr.offcenter_ap
    dcm.offcenter_fh = last_mr.offcenter_fh
    dcm.offcenter_rl = last_mr.offcenter_rl
    shared_fgs = ds.SharedFunctionalGroupsSequence
    num_stack_slices = last_mr.num_stack_slices
    num_fg = len(ds.PerFrameFunctionalGroupsSequence)
    dcm.num_stacks = num_fg // num_stack_slices
    dcm.stacks = [[] for i in range(dcm.num_stacks)]
    for fg in ds.PerFrameFunctionalGroupsSequence:
        _read_functional_group(dcm, fg, shared_fgs)
    dcm.raw_data = ds[(0x7fe00010)].value
    return dcm
