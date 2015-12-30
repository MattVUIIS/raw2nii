from __future__ import division
import nibabel
import nibabel.dft
import numpy as np
import pprint
from decimal import Decimal

from DICOMFile import DICOMFile, ImageFrame, MRFrame

__all__ = ['read_dicom']


def _read_frame(dcm, index, sl):
    frame = ImageFrame()
    #Main frame data
    #Phillips-specific tag data is in tag 0x2005, 0x140f
    private_tag = sl[(0x2005, 0x140f)][0]
    frame.slice_num = np.int(private_tag[(0x2001, 0x100a)].value)
    frame.echo_num = np.int(private_tag[(0x0018, 0x0086)].value)
    frame.dynamic = np.int(private_tag[(0x0020, 0x0100)].value)
    frame.cardiac_phase = np.int(private_tag[(0x2001, 0x1008)].value)
    frame.bvalue = np.int(private_tag[(0x2005, 0x1412)].value)
    frame.grad_orient = np.int(private_tag[(0x2005, 0x1413)].value)
    #If label type is not present, default to 1
    label_type = private_tag.get((0x2005, 0x1429))
    if label_type:
        frame.label_type = label_type.value
    else:
        frame.label_type = '-'
    frame.img_type = private_tag[(0x2005, 0x1011)].value
    frame.img_seq = private_tag[(0x2005, 0x106e)].value
    frame.index = index
    frame.px_size = dcm.px_size
    frame.scan_prct = Decimal(sl.MRFOVGeometrySequence[0].PercentSampling)
    frame.res_x = dcm.res_x
    frame.res_y = dcm.res_y
    pvts = sl.PixelValueTransformationSequence[0]
    frame.rescale_interc = Decimal(pvts.RescaleIntercept)
    frame.rescale_slope = Decimal(pvts.RescaleSlope)
    frame.scale_slope = Decimal(private_tag[(0x2005, 0x100e)].value)
    fvls = sl.FrameVOILUTSequence[0]
    frame.window_center = np.int(fvls.WindowCenter)
    frame.window_width = np.int(fvls.WindowWidth)
    pms = sl.PixelMeasuresSequence[0]
    frame.sl_thickness = Decimal(pms.SliceThickness)
    space_bt_slices = Decimal(private_tag[(0x0018, 0x0088)].value)
    frame.sl_gap = space_bt_slices - frame.sl_thickness
    frame.dsply_orient = private_tag[(0x2005, 0x1004)].value
    frame.fmri_stat_indication = np.int(private_tag[(0x2005, 0x1063)].value)
    frame.img_type_ed_es = private_tag[(0x2001, 0x1007)].value
    frame.px_spacing = np.array(pms.PixelSpacing, Decimal)
    frame.echo_t = Decimal(private_tag[(0x0018, 0x0081)].value)
    frame.dyn_scan_begin_t = Decimal(private_tag[(0x2005, 0x10a0)].value)
    #If trigger time is not present, default to 0.0
    trigger_t = private_tag.get((0x0018, 0x1060))
    if trigger_t:
        frame.trigger_t = Decimal(trigger_t.value)
    else:
        frame.trigger_t = 0.0
    frame.diff_b_factor = Decimal(private_tag[(0x2001, 0x1003)].value)
    frame.no_avgs = np.int(private_tag[(0x0018, 0x0083)].value)
    frame.img_flip_ang = Decimal(private_tag[(0x0018, 0x1314)].value)
    frame.cardiac_freq = np.int(private_tag[(0x0018, 0x1088)].value)
    frame.min_rr_interval = np.int(private_tag[(0x0018, 0x1081)].value)
    frame.max_rr_interval = np.int(private_tag[(0x0018, 0x1082)].value)
    frame.turbo_factor = np.int(private_tag[(0x0018, 0x0091)].value)
    frame.inversion_delay = Decimal(private_tag[(0x0018, 0x0082)].value)
    frame.diff_anisotropy_type = private_tag[(0x0018, 0x9147)].value
    frame.contrast_type = sl.MRImageFrameTypeSequence[0].AcquisitionContrast
    frame.diff_ap = Decimal(private_tag[(0x2005, 0x10b1)].value)
    frame.diff_fh = Decimal(private_tag[(0x2005, 0x10b2)].value)
    frame.diff_rl = Decimal(private_tag[(0x2005, 0x10b0)].value)
    fcs = sl.FrameContentSequence[0]
    frame.stack_id = np.int(fcs[(0x0020, 0x9056)].value) - 1
    mr_frame = dcm.mr_stack[frame.stack_id]
    frame.ang_ap = mr_frame.ang_ap
    frame.ang_fh = mr_frame.ang_fh
    frame.ang_rl = mr_frame.ang_rl
    img_pos = np.array(sl.PlanePositionSequence[0].ImagePositionPatient,
        Decimal)
    img_orient = np.array(sl.PlaneOrientationSequence[0]
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
    return frame

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
    dcm.series_data_type = ds[(0x2005, 0x1035)].value
    dcm.acquisition_nr = np.int(ds.AcquisitionNumber)
    dcm.recon_nr = np.int(ds[(0x2001, 0x101d)].value)
    dcm.acquisition_dur = Decimal(ds.AcquisitionDuration)
    dcm.max_n_phases_mr = np.int(ds[(0x2001, 0x1017)].value)
    dcm.max_n_echoes = np.int(ds[(0x2001, 0x1014)].value)
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
    #Read MR stack
    dcm.mr_stack = []
    for item in ds[(0x2001, 0x105f)]:
        mr = MRFrame()
        mr.preparation_dir = item[0x2005107b].value
        mr.fov_ap = Decimal(item[0x20051074].value)
        mr.fov_fh = Decimal(item[0x20051075].value)
        mr.fov_rl = Decimal(item[0x20051076].value)
        mr.ang_ap = Decimal(item[(0x2005, 0x1071)].value)
        mr.ang_fh = Decimal(item[(0x2005, 0x1072)].value)
        mr.ang_rl = Decimal(item[(0x2005, 0x1073)].value)
        mr.offcenter_ap = Decimal(item[0x20051078].value)
        mr.offcenter_fh = Decimal(item[0x20051079].value)
        mr.offcenter_rl = Decimal(item[0x2005107a].value)
        mr.view_axis = item[(0x2005, 0x1081)].value
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
    dcm.frames = []
    for i, field in enumerate(ds.PerFrameFunctionalGroupsSequence):
        frame = _read_frame(dcm, i, field)
        dcm.frames.append(frame)
    dcm.raw_data = ds[(0x7fe00010)].value
    return dcm
