from __future__ import division
import numpy as np
import pprint
import re
from decimal import Decimal

import d2p_defines
import raw2nii_version


__all__ = ['write_parrec_from_dicom']

def _get_header(dataset_name):
    info = {
        'dataset_name': dataset_name,
        'tool_version': raw2nii_version.VERSION,
        'par_version': 'V4.2',
    }
    return d2p_defines.PAR_HEADER.format(**info)

def _get_general_info(dcm):
    gen_info = {
#.    Patient name                       :   {patient_name}
        'patient_name': dcm.patient_name,
#.    Examination name                   :   {exam_name}
        'exam_name': dcm.exam_name,
#.    Protocol name                      :   {protocol_name}
        'protocol_name': dcm.protocol_name,
#.    Examination date/time              :   {exam_datetime}
        'exam_date': '{0}.{1}.{2}'.format(dcm.exam_date[:4],
            dcm.exam_date[4:6], dcm.exam_date[6:8]),
        'exam_time': '{0}:{1}:{2}'.format(dcm.exam_time[:2],
            dcm.exam_time[2:4], dcm.exam_time[4:6]),
#.    Series Type                        :   Image   {series_type}
        'series_type': d2p_defines.SERIES_DATA_TYPE_ENUM[dcm.series_data_type],
#.    Acquisition nr                     :   {acquisition_nr}
        'acquisition_nr': dcm.acquisition_nr,
#.    Reconstruction nr                  :   {recon_nr}
        'recon_nr': dcm.recon_nr,
#.    Scan Duration [sec]                :   {scan_duration}
        'scan_duration': '{0:.1f}'.format(dcm.acquisition_dur),
#.    Max. number of cardiac phases      :   {max_n_cardiac_phases}
        'max_n_cardiac_phases': dcm.max_n_phases_mr,
#.    Max. number of echoes              :   {max_n_echoes}
        'max_n_echoes': dcm.max_n_echoes,
#.    Max. number of slices/locations    :   {max_n_slices}
        'max_n_slices': dcm.max_n_slices,
#.    Max. number of dynamics            :   {max_n_dynamics}
        'max_n_dynamics': dcm.max_n_dyn,
#.    Max. number of mixes               :   {max_n_mixes}
        'max_n_mixes': dcm.max_n_mixes,
#.    Patient position                   :   {patient_pos}
        'patient_pos': d2p_defines.PATIENT_POSITION_ENUM[dcm.patient_pos],
#.    Preparation direction              :   {preparation_dir}
        'preparation_dir': d2p_defines.PREPARATION_DIR_ENUM[
            dcm.preparation_dir],
#.    Technique                          :   {technique}
        'technique': dcm.technique,
#.    Scan resolution  (x, y)            :   {scan_res_x}  {scan_res_y}
        'scan_res_x': dcm.scan_res_x,
        'scan_res_y': dcm.scan_res_y,
#.    Scan mode                          :   {scan_mode}
        'scan_mode': dcm.scan_mode,
#.    Repetition time [ms]               :   {rep_time}
        'rep_time': '{0:.3f}'.format(dcm.rep_time),
#.    FOV (ap,fh,rl) [mm]                :   {fov_ap}  {fov_fh}  {fov_rl}
        'fov_ap': '{0:.3f}'.format(dcm.fov_ap),
        'fov_fh': '{0:.3f}'.format(dcm.fov_fh),
        'fov_rl': '{0:.3f}'.format(dcm.fov_rl),
#.    Water Fat shift [pixels]           :   {water_fat_shift}
        'water_fat_shift': '{0:.3f}'.format(dcm.water_fat_shift),
#.    Angulation midslice(ap,fh,rl)[degr]:   {ang_midslice_ap}  {ang_midslice_fh}  {ang_midslice_rl}
        'ang_midslice_ap': '{0:.3f}'.format(dcm.ang_ap),
        'ang_midslice_fh': '{0:.3f}'.format(dcm.ang_fh),
        'ang_midslice_rl': '{0:.3f}'.format(dcm.ang_rl),
#.    Off Centre midslice(ap,fh,rl) [mm] :   {offcenter_midslice_ap}  {offcenter_midslice_fh}  {offcenter_midslice_rl}
        'offcenter_midslice_ap': '{0:.3f}'.format(dcm.offcenter_ap),
        'offcenter_midslice_fh': '{0:.3f}'.format(dcm.offcenter_fh),
        'offcenter_midslice_rl': '{0:.3f}'.format(dcm.offcenter_rl),
#.    Flow compensation <0=no 1=yes> ?   :   {flow_compensation}
        'flow_compensation': d2p_defines.BOOLEAN_ENUM[dcm.flow_compensation],
#.    Presaturation     <0=no 1=yes> ?   :   {presaturation}
        'presaturation': d2p_defines.BOOLEAN_ENUM[dcm.presaturation],
#.    Phase encoding velocity [cm/sec]   :   {phase_encoding_velocity_0}  {phase_encoding_velocity_1}  {phase_encoding_velocity_2}
        'phase_encoding_velocity_0': '{0:.6f}'.format(dcm.phase_encoding_velocity[0]),
        'phase_encoding_velocity_1': '{0:.6f}'.format(dcm.phase_encoding_velocity[1]),
        'phase_encoding_velocity_2': '{0:.6f}'.format(dcm.phase_encoding_velocity[2]),
#.    MTC               <0=no 1=yes> ?   :   {mtc}
        'mtc': d2p_defines.BOOLEAN_ENUM[dcm.mtc],
#.    SPIR              <0=no 1=yes> ?   :   {spir}
        'spir': d2p_defines.BOOLEAN_ENUM[dcm.spir],
#.    EPI factor        <0,1=no EPI>     :   {epi_factor}
        'epi_factor': dcm.epi_factor,
#.    Dynamic scan      <0=no 1=yes> ?   :   {dynamic_scan}
        'dynamic_scan': d2p_defines.BOOLEAN_ENUM[dcm.dyn_scan],
#.    Diffusion         <0=no 1=yes> ?   :   {diffusion}
        'diffusion': d2p_defines.BOOLEAN_ENUM[dcm.diffusion],
#.    Diffusion echo time [ms]           :   {diff_echo_time}
        'diff_echo_time': '{0:.4f}'.format(dcm.diff_echo_time),
#.    Max. number of diffusion values    :   {max_n_diff_values}
        'max_n_diff_values': dcm.max_n_bvalues,
#.    Max. number of gradient orients    :   {max_n_grad_orients}
        'max_n_grad_orients': dcm.max_n_grad_orients,
#.    Number of label types   <0=no ASL> :   {n_label_types}
        'n_label_types': dcm.num_label_types,
    }
    return d2p_defines.PAR_GEN_INFO.format(**gen_info)


def _get_image_def(frame):
    return ' '.join([
    #  slice number                             (integer)
        '{0:3d}'.format(frame.slice_num),
    #  echo number                              (integer)
        '{0:3d}'.format(frame.echo_num),
    #  dynamic scan number                      (integer)
        '{0:4d}'.format(frame.dynamic),
    #  cardiac phase number                     (integer)
        '{0:2d}'.format(frame.cardiac_phase),
    #  image_type_mr                            (integer)
        str(d2p_defines.IMG_TYPE_ENUM[frame.img_type]),
    #  scanning sequence                        (integer)
        str(d2p_defines.IMG_SEQ_ENUM[frame.img_seq]),
    #  index in REC file (in images)            (integer)
        '{0:5d}'.format(frame.index),
    #  image pixel size (in bits)               (integer)
        '{0:3d}'.format(frame.px_size),
    #  scan percentage                          (integer)
        '{0:5d}'.format(np.int(frame.scan_prct)),
        #str(frame.scan_prct),
    #  recon resolution (x y)                   (2*integer)
        '{0:4d}'.format(frame.res_x), '{0:4d}'.format(frame.res_y),
    #  rescale intercept                        (float)
        '{0:11.5f}'.format(frame.rescale_interc),
    #  rescale slope                            (float)
        '{0:9.5f}'.format(frame.rescale_slope),
    #  scale slope                              (float)
        '{0:.5e}'.format(np.float(frame.scale_slope)),
    #  window center                            (integer)
        '{0:5d}'.format(frame.window_center),
    #  window width                             (integer)
        '{0:5d}'.format(frame.window_width),
    #  image angulation (ap,fh,rl in degrees )  (3*float)
        '{0:6.2f}'.format(frame.ang_ap), '{0:6.2f}'.format(frame.ang_fh),
        '{0:6.2f}'.format(frame.ang_rl),
    #  image offcentre (ap,fh,rl in mm )        (3*float)
        '{0:7.2f}'.format(frame.offcenter_ap),
        '{0:7.2f}'.format(frame.offcenter_fh),
        '{0:7.2f}'.format(frame.offcenter_rl),
    #  slice thickness (in mm )                 (float)
        '{0:.3f}'.format(frame.sl_thickness),
    #  slice gap (in mm )                       (float)
        '{0:.3f}'.format(frame.sl_gap),
    #  image_display_orientation                (integer)
        str(d2p_defines.DISPLAY_ORIENT_ENUM.get(frame.dsply_orient)),
    #  slice orientation ( TRA/SAG/COR )        (integer)
        str(d2p_defines.SLICE_ORIENT_ENUM[frame.slice_orient]),
    #  fmri_status_indication                   (integer)
        str(frame.fmri_stat_indication),
    #  image_type_ed_es  (end diast/end syst)   (integer)
        str(d2p_defines.IMG_TYPE_ED_ES_ENUM[frame.img_type_ed_es]),
    #  pixel spacing (x,y) (in mm)              (2*float)
        '{0:6.3f}'.format(frame.px_spacing[0]),
        '{0:6.3f}'.format(frame.px_spacing[1]),
    #  echo_time                                (float)
        '{0:6.2f}'.format(frame.echo_t),
    #  dyn_scan_begin_time                      (float)
        '{0:7.2f}'.format(frame.dyn_scan_begin_t),
    #  trigger_time                             (float)
        '{0:8.2f}'.format(frame.trigger_t),
    #  diffusion_b_factor                       (float)
        '{0:7.2f}'.format(frame.diff_b_factor),
    #  number of averages                       (integer)
        '{0:3}'.format(frame.no_avgs),
    #  image_flip_angle (in degrees)            (float)
        '{0:7.2f}'.format(frame.img_flip_ang),
    #  cardiac frequency   (bpm)                (integer)
        '{0:5d}'.format(frame.cardiac_freq),
    #  minimum RR-interval (in ms)              (integer)
        '{0:4d}'.format(frame.min_rr_interval),
    #  maximum RR-interval (in ms)              (integer)
        '{0:4d}'.format(frame.max_rr_interval),
    #  TURBO factor  <0=no turbo>               (integer)
        '{0:5d}'.format(frame.turbo_factor),
    #  Inversion delay (in ms)                  (float)
        '{0:5.1f}'.format(frame.inversion_delay),
    #  diffusion b value number    (imagekey!)  (integer)
        '{0:2d}'.format(frame.bvalue),
    #  gradient orientation number (imagekey!)  (integer)
        '{0:3d}'.format(frame.grad_orient),
    #  contrast type                            (string)
        '{0:4d}'.format(d2p_defines.CONTRAST_TYPE_ENUM[frame.contrast_type]),
    #  diffusion anisotropy type                (string)
        d2p_defines.DIFF_ANISOTROPY_TYPE_ENUM[frame.diff_anisotropy_type].rjust(
            4),
    #  diffusion (ap, fh, rl)                   (3*float)
        '{0:7.3f}'.format(frame.diff_ap), '{0:8.3f}'.format(frame.diff_fh),
        '{0:8.3f}'.format(frame.diff_rl),
    #  label type (ASL)            (imagekey!)  (integer)
        '{0:2d}'.format(d2p_defines.LABEL_TYPE_ENUM[frame.label_type]),
    ]) + '\n'

def _sanitize_field_names(*field_names):
    new_fields = list(field_names)
    for i, field in enumerate(field_names):
        field = re.sub(' ', '_', str(field))
        field = re.sub('[\s\\\/\*\?"<>\|]+', '_', field)
        field = re.sub('[:]+', '.', field)
        new_fields[i] = field
    return '_'.join(new_fields)

def write_parrec_from_dicom(par_fname, rec_fname, dcm):
    series_time = '.'.join(dcm.series_time[i:i+2] for i in range(0, 6, 2))
    dataset_name = _sanitize_field_names(dcm.patient_name,
        '{0:02d}'.format(dcm.acquisition_nr), '{0:02d}'.format(dcm.recon_nr),
        series_time, '({0})'.format(dcm.protocol_name))
    with open(par_fname, 'wb') as f:
        f.write(_get_header(dataset_name))
        f.write(_get_general_info(dcm))
        f.write(d2p_defines.PAR_MIDDLE_SECTION)
        f.writelines(_get_image_def(frame) for frame in dcm.frames)
        f.write(d2p_defines.PAR_FOOTER)
    with open(rec_fname, 'wb') as f:
        f.write(dcm.raw_data)
