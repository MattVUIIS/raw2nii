LABEL_TYPE_ENUM = {
	'-': 1,
	'CONTROL': 0,
	'LABEL': 1,
	'USERDEF2': 2,
	'USERDEF3': 3,
	'USERDEF4': 4,
	'USERDEF5': 5,
	'USERDEF6': 6,
	'USERDEF7': 7,
	'USERDEF8': 8,
	'USERDEF9': 9,
}

IMG_TYPE_ENUM = {
    'M': 0,
	'R': 1,
	'I': 2,
	'P': 3,
	'CR': 4,
	'T0': 5,
	'T1': 6,
	'T2': 7,
	'RHO': 8,
	'SPECTRO': 9,
	'DERIVED': 10,
	'ADC': 11,
	'RCBV': 12,
	'RCBF': 13,
	'MTT': 14,
	'TTP': 15,
	'FA': 16,
	'EADC': 17,
	'B0': 18,
	'DELAY': 19,
	'MAXRELENH': 20,
	'RELENH': 21,
	'MAXENH': 22,
	'WASHIN': 23,
	'WASHOUT': 24,
	'BREVENH': 25,
	'AREACURV': 26,
	'ANATOMIC': 27,
	'T_TEST': 28,
	'STD_DEVIATION': 29,
	'PERFUSION': 30,
	'T2_STAR': 31,
	'R2': 32,
	'R2_STAR': 33,
	'W': 34,
	'IP': 35,
	'OP': 36,
	'F': 37,
	'SPARE1': 38,
	'SPARE2': 39,
}

IMG_SEQ_ENUM = {
	'IR': 0,
	'SE': 1,
	'FFE': 2,
	'DERIVED': 3,
	'PCA': 4,
	'UNSPECIFIED': 5,
	'SPECTRO': 6,
	'SI': 7,
	'B0': 8,
	'B1': 9,
}

DISPLAY_ORIENT_ENUM = {
	'NONE': 0,
	'RIGHT90': 1,
	'RIGHT180': 2,
	'LEFT90': 3,
	'VM': 4,
	'RIGHT90VM': 5,
	'RIGHT180VM': 6,
	'LEFT90VM': 7,
}

IMG_TYPE_ED_ES_ENUM = {
	'Ed': 0, # ED
	'Es': 1, # ES
	'U': 2, # UNKNOWN
}
"""
SLICE_ORIENT_ENUM = {
	'Undefined': 0,
	'Transversal': 1,
	'Sagittal': 2,
	'Coronal': 3,
	'UNDEFINED': 0,
	'TRANSVERSAL': 1,
	'SAGITTAL': 2,
	'CORONAL': 3,
}
"""
SLICE_ORIENT_ENUM = {
	'FH': 1,
	'RL': 2,
	'AP': 3,
}

CONTRAST_TYPE_ENUM = {
	'DIFFUSION': 0,
	'FLOW_ENCODED': 1,
	'FLUID_ATTENUATED': 2,
	'PERFUSION': 3,
	'PROTON_DENSITY': 4,
	'STIR': 5,
	'TAGGING': 6,
	'T1': 7,
	'T2': 8,
	'T2_STAR': 9,
	'TOF': 10,
	'UNKNOWN': 11,
	'MIXED': 12,
}

DIFF_ANISOTROPY_TYPE_ENUM = {
	'-': 0,
	'': '-',
	'FRACTIONAL': 0,
	'RELATIVE': 1,
	'VOLUME_RATIO': 2,
}

SERIES_DATA_TYPE_ENUM = {
	'PIXEL': 'MRSERIES',
	'RAW': 'UNDEFINED',
	'COMPLEX': 'UNDEFINED',
	'SPECTRA': 'UNDEFINED',
}

PATIENT_POSITION_ENUM = {
	'HFP': 'HeadFirstProne',
	'HFS': 'HeadFirstSupine',
	'HFDR': 'HeadFirstDecubitusRight',
	'HFDL': 'HeadFirstDecubitusLeft',
	'FFP': 'FeetFirstProne',
	'FFS': 'FeetFirstSupine',
	'FFDR': 'FeetFirstDecubitusRight',
	'FFDL': 'FeetFirstDecubitusLeft',
}

PREPARATION_DIR_ENUM = {
	'AP': 'Anterior-Posterior',
	'RL': 'Right-Left',
	'FH': 'Feet-Head',
}

BOOLEAN_ENUM = {
	'Y': 1,
	'N': 0,
}

PAR_HEADER = """# === DATA DESCRIPTION FILE ======================================================
#
# CAUTION - Investigational device.
# Limited by Federal Law to investigational use.
#
# Dataset name: {dataset_name} (raw2nii {tool_version})
#
# CLINICAL TRYOUT             Research image export tool     {par_version}
#
# === GENERAL INFORMATION ========================================================
#
"""

PAR_GEN_INFO = """.    Patient name                       :   {patient_name}
.    Examination name                   :   {exam_name}
.    Protocol name                      :   {protocol_name}
.    Examination date/time              :   {exam_date} / {exam_time}
.    Series Type                        :   Image   {series_type}
.    Acquisition nr                     :   {acquisition_nr}
.    Reconstruction nr                  :   {recon_nr}
.    Scan Duration [sec]                :   {scan_duration}
.    Max. number of cardiac phases      :   {max_n_cardiac_phases}
.    Max. number of echoes              :   {max_n_echoes}
.    Max. number of slices/locations    :   {max_n_slices}
.    Max. number of dynamics            :   {max_n_dynamics}
.    Max. number of mixes               :   {max_n_mixes}
.    Patient position                   :   {patient_pos}
.    Preparation direction              :   {preparation_dir}
.    Technique                          :   {technique}
.    Scan resolution  (x, y)            :   {scan_res_x}  {scan_res_y}
.    Scan mode                          :   {scan_mode}
.    Repetition time [ms]               :   {rep_time}
.    FOV (ap,fh,rl) [mm]                :   {fov_ap}  {fov_fh}  {fov_rl}
.    Water Fat shift [pixels]           :   {water_fat_shift}
.    Angulation midslice(ap,fh,rl)[degr]:   {ang_midslice_ap}  {ang_midslice_fh}  {ang_midslice_rl}
.    Off Centre midslice(ap,fh,rl) [mm] :   {offcenter_midslice_ap}  {offcenter_midslice_fh}  {offcenter_midslice_rl}
.    Flow compensation <0=no 1=yes> ?   :   {flow_compensation}
.    Presaturation     <0=no 1=yes> ?   :   {presaturation}
.    Phase encoding velocity [cm/sec]   :   {phase_encoding_velocity_0}  {phase_encoding_velocity_1}  {phase_encoding_velocity_2}
.    MTC               <0=no 1=yes> ?   :   {mtc}
.    SPIR              <0=no 1=yes> ?   :   {spir}
.    EPI factor        <0,1=no EPI>     :   {epi_factor}
.    Dynamic scan      <0=no 1=yes> ?   :   {dynamic_scan}
.    Diffusion         <0=no 1=yes> ?   :   {diffusion}
.    Diffusion echo time [ms]           :   {diff_echo_time}
.    Max. number of diffusion values    :   {max_n_diff_values}
.    Max. number of gradient orients    :   {max_n_grad_orients}
.    Number of label types   <0=no ASL> :   {n_label_types}
"""

PAR_MIDDLE_SECTION = """#
# === PIXEL VALUES =============================================================
#  PV = pixel value in REC file, FP = floating point value, DV = displayed value on console
#  RS = rescale slope,           RI = rescale intercept,    SS = scale slope
#  DV = PV * RS + RI             FP = DV / (RS * SS)
#
# === IMAGE INFORMATION DEFINITION =============================================
#  The rest of this file contains ONE line per image, this line contains the following information:
#
#  slice number                             (integer)
#  echo number                              (integer)
#  dynamic scan number                      (integer)
#  cardiac phase number                     (integer)
#  image_type_mr                            (integer)
#  scanning sequence                        (integer)
#  index in REC file (in images)            (integer)
#  image pixel size (in bits)               (integer)
#  scan percentage                          (integer)
#  recon resolution (x y)                   (2*integer)
#  rescale intercept                        (float)
#  rescale slope                            (float)
#  scale slope                              (float)
#  window center                            (integer)
#  window width                             (integer)
#  image angulation (ap,fh,rl in degrees )  (3*float)
#  image offcentre (ap,fh,rl in mm )        (3*float)
#  slice thickness (in mm )                 (float)
#  slice gap (in mm )                       (float)
#  image_display_orientation                (integer)
#  slice orientation ( TRA/SAG/COR )        (integer)
#  fmri_status_indication                   (integer)
#  image_type_ed_es  (end diast/end syst)   (integer)
#  pixel spacing (x,y) (in mm)              (2*float)
#  echo_time                                (float)
#  dyn_scan_begin_time                      (float)
#  trigger_time                             (float)
#  diffusion_b_factor                       (float)
#  number of averages                       (integer)
#  image_flip_angle (in degrees)            (float)
#  cardiac frequency   (bpm)                (integer)
#  minimum RR-interval (in ms)              (integer)
#  maximum RR-interval (in ms)              (integer)
#  TURBO factor  <0=no turbo>               (integer)
#  Inversion delay (in ms)                  (float)
#  diffusion b value number    (imagekey!)  (integer)
#  gradient orientation number (imagekey!)  (integer)
#  contrast type                            (string)
#  diffusion anisotropy type                (string)
#  diffusion (ap, fh, rl)                   (3*float)
#  label type (ASL)            (imagekey!)  (integer)
#
# === IMAGE INFORMATION ==========================================================
#  sl ec  dyn ph ty    idx pix scan% rec size                (re)scale              window        angulation              offcentre        thick   gap   info      spacing     echo     dtime   ttime    diff  avg  flip    freq   RR-int  turbo delay b grad cont anis         diffusion       L.ty

"""

PAR_FOOTER = """
# === END OF DATA DESCRIPTION FILE ===============================================
"""
