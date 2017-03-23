# raw2nii
An open tool for converting Phillips PAR/REC files into the NIfTI format.
Dependencies: numpy, nibabel

### Usage
From the project folder:
```bash
./raw2nii.py img.PAR img.NII
```
It also supports DICOM to PAR v4.2 conversion:
```bash
./raw2nii.py img.DCM img.PAR
```
DICOM to PARREC conversion is still in the experimental phase. Don't rely on it
for any purpose other than testing.
