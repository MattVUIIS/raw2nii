from __future__ import division
import logging
import numpy as np
import os
import sys


class PARInfo:
    def __repr__(self):
        s = []
        for key, val in self.__dict__.items():
            s.append('{0}="{1}"'.format(key, val))
        return "<PARInfo {0}>".format(" ".join(s))

def _skip_lines(file_obj, n):
    for i in range(n):
        file_obj.readline()

def _read_general_info(file_obj):
    """ Reads general info parameters from PAR file and returns a dict """
    logger = logging.getLogger('raw2nii')
    gen_info = {}
    line = None
    while line != "":
        pos = file_obj.tell()
        line = file_obj.readline()
        if ':' not in line:
            file_obj.seek(pos)
            break
        #logger.debug("Line = '{0}'".format(line))
        first, second = line.split(":", 1)
        key = first[1:].strip()
        val = second.strip()
        #logger.debug("Key = '{0}' Val = '{1}'".format(key, val))
        gen_info[key] = val
    return gen_info

def _skip_comment_lines(file_obj):
    logger = logging.getLogger('raw2nii')
    line = None
    while line != "":
        pos = file_obj.tell()
        line = file_obj.readline()
        if not line.startswith('#'):
            file_obj.seek(pos)
            break

def read_par(parfilename):
    """
        parfilename: string with complete par-file name (with path)
        par: MATLAB structure with many important fields

        Philips PAR file interpreter. Reads several important fields from PAR files, and returns them in a structure.
        Reads version number of Philips research tools and interpretes
        accordingly.
        Research tools are used to extract data from database; dataformats differ considerably
        between versions. R2AGUI now handles V3 and V4
    """
    logger = logging.getLogger('raw2nii')
    par = PARInfo()
    with open(parfilename, 'rb') as parfile:
        par.problemreading = 0
        _skip_lines(parfile, 7)  # Skip first 7 lines
        par.ResToolsVersion = parfile.readline().split()[-1]
        logger.debug("ResToolVersion: {0}".format(par.ResToolsVersion))
        if 'V3' == par.ResToolsVersion:
            """
            TODO: Convert this code
            parameter=textread(parfile,char('%s'),5,'delimiter',char(':'),'headerlines',13)
            par.name = copy(parameter[2])
            parameter=textread(parfile,char('%u'),16,'delimiter',char('.:Acquisitionr'),'headerlines',15)
            par.scno = copy(parameter[16])
            parameter=textread(parfile,char('%u'),31,'delimiter',char('.:Max.numberofslices/locations'),'headerlines',20)
            par.slice = copy(parameter[31])
            parameter=textread(parfile,char('%u'),23,'delimiter',char('.:Max.numberofdynamics'),'headerlines',21)
            par.dyn = copy(parameter[23])
            parameter=textread(parfile,char('%s'),25,'delimiter',char('.:Imagepixelsize[orbits]'),'headerlines',23)
            par.bit = copy(parameter[25])
            parameter=textread(parfile,char('%u'),24,'delimiter',char('.:Reconresolution(x,y)'),'headerlines',28)
            x=parameter[23]
            y=parameter[24]
            z=par.slice
            par.dim = copy(cat(x,y,z))
            parameter=textread(parfile,char('%s'),20,'headerlines',90)
            par.sliceorient = copy(parameter[20])
            parameter=textread(parfile,char('%f'),22,'delimiter',char('.:FOV(ap,fh,rl)[mm]'),'headerlines',31)
            if strcmp(par.sliceorient,'1'):
                fovx=parameter[20]
                fovy=parameter[22]
                par.sliceorient
            if strcmp(par.sliceorient,'2'):
                fovx=parameter[21]
                fovy=parameter[20]
                par.sliceorient
            if strcmp(par.sliceorient,'3'):
                fovx=parameter[22]
                fovy=parameter[21]
                par.sliceorient
            par.fov_apfhrl = copy(cat(parameter[20],parameter[21],parameter[22]))
            parameter=textread(parfile,char('%f'),21,'delimiter',char('.:Slicethickness[mm]'),'headerlines',32)
            par.slth = copy(parameter[21])
            parameter=textread(parfile,char('%f'),15,'delimiter',char('.:Slicegap[mm]'),'headerlines',33)
            par.gap = copy(parameter[15])
            fovz=(par.gap + par.slth) * par.slice
            par.fov = copy(cat(fovx,fovy,fovz))
            parameter=textread(parfile,char('%f'),39,'delimiter',char('.:Angulationmidslice(ap,fh,rl)[degr]'),'headerlines',35)
            par.angAP = copy(parameter[37])
            par.angFH = copy(parameter[38])
            par.angRL = copy(parameter[39])
            parameter=textread(parfile,char('%f'),36,'delimiter',char('.:OffCentremidslice(ap,fh,rl)[mm]'),'headerlines',36)
            par.offAP = copy(parameter[34])
            par.offFH = copy(parameter[35])
            par.offRL = copy(parameter[36])
            parameter=textread(parfile,char('%s'),24,'headerlines',88)
            voxx=str2num(parameter[23])
            voxy=str2num(parameter[24])
            voxz=par.slth + par.gap
            par.vox = copy(cat(voxx,voxy,voxz))
            par.rescale_slope = copy(str2num(parameter[9]))
            par.rescale_interc = copy(str2num(parameter[8]))
            parameternextline=textread(parfile,char('%s'),24,'headerlines',90)
            if (parameternextline[1] - parameter[1]) > 0:
                par.slicessorted = copy(1)
            else:
                par.slicessorted = copy(2)
            slice_index=textread(parfile,'','delimiter',char(' '),'headerlines',88,'commentstyle','shell')
            par.RT = copy((slice_index[end(),26] - slice_index[1,26]) / (par.dyn - 1))
            clear('parameter')
            """
        elif par.ResToolsVersion in ('V4', 'V4.1', 'V4.2'):
            _skip_lines(parfile, 5)
            gen_info = _read_general_info(parfile)
            logger.debug("Parameters name: {0}".format(gen_info))
            par.name = gen_info['Protocol name']
            par.scno = int(gen_info['Acquisition nr'])
            par.slice = int(gen_info['Max. number of slices/locations'])
            par.dyn = int(gen_info['Max. number of dynamics'])
            _skip_comment_lines(parfile)
            rows = []
            line = None
            first_row = None
#1   1    1  1 0 2     0  16    50  256  256     0.00000   7.11331 8.67964e-03  1070  1860   0.00   0.00   0.00  -20.00   20.00   60.00 10.000 10.000 0 2 0 2  0.977  0.977   4.60    0.00     0.00    0.00   1   15.00     0    0    0    64   0.0  1   1    7    0   0.000    0.000    0.000  1
            #col_types = (int, int, int, int, int, int, int, int, int, int, int,
            #    float, float, float, int, int, float, float, float, float,
            #    float, float, float, float, int, int, int, int, float,
            #    float, float, float, float, float, int, float, int, int,
            #    int, int, float, int, int, str, str, float, float, float,
            #    int)
            while line != "":
                line = parfile.readline()
                if line.strip() == "":  # Empty line
                    continue
                if line.startswith("#"):  # End of image info
                    break
                cols = line.strip().split()
                #Grab first 11 cols, 43rd for diffusion gradient number,
                #and then 42nd for inversion delay. We use [::-1] to reverse
                #the slice containing indices element 42 and 43
                row = [int(i) for i in cols[:11] + cols[41:43][::-1]]
                #row = [col_types[i](x) for i, x in enumerate(cols, 0)]
                if first_row is None:
                    first_row = cols
                rows.append(row)
            last_row = cols
            """
            #All columns
            rectypes = [('slice_number', int), ('echo_number', int),
                ('dynamic_scan_number', int), ('cardiac_phase_number', int),
                ('image_type_mr', int), ('scanning_sequence', int),
                ('index_in_REC_file', int), ('image_pixel_size', int),
                ('scan_percentage', int), ('recon_resolution_x', int),
                ('recon_resolution_y', int), ('rescale_intercept', float),
                ('rescale_slope', float), ('scale_slope', float),
                ('window_center', int), ('window_width', int),
                ('image_angulation_ap', float), ('image_angulation_fh', float),
                ('image_angulation_rl', float), ('image_offcentre_ap', float),
                ('image_offcentre_fh', float), ('image_offcentre_rl', float),
                ('slice_thickness', float), ('slice_gap', float),
                ('image_display_orientation', int), ('slice_orientation', int),
                ('fmri_status_indication', int), ('image_type_ed_es', int),
                ('pixel_spacing_x', float), ('pixel_spacing_y', float),
                ('echo_time', float), ('dyn_scan_begin_time', float),
                ('trigger_time', float), ('diffusion_b_factor', float),
                ('number_of_averages', int), ('image_flip_angle', float),
                ('cardiac_frequency', int), ('minimum_rr_interval', int),
                ('maximum_rr_interval', int), ('turbo_factor', int),
                ('inversion_delay', float), ('diffusion_b_value_number', int),
                ('gradient_orientation_number', int),
                ('contrast_type', str), ('diffusion_anisotropy_type', str),
                ('diffusion_ap', float), ('diffusion_fh', float),
                ('diffusion_rl', float), ('label_type', int)]
            """
            slice_index = np.array(rows)
            #If there is more than one slice, check the order of the
            #slice numbers. 1 = ascending order, 2 = descending order
            if slice_index.shape[0] > 1:
                par.slices_sorted = (2, 1)[slice_index[1,0] > slice_index[0,0]]
            else:
                par.slices_sorted = 2
            #from V 4.2 and higher, par file slice index lines contain
            #more info on diffusion gradients. For data versions below V4.2,
            #add zero columns.
            if par.ResToolsVersion in ('V4', 'V4.1'):
                nrows = slice_index.shape[0]
                zeros = np.zeros((nrows, 2), dtype=np.int)
                slice_index = np.concatenate((slice_index, zeros), axis=1)
            if par.dyn > 1:
                #estimate scan-duration from dtime PAR file row
                par.RT = (float(last_row[31]) - float(first_row[31])) / (par.dyn
                    - 1)
            else:
                par.RT = np.nan
            r_s, r_i, s_s, issue = iterate_loadpar(parfilename)
            par.slice_index = slice_index[:,0:13]
            #rectypes = [('slice_number', int), ('echo_number', int),
            #    ('dynamic_scan_number', int), ('cardiac_phase_number', int),
            #    ('image_type_mr', int), ('scanning_sequence', int),
            #    ('index_in_REC_file', int), ('image_pixel_size', int),
            #    ('scan_percentage', int), ('recon_resolution_x', int),
            #    ('recon_resolution_y', int), ('diffusion_b_value_number', int)]
            #par.slice_index_as_rec = np.rec.fromrecords(rows, dtype=rectypes)
            if len(first_row) < 15:
                logger.error("Problem: No actual slices measured for volume "
                    "{0}. Volume might be corrupt.".format(parfilename))
                par.problemreading = 1
            else:
                par.sliceorient = int(first_row[25])
                x = int(first_row[9])
                y = int(first_row[10])
                z = par.slice
                par.dim = np.array([x, y, z])
                if issue:
                    par.rescale_slope = r_s
                    par.rescale_interc = r_i
                    par.issue = True
                    par.scale_slope = s_s
                else:
                    par.rescale_slope = 1 / float(first_row[13])
                    par.rescale_interc = float(first_row[11])
                    par.issue = False
                par.bit = int(first_row[7])
                par.slth = float(first_row[22])
                par.gap = float(first_row[23])
                voxx = float(first_row[28])
                voxy = float(first_row[29])
                voxz = par.slth + par.gap
                par.vox = np.array([voxx, voxy, voxz])
                fovz, fovx, fovy = (float(f) for f in
                    gen_info['FOV (ap,fh,rl) [mm]'].split())
                par.fov = np.array([fovx, fovy, fovz])
                par.fov_apfhrl = np.array([fovz, fovx, fovy])
                par.angAP, par.angFH, par.angRL = (float(f) for f in
                    gen_info['Angulation midslice(ap,fh,rl)[degr]'].split())
                par.offAP, par.offFH, par.offRL = (float(f) for f in
                    gen_info['Off Centre midslice(ap,fh,rl) [mm]'].split())
    logger.debug("PARInfo {0}".format(par))
    return par

def iterate_loadpar(parfilename):
    """
        Parameters:
            parfilename: fullpath to the par file
        Returns:
           rescaleSlope: an array of the rescaleSlope values
           rescaleIntercept: an array of the rescaleIntercept values
           scaleSlope: an array of the scalSlope values
           issue: boolean of true if any of the arrays do not have 1 unique
            element
    """
    logger = logging.getLogger('raw2nii')
    scaleSlope = []
    rescaleSlope = []
    issue = False
    rescaleIntercept = []
    with open(parfilename, 'rb') as parfile:
        line = None
        while line != "":
            line = parfile.readline()
            if line.strip() == "":  # Empty line
                continue
            if line[0] != '.' and line[0] != '#':
                test = line.split()[:20]
                scaleSlope.append(float(test[13]))
                rescaleSlope.append(float(test[12]))
                rescaleIntercept.append(float(test[11]))
    issue = (np.product(np.unique(scaleSlope).shape) != 1
        or np.product(np.unique(rescaleIntercept).shape) != 1
        or np.product(np.unique(rescaleSlope).shape) != 1)
    if issue:
        logger.warning('Multiple scaling factors detected. Switching to float '
            '32 nifti and rescaling')
    return rescaleSlope, rescaleIntercept, scaleSlope, issue

if __name__ == "__main__":
    logger = logging.getLogger('raw2nii')
    logger.setLevel(logging.DEBUG)
    read_par("/Users/matthew/Downloads/abc.DICOM.06-03-2015_11.30.53AM/XMLPARREC/Morgan_211536_01_01_15.27.16_(WIP_Survey_SHC32).PAR")
