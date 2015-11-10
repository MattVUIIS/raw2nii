from __future__ import division
import collections
import logging
import numpy as np
import os
import re

from default_par_vars import get_default_parvars

_ImgTag = collections.namedtuple('ImgTag', ('len', 'mem', 'var'))
IMAGE_INFORMATION_DEFINTION_LINE = ('# === IMAGE INFORMATION DEFINITION ======='
    '======================================')
IMAGE_INFORMATION_LINE = ('# === IMAGE INFORMATION ============================'
    '==============================')
OK_COMMENT = ('#  The rest of this file contains ONE line per image, this line '
    'contains the following information:')

def strtok(s, delim=' '):
    begin, middle, end = s.lstrip(delim).partition(delim)
    return begin, middle + end

def strmatch(s, arr):
    for idx, item in enumerate(arr):
        i = item.find(s)
        if i > -1:
            return idx
    return -1

def strfind(s, item):
    return [match.start() for match in re.finditer(item, s)]

_READ_STATE_INFO_DEF = 0
_READ_STATE_DEFAULT = 1

def load_par(filename):
    logger = logging.getLogger('raw2nii')
    if not os.path.exists(filename):
        filename += 'v2'
    if not os.path.exists(filename):
        raise IOError('{0} does not exist'.format(filename))
    par, vars, slicevars = get_default_parvars()
    imgTagFormat = collections.defaultdict(list)
    doneslices = False
    read_state = _READ_STATE_DEFAULT
    with open(filename, 'rb') as fp:
        for line in fp:
            line = line.strip()
            if not line:
                par.fmt.append(('comment', line))
                continue
            if read_state == _READ_STATE_DEFAULT:
                c = line[0]
                if c == '#':  # comment
                    par.fmt.append(('comment', line))
                    if line == IMAGE_INFORMATION_DEFINTION_LINE:
                        read_state = _READ_STATE_INFO_DEF
                        imgTagFormat.clear()
                elif c == '.':  # scan file variable
                    tag, key = strtok(line, ':')
                    idx = strmatch(tag, vars[:,0])
                    if idx > -1:
                        if int(vars[idx,3]):
                            par[vars[idx,1]][vars[idx,2]] = np.fromstring(
                                key[1:].strip(), sep=' ')
                        else:
                            par[vars[idx,1]][vars[idx,2]] = key[1:].strip()
                        par.fmt.append(('var', vars[idx,1], vars[idx,2]))
                    else:
                        logger.warning('Unknown parameter '
                            'declaration: {0}'.format(line))
                        par.fmt.append(('comment', line))
                else:  # parse as image slice information
                    par = parseScanImgLine(par, line, imgTagFormat)
                    if not doneslices:
                        doneslices = True
                        par.fmt.append(('SLICES',))
            elif read_state == _READ_STATE_INFO_DEF:
                par.fmt.append(('comment', line))
                if line == IMAGE_INFORMATION_LINE:
                    read_state = _READ_STATE_DEFAULT
                else:
                    if len(line) > 1:
                        idx = strmatch(line, slicevars[:,0])
                        if idx == -1:
                            if line != OK_COMMENT:
                                logger.warning('Cannot interpret slice '
                                    'variable: {0}'.format(line))
                        else:
                            _len, mem, var = slicevars[idx][1:]
                            imgTagFormat['len'].append(int(_len))
                            imgTagFormat['mem'].append(mem)
                            imgTagFormat['var'].append(var)
                            par.scnfmt.append((line, mem, var))
    dynamic_scan_num = []
    slice_num = []
    for image in par.img:
        dynamic_scan_num.extend(image['info']['dynamic_scan_num'])
        slice_num.extend(image['info']['slice_num'])
    par.NumberOfVolumes = len(np.unique(dynamic_scan_num))
    NoV = len(dynamic_scan_num)//len(np.unique(slice_num))
    if NoV != par.NumberOfVolumes:
        logger.warning('Dynamic Scan Number does not match number of slices. '
            'Assuming slices are ordered.')
        cnt = [0] * max(slice_num)
        for s in par.img:
            slice_index = int(s['info']['slice_num'][0]) - 1
            cnt[slice_index] += 1
            s['info']['dynamic_scan_num'] = [cnt[slice_index]]
        dynamic_scan_num = []
        for image in par.img:
            dynamic_scan_num.extend(image['info']['dynamic_scan_num'])
        par.NumberOfVolumes = len(np.unique(dynamic_scan_num))
    orient = par.img[0]['orient']['orient']
    scn = par['scn']
    fov = scn['fov']
    if orient == 'TRA':
        if not np.allclose(fov[0], fov[2]):
            logger.warning('AXIAL (TRA): par.scn.fov[0]~=par.scn.fov[2]. '
                'Setting to max')
            fov = fov.astype('float')
            fov[[0, 2]] = np.max(fov[[0, 2]])
        fov_div_slices = fov[1] / par['max']['num_slices']
    elif orient == 'COR':
        if not np.allclose(fov[1], fov[2]):
            logger.warning('AXIAL (COR): par.scn.fov[1]~=par.scn.fov[2]. '
                'Setting to max')
            fov = fov.astype('float')
            fov[[1, 2]] = np.max(fov[[1, 2]])
        fov_div_slices = fov[0] / par['max']['num_slices']
    elif orient == 'SAG':
        if not np.allclose(fov[1], fov[0]):
            logger.warning('AXIAL (SAG): par.scn.fov[1]~=par.scn.fov[0]. '
                'Setting to max')
            fov = fov.astype('float')
            fov[[0, 1]] = np.max(fov[[0, 1]])
        fov_div_slices = fov[2] / par['max']['num_slices']
    if not np.allclose(scn['slicethk'], fov_div_slices):
        logger.warning('Slice Thickness does not match fov/num_slices. '
            'ADJUSTING!!!')
        scn['slicethk'] = fov_div_slices
    if not np.allclose(scn['slicegap'], 0):
        logger.warning('Non-zero slice gap: adjusting slice thickness')
        #Read the header variables
        scn['slicethk'] += scn['slicegap']
        scn['slicegap'] = 0
    #Determine File Volume/Slice Order (inputVolumeSliceOrder):
    # a) volume - all slices are listed (in order) for each volume before
    #   the next volume
    # b) slices - the same slice is listed for all volumes before the next
    # slice of the first volume (volumes are ordered)
    # c) other - some other ordering (any ordering of volumes/slices is
    # supported in the PAR file format)
    # Procedure: Build a matrix with each row having: [VOLUME SLICE IDX]
    N = len(par.img)
    SRT = np.zeros((N, 3), dtype=int)
    for j, s in enumerate(par.img):
        SRT[j,] = [s['info']['dynamic_scan_num'][0],
            s['info']['slice_num'][0], j]
    SRT = SRT[np.lexsort((SRT[:,2], SRT[:,1], SRT[:,0]))]
    is_in_volume_order = np.all(SRT[:,2] == np.arange(N))
    if is_in_volume_order:
        par.inputVolumeSliceOrder = 'volume'
    else:
        SRT = SRT[np.lexsort((SRT[:,2], SRT[:,0], SRT[:,1]))]
        is_in_volume_order = np.all(SRT[:,2] == np.arange(N))
        if is_in_volume_order:
            par.inputVolumeSliceOrder = 'slice'
        else:
            par.inputVolumeSliceOrder = 'unknown'
            logger.warning('Slice ordering is not a predefined type.')
            logger.info('This toolbox is compatible with arbitrary '
                'ordering of slices in PAR/REC files.')
            logger.info('However, other toolboxes or REC readers may '
                'assume a specific ordering.')
    return par

# Read
#
# === PIXEL VALUES =============================================================
#  PV = pixel value in REC file     FP = floating point value
#  DV = displayed value on console  RS = rescale slope
#  RI = rescale intercept           SS = scale slope
#  DV = PV * RS + RI                FP = DV / (RS * SS)
#
def parseScanImgLine(par, line, imgTagFormat):
    logger = logging.getLogger('raw2nii')
    h = np.fromstring(line, sep=' ')
    if h.shape[0] != np.sum(imgTagFormat['len']):
        raise ValueError('Slice tag format does not match the number of '
            'entries')
    img = collections.defaultdict(dict)
    k = 0
    for j in range(len(imgTagFormat['mem'])):
        mem = imgTagFormat['mem'][j]
        var = imgTagFormat['var'][j]
        _len = imgTagFormat['len'][j]
        img[mem][var] = h[k:k + _len]
        k += _len
    img['orient']['orient'] = 'UNK'
    if img['orient']['slice_orientation'] == 1:
        img['orient']['orient'] = 'TRA'
    elif img['orient']['slice_orientation'] == 2:
        img['orient']['orient'] = 'SAG'
    elif img['orient']['slice_orientation'] == 3:
        img['orient']['orient'] = 'COR'
    if 'pix_bits' not in par['scn']:
        par['scn']['pix_bits'] = img['info']['pix_bits'][:]
    elif 'pix_bits' in img['info']:
        if not np.allclose(par['scn']['pix_bits'], img['info']['pix_bits']):
            logger.warning('REC file contains image slices with variable bits '
                'per pixel.')
    if 'recon_res' not in par['scn']:
        par['scn']['recon_res'] = img['info']['recon_res'][:].T
    elif 'recon_res' in img['info']:
        if not np.allclose(par['scn']['recon_res'], img['info']['recon_res']):
            logger.warning('REC file contains image slices with variable '
                'reconstruction sizes.')
    if 'slicethk' not in par['scn']:
        par['scn']['slicethk'] = img['info']['slicethk'][:]
    elif 'slicethk' in img['info']:
        if not np.allclose(par['scn']['slicethk'], img['info']['slicethk']):
            logger.warning('REC file contains image slices with variable '
                'slice thickness.')
    if 'slicegap' not in par['scn']:
        par['scn']['slicegap'] = img['info']['slicegap'][:]
    elif 'slicegap' in img['info']:
        if not np.allclose(par['scn']['slicegap'], img['info']['slicegap']):
            logger.warning('REC file contains image slices with variable '
                'slice gap.')
    par.img.append(img)
    m = re.search(r'V(\d\.\d)', par.fmt[7][1])
    if m:
        par.version = float(m.group(1))
    else:
        par.version = -1
    return par
