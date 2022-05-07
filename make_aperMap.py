#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
"""
This script runs make_aperMap
"""

import os
import numpy as np

def readString_symbol(fileName, column1, symbol): #column number starts from 0
    x0 = []
    for line in open(fileName, 'r'):
        #skip over empty lines or lines starting with spaces or spaces+#
        tempString=line.strip()
        if (tempString[0] == '#' or tempString[0]==''):
            continue

        line = line.split(symbol)
        if line[column1]=='-':
            line[column1]=float('nan')

        x = line[column1].strip()
        x0.append(x)
    return np.array(x0)

def str2bool(s):
    flag = False
    true_synonym = ['1','yes','Yes','YES','true','True','TRUE']

    if np.sum(true_synonym==s)>0:
        flag = True

    return flag

def make_aperMap_usage():
    """
    Print make_aperMap usage description.
    """

    import textwrap

    #spclist = ', '.join(available_spectrographs)
    #spcl = textwrap.wrap(spclist, width=70)
    descs = '##  '
    descs += '\x1B[1;37;42m' + 'make_aperMap : '
    descs += 'The Python code for creating IFUM AperMap file v{0:s}'.format('1.0') \
              + '\x1B[' + '0m' + '\n'
    descs += '##  '
    #descs += '\n##  Available spectrographs include:'
    #for ispcl in spcl:
    #    descs += '\n##   ' + ispcl
    return descs

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description=make_aperMap_usage(),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('input_file', type=str,
                        help='File contains all the input parameters')
    #parser.add_argument('IFU', type=str,
    #                    help='IFU type: LSB, STD, HR')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):
    import time
    from utils_ifum import pack_4fits, write_aperMap

    print('Hello Yingyi')

    homeDir = os.getcwd()
    current_time = time.localtime()
    file_date = '%02d%02d%02d'%(current_time.tm_year%100,current_time.tm_mon,current_time.tm_mday)

    #### input parameters
    input_file = args.input_file
    keys = readString_symbol(input_file, 0, '=')
    strs = readString_symbol(input_file, 1, '=')

    IFU_type  = strs[np.where(keys=='IFU type')][0] # args.IFU # 'LSB'
    file_num  = strs[np.where(keys=='File number')][0] # args.file_num #'9991'
    Config  = strs[np.where(keys=='Config')][0] # args.file_num #'9991'
    path_work = strs[np.where(keys=='Working path')][0] # homeDir+'/data_raw'
    dir_raw = os.path.join(path_work, strs[np.where(keys=='Raw data dir')][0]) # homeDir+'/data_raw'
    dir_new = os.path.join(path_work, strs[np.where(keys=='Packed data dir')][0]) # homeDir+'/data_packed'
    dir_aperMap = os.path.join(path_work, strs[np.where(keys=='AperMap dir')][0]) # homeDir+'/AperMap'

    #### have an image mask?
    img_mask_flag = str2bool(strs[np.where(keys=='File number')][0]) # False #True
    if img_mask_flag:
        img_mask_path = strs[np.where(keys=='Image mask path')][0] # homeDir+'/img_mask'
        img_mask_name = strs[np.where(keys=='Image mask name')][0] # 'Config1_m2fs'
    else:
        img_mask_path = []
        img_mask_name = []

    #### have badfibers?
    add_badfiber_flag_b = str2bool(strs[np.where(keys=='Badfibers b-side')][0]) # False #True
    if add_badfiber_flag_b:
        add_badfiber_spat_id_b = np.array([np.int32(temp) for temp in strs[np.where(keys=='Badfibers b-side positions')][0].split(",")]) # np.array([982, 1173, 1960])
    else:
        add_badfiber_spat_id_b = np.array([])

    add_badfiber_flag_r = str2bool(strs[np.where(keys=='Badfibers r-side')][0]) # False #True
    if add_badfiber_flag_r:
        add_badfiber_spat_id_r = np.array([np.int32(temp) for temp in strs[np.where(keys=='Badfibers r-side positions')][0].split(",")]) # np.array([1395, 1414, 1972, 1982])
    else:
        add_badfiber_spat_id_r = np.array([])

    #### IFU parameters
    if IFU_type=='LSB':
        N_xx, N_yy = 18, 20
        width_y = 17.3
    elif IFU_type=='STD':
        N_xx, N_yy = 23, 24
        width_y = 10.0
    elif IFU_type=='HR':
        N_xx, N_yy = 27, 32
        width_y =  5.0
    elif IFU_type=='M2FS':
        N_xx, N_yy = 16, 16
        width_y =  5.0

    for Shoe in ['b','r']:
        #if Shoe == 'r':
        #    return 0

        #### step 1
        file_name_temp = Shoe+file_num

        print('++++ Start Step 1: %s'%file_name_temp)
        pack_4fits(file_name_temp,dir_raw,dir_new,img_mask_flag,img_mask_path,img_mask_name)
        print('++++ Done Step 1: %s\n'%file_name_temp)

        #### step 2
        traceFile = dir_new+'/'+file_name_temp+'.fits'

        path_MasterEdges = dir_new+'/Masters/MasterEdges_%s.fits.gz'%(file_name_temp)
        path_MasterSlits = dir_new+'/Masters/MasterSlits_%s.fits.gz'%(file_name_temp)

        print('++++ Start Step 2')
        if os.path.exists(path_MasterSlits):
            print("Caution: MasterSlits file exists!!! Continue to Step 3.")
        else:
            os.system('pypeit_trace_edges -t '+traceFile+' -s magellan_m2fs_blue')

            #path_MasterEdges_default = dir_new+'/Masters/MasterEdges_A_1_01.fits.gz'
            path_MasterEdges_default = dir_new+'/Masters/MasterEdges_A_1_DET01.fits.gz'
            os.system('mv '+path_MasterEdges_default+' '+path_MasterEdges)
            path_MasterSlits_default = dir_new+'/Masters/MasterSlits_A_1_DET01.fits.gz'
            os.system('mv '+path_MasterSlits_default+' '+path_MasterSlits)
        print('++++ Done Step 2 \n')

        #### step 3
        file_name = dir_aperMap+'/ap%s_%s_%s'%(Shoe,IFU_type,Config)

        if Shoe=='b':
            add_badfiber_flag, add_badfiber_spat_id = add_badfiber_flag_b, add_badfiber_spat_id_b
        elif Shoe=='r':
            add_badfiber_flag, add_badfiber_spat_id = add_badfiber_flag_r, add_badfiber_spat_id_r

        print('++++ Start Step 3')
        write_aperMap(path_MasterSlits, IFU_type,  Shoe, file_name, file_date, N_xx, N_yy, img_mask_flag,img_mask_path,img_mask_name, add_badfiber_flag, add_badfiber_spat_id)
        print('++++ Done Step 3\n')

    return 0


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()

