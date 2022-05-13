import os
import sys
import numpy as np
from utils_ifum import pack_4fits

def main(argv, arc):
    file_num = argv[1] #"0164"
    work_dir = "/Users/yysong/git2/ifum_aperMap_example"
    dir_raw = os.path.join(work_dir,"data_raw/ut20220512")
    dir_new = os.path.join(work_dir,"data_packed")
    img_mask_flag = False
    img_mask_path = 'None'
    img_mask_name = 'None'

    for Shoe in ['b','r']:
        file_name_temp = Shoe+file_num
        pack_4fits(file_name_temp,dir_raw,dir_new,img_mask_flag,img_mask_path,img_mask_name)

if __name__ == '__main__':
    main(sys.argv, len(sys.argv))
