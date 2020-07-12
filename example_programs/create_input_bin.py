'''
The "domain" will be 40x10,
A 40x10 array of ints will define the 'mask',  all zeros except for 1s through the 5th row, columns 5:25
A 40x10 array of floats defining the 'temperature' will be read in and distributed to the processors,
'''

import numpy as np
from pathlib import Path
import struct

def load_domain_mask(filename):

    mask = np.zeros([10, 40])

    # initializing mask:
    for col in range(5, 26):
        mask[4][col] = 1

    mask.ravel(order='F').astype('int32').tofile(filename)
    print("mask column ordered: ", np.fromfile(filename, dtype=np.int32))
    print("mask row ordered: ", mask.ravel(order='C'))

def load_mitgcm_domain_mask(filename):

    mask = np.zeros([40, 90])

    # initializing mask:
    for col in range(90):
        mask[12][col] = col+1.0

    #file = open(filename, 'wb')
    #file.write(mask)
    #file.close()

    #dt = np.dtype(('>f4'))
    #mask.ravel(order='F').astype('>f4').tofile(filename)
    mask.ravel(order='C').astype('>f4').tofile(filename)
    #print("mask column ordered: ", np.fromfile(filename, dtype='>f4'))
    print("mask row ordered: ", np.fromfile(filename, dtype='>f4'))
    #print("mask row ordered: ", mask.ravel(order='C'))

def load_domain_temp(filename):
    temp  = np.zeros([10,40])
    row_len, col_len = temp.shape

    for row in range(len(temp)):
        for col in range(col_len):
            temp[row][col] = (row+1)*(col+1)
    temp.ravel(order='F').astype('float32').tofile(filename)
    print("temp column ordered: ", np.fromfile(filename, dtype=np.float32))
    print("temp row ordered: ", temp.ravel(order='C'))

def load_mitgcm_domain_temp(filename):
    temp  = np.zeros([40,90])
    row_len, col_len = temp.shape

    for row in range(len(temp)):
        for col in range(col_len):
            temp[row][col] = (row+1)*(col+1)


    #temp.ravel(order='F').astype('>f4').tofile(filename)
    temp.ravel(order='C').astype('>f4').tofile(filename)
    print("temp row ordered: ", np.fromfile(filename, dtype='>f4'))

    #print("temp column ordered: ", np.fromfile(filename, dtype='>f4'))
    #print("temp row ordered: ", temp.ravel(order='C'))

def create_two_masks(filename1, filename2):
    mask1 = np.zeros([10, 40])
    mask2 = np.zeros([10, 40])

    for col in range(20):
        mask1[0][col] = col+1
        mask2[0][col] = 20 - col

    mask1.ravel(order='F').astype('int32').tofile(filename1)
    mask2.ravel(order='F').astype('int32').tofile(filename2)
    print("mask1 column ordered: ", np.fromfile(filename1, dtype=np.int32))
    print("mask1 row ordered: ", mask1.ravel(order='C'))
    print("mask2 column ordered: ", np.fromfile(filename2, dtype=np.int32))
    print("mask2 row ordered: ", mask2.ravel(order='C'))

def read_mitgcm_bin():
    filename = "/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/bathymetry.bin"
    with open(filename, 'rb') as f:
        data = np.fromfile(f, dtype=np.float32, count=-1)
        print(data)
    #print("bathy column ordered: ", np.fromfile(filename, dtype=np.float32, count=-1)


def main():
    filename_mask = "./input_domains/domain_mask.bin"
    filename_temp = "./input_domains/domain_temp.bin"
    filename1 = "./input_domains/domain_mask1.bin"
    filename2 = "./input_domains/domain_mask2.bin"

    fname_ob_mitgcm = "/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/domain_flt32_mask.bin"
    fname_temp_mitgcm = "/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/domain_temp.bin"

    #load_domain_mask(filename_mask)
    #load_domain_temp(filename_temp)
    #create_two_masks(filename1, filename2)

    load_mitgcm_domain_temp(fname_temp_mitgcm)
    #read_mitgcm_bin()
    load_mitgcm_domain_mask(fname_ob_mitgcm)

if __name__ == "__main__":
    main()
