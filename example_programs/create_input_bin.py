'''
The "domain" will be 40x10,
A 40x10 array of ints will define the 'mask',  all zeros except for 1s through the 5th row, columns 5:25
A 40x10 array of floats defining the 'temperature' will be read in and distributed to the processors,
'''

import numpy as np

def load_domain_mask(filename):

    mask = np.zeros([10, 40])

    # initializing mask:
    for col in range(5, 26):
        mask[4][col] = 1

    mask.ravel(order='F').astype('int32').tofile(filename)
    print("mask column ordered: ", np.fromfile(filename, dtype=np.int32))
    print("mask row ordered: ", mask.ravel(order='C'))

def load_domain_temp(filename):
    temp  = np.zeros([10,40])
    row_len, col_len = temp.shape

    for row in range(len(temp)):
        for col in range(col_len):
            temp[row][col] = (row+1)*(col+1)
    temp.ravel(order='F').astype('float32').tofile(filename)
    print("temp column ordered: ", np.fromfile(filename, dtype=np.float32))
    print("temp row ordered: ", temp.ravel(order='C'))

def main():
    filename_mask = "./input_domains/domain_mask.bin"
    filename_temp = "./input_domains/domain_temp.bin"

    load_domain_mask(filename_mask)
    load_domain_temp(filename_temp)

if __name__ == "__main__":
    main()
