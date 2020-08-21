import numpy as np
from pathlib import Path
import struct
import matplotlib.pyplot as plt
from pprint import pprint


def load_mitgcm_single_mask(filename):

    mask = np.zeros([16, 20])

    # Initializing mask:
    for col in range(5,15):
        mask[8][col] = col+1.0

    #mask.ravel(order='F').astype('>f4').tofile(filename)
    mask.ravel(order='C').astype('>f4').tofile(filename)
    #print("mask column ordered: ", np.fromfile(filename, dtype='>f4'))
    print("mask row ordered: ", np.fromfile(filename, dtype='>f4'))
    #print("mask row ordered: ", mask.ravel(order='C'))


def load_mitgcm_boxregion_masks(dir, bathy_fname, xs, xe, ys, ye, dimx, dimy):
    '''
    Creates a open boundary masks to surround a box region within the dimx by dimy global array

    params:
    xs :: x index start wrt dimx by dimy global array
    xe :: x index end wrt dimx by dimy global array
    ys :: y index start wrt dimx by dimy global array
    ye :: y index end wrt dimx by dimy global array
    dimx :: x dimension of global array
    dim :: y dimension of global array
    '''

    mask1 = np.zeros([dimy, dimx])
    mask2 = np.zeros([dimy, dimx])
    mask3 = np.zeros([dimy, dimx])
    mask4 = np.zeros([dimy, dimx])
    mask_combined = np.zeros([dimy, dimx])

    mask1_count = 1
    mask2_count = 1
    mask3_count = 1
    mask4_count = 1
    mask_comb_count = 1

    # Creating mask:
    for col in range(xs, xe+1):
        mask1[ys][col] = mask1_count
        mask_combined[ys][col] = mask1_count

        mask2[ye][col] = mask2_count
        mask_combined[ye][col] = mask2_count
        mask1_count += 1
        mask2_count += 1

    for row in range(ys, ye+1):
        mask3[row][xs] = mask3_count
        mask_combined[row][xs] = mask3_count

        mask4[row][xe] = mask4_count
        mask_combined[row][xe] = mask4_count

        mask3_count += 1
        mask4_count += 1

    # Creating mask that combines all above
    '''
    for col in range(xs, xe):
        mask_combined[ys][col] = mask_comb_count
        mask_comb_count += 1

    for row in range(ys, ye):
        mask_combined[row][xe] = mask_comb_count
        mask_comb_count += 1

    for col in range(xe, xs, -1):
        mask_combined[ye][col] = mask_comb_count
        mask_comb_count += 1

    for row in range(ye, ys, -1):
        mask_combined[row][xs] = mask_comb_count
        mask_comb_count += 1
    '''


    fname1 = dir + "flt32_mask1.bin"
    fname2 = dir + "flt32_mask2.bin"
    fname3 = dir + "flt32_mask3.bin"
    fname4 = dir + "flt32_mask4.bin"
#fname_comb = dir + "flt32_mask_comb.bin"

    #Writing to file in big endian format for float32 values:
    mask1.ravel(order='C').astype('>f4').tofile(fname1)
    mask2.ravel(order='C').astype('>f4').tofile(fname2)
    mask3.ravel(order='C').astype('>f4').tofile(fname3)
    mask4.ravel(order='C').astype('>f4').tofile(fname4)
#    mask_combined.ravel(order='C').astype('>f4').tofile(fname_comb)

    #mask.ravel(order='F').astype('>f4').tofile(filename)


    #print("mask column ordered: ", np.fromfile(filename, dtype='>f4'))
    print("mask1 row ordered printing row ",ys,":", np.fromfile(fname1, dtype='>f4').reshape(dimy,dimx)[ys])
    print("mask2 row ordered printing row ",ye,":", np.fromfile(fname2, dtype='>f4').reshape(dimy,dimx)[ye])
    print("mask3 row ordered printing col ",xs,":", np.fromfile(fname3, dtype='>f4').reshape(dimy,dimx)[:,xs])
    print("mask4 row ordered printing col ",xe,":", np.fromfile(fname4, dtype='>f4').reshape(dimy,dimx)[:,xe])
    #print("mask row ordered: ", mask.ravel(order='C'))

    # Graphing mask combined
    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(mask_combined, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title("Region: Mask combined")
    plt.show()

    bathy = np.fromfile(bathy_fname, dtype='>f4').reshape(dimy,dimx)
    plt.figure(num=2, figsize=(7,6))
    plt.imshow(bathy, origin='lower')
    plt.show()

    land_mask = np.where(bathy ==0, -1, 0)
    plt.figure(num=3, figsize=(7,6))
    plt.imshow(land_mask, origin='lower')
    plt.show()

    land_mask_ob = land_mask * 0
    land_mask_ob = land_mask_ob + mask1
    land_mask_ob = land_mask_ob + mask2
    land_mask_ob = land_mask_ob + mask3
    land_mask_ob = land_mask_ob + mask4

    plt.figure(num=4, figsize=(7,6))
    plt.imshow(land_mask_ob, origin='lower')
    plt.show()

    land_mask_ob = land_mask * 0
    land_mask_ob = np.where(mask1 > 0, land_mask_ob + 1, land_mask_ob)
    land_mask_ob = np.where(mask2 > 0, land_mask_ob + 2, land_mask_ob)
    land_mask_ob = np.where(mask3 > 0, land_mask_ob + 3, land_mask_ob)
    land_mask_ob = np.where(mask4 > 0, land_mask_ob + 4, land_mask_ob)

    plt.figure(num=5, figsize=(7,6))
    plt.imshow(land_mask_ob, origin='lower')
    plt.show()



if __name__ == "__main__":
    dir_ob_mitgcm = "/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/input_ob/"
    bathy_fname = "bathy.labsea1979"
    '''
    Will create a 30 by 10 box region in the pacific ocean
    Within the 90 by 40 array:
    x: 40 - 70
    y: 15 - 25
    '''
    load_mitgcm_boxregion_masks(dir_ob_mitgcm, bathy_fname, 5, 15, 4, 12, 20, 16)
