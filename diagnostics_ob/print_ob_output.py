import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/home/mitgcm/Work/MITgcm/utils/python/MITgcmutils/')
import MITgcmutils as mitgcm


def test_ob_outputs3D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter, num_obPnts):
    '''
    Some params:
    ob_mask :: filename of file containing open boundary mask
    ob_output :: filename of file containing open boundary output
    fname :: field name
    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
    filePrec :: file precision
    myIter :: iter number at which file was outputted
    '''
    
    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)
    full_field = np.zeros([40,90])

    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=myIter)#np.NaN)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field[0][0], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.show()

    if (filePrec == 32):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
        depth = int(len(ob_out)/num_obPnts)
        output = ob_out.reshape(depth,num_obPnts)
        depth = int(depth/2)
    elif (filePrec == 64):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
        depth = int(len(ob_out)/num_obPnts)
        output = ob_out.reshape(depth,num_obPnts)

    print(output.shape)
    print("depth:",depth)

    newFieldOnMask = [[0 for x in range(num_obPnts)] for y in range(depth)]
    diff = [[0 for x in range(num_obPnts)] for y in range(depth)]
    print("shape of newFieldOnMask:",len(newFieldOnMask[0]))

    counter = 0

    for k in range(depth):
        print("K",k)
    for row in range(len(mask)):
        for col in range(len(mask[0])):
            if(mask[row][col] > 0):
                print("counter:",counter,"field:",field[fieldNum][0][row][col])
                for k in range(depth):
                    if (field[fieldNum][k][row][col] == 0):
                        newFieldOnMask[k][counter] = None
                    else:
                        newFieldOnMask[k][counter] =  field[fieldNum][k][row][col]
                        full_field[row][col] = field[fieldNum][k][row][col]

                counter += 1

    for i in range(num_obPnts):
        for k in range(depth):
            if(output[k][i] == 0):
                output[k][i] = None


    print("output shape",output.shape)
    print("output:",output)
    print("newFieldOnMask last element",newFieldOnMask[0][-1] )
    print("newFieldOnMask:",newFieldOnMask[0])

    for k in range(depth):
        for i in range(len(newFieldOnMask[0])):
            if (newFieldOnMask[k][i] == None):
                diff[k][i] = None
            else:
                print(output[k][i], newFieldOnMask[k][i])
                diff[k][i] = abs(output[k][i] - newFieldOnMask[k][i])
        print("abs difference between field and output:",diff[k])

        plt.figure(num=1,clear=True, figsize=(7,6))
        plt.subplot(211)
        plt.imshow(full_field, origin='lower')#, vmin=-2, vmax=16)
        plt.colorbar()
        plt.title("field values on OB points")
        plt.show()

        plt.subplot(211)
        plt.plot(output[k], 'r--', linewidth=2.5, label="ob output")
        plt.plot(newFieldOnMask[k], 'b', label="original field")
        plt.title("OB Output VS. Original Field")
        plt.legend(loc="center right")

        plt.subplot(212)
        plt.plot(diff[k])
        plt.title("Abs Difference between Output and Original field:")
        plt.show()



def test_ob_outputs2D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter):
# NOTE: the array field may or may not have a third dimension if there are multiple fields outputted from the diagnostics package.
# Change accordingly:
# If yes:  field[fieldNum][y][x]
# Otherwise: field[y]][x]

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)
    full_field = np.zeros([40,90])


    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=myIter)#np.NaN)

    print(field.shape)
    print(itrs)
    print(meta)


    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field, origin='lower')#, vmin=-2, vmax=16)
#    plt.imshow(field[fieldNum], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.show()

    if (filePrec == 32):
        output = np.fromfile(str(output_dir / ob_output), dtype='>f4')
    elif (filePrec == 64):
        output = np.fromfile(str(output_dir / ob_output), dtype='>f8')


    print(output.shape)
    newFieldOnMask = np.zeros(len(output))
    diff = np.zeros(len(output))

    counter = 0

    for row in range(len(mask)):
        for col in range(len(mask[0])):
            if(mask[row][col] > 0):
                if (field[row][col] == 0):
                    newFieldOnMask[counter] = None
                else:
#                newFieldOnMask[counter] =  field[fieldNum][row][col]
                    newFieldOnMask[counter] =  field[row][col]
                    full_field[row][col] = field[row][col]

                counter += 1

    for i in range(len(output)):
            if(output[i] == 0):
                output[i] = None

    print("len:", len(newFieldOnMask))
    diff = abs(output - newFieldOnMask)
    print("abs difference between field and output:",diff)
    print("output",output)
    print("newFieldOnMask",newFieldOnMask)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(full_field, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title("field values on OB points")
    plt.show()

    plt.subplot(211)
    plt.plot(output, 'r--', linewidth=2.5, label="ob output")
    plt.plot(newFieldOnMask, 'b', label="original field")
    plt.title("OB Output VS. Original Field")
    plt.legend()

    plt.subplot(212)
    plt.plot(diff)
    plt.title("Abs Difference between Output and Original field:")
    plt.show()



if __name__ == "__main__":



    fld_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/diags')
    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run')
    mask_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/input')

#    Plotting and comparing 2D field outputs:
###PARAMS: fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter
#    test_ob_outputs2D(fld_dir, output_dir, mask_dir, "flt32_mask4.bin", "MASK_04_ETAN    _00036007.bin", 'ETAN', 0, 64, 36007)

#    Plotting and comparing 3D field outputs:
### PARAMS: fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter, num_obPnts
#    ob_mask :: filename of file containing open boundary mask
#    ob_output :: filename of file containing open boundary output
#    fname :: field name
#    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
#    filePrec :: file precision
#    myIter :: iter number at which file was outputted

#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_04_SALT    _00036007.bin", 'SALT', 1, 64, 36007, 10)
    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask2.bin", "MASK_02_THETA   _00036007.bin", 'THETA', 0, 64, 36007, 20)
#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask2.bin", "MASK_02_SALT    _00036007.bin", 'SALT', 1, 64, 36007, 30)
#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask3.bin", "MASK_03_SALT    _00036007.bin", 'SALT', 1, 64, 36007, 10)
#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask4.bin", "MASK_04_SALT    _00036007.bin", 'SALT', 1, 64, 36007, 10)
