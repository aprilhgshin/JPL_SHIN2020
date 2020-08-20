import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/home/mitgcm/Work/MITgcm/utils/python/MITgcmutils/')
import MITgcmutils as mitgcm

def test_ob_outputs3D(output_dir, mask_dir, ob_mask, ob_output, fname, depth):

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)
#    field = np.fromfile(str(output_dir / field), dtype='>f4')

    fname1 = fname+".0000036000.001.001"
    fname2 = fname + ".0000036000.002.001"

    field1 = mitgcm.rdmds(str(output_dir / fname1), itrs=-1)#.reshape(40,90)
    field2 = mitgcm.rdmds(str(output_dir / fname2), itrs=-1)#.reshape(40,90)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field1[0], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.subplot(212)
    plt.imshow(field2[0], origin='lower')
    plt.colorbar()
    plt.show()

    plt.title("mask")
    plt.imshow(mask, origin='lower')
    plt.colorbar()
    plt.show()

    print(mask.shape)

    ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
    depth = int(len(ob_out)/90)
    output = ob_out.reshape(depth,90)

    newFieldOnMask = [[0 for x in range(90)] for y in range(depth)]
    diff = [[0 for x in range(90)] for y in range(depth)]

    print(len(field1))
    print(field1[0][2][46])
    print(len(mask[0])-1)
    print("output shape:",output.shape)



    counter = 0

    for row in range(len(mask)-1):
        for col in range(len(mask[0])-1):
#                print("mask >= 1:", mask[row][col], "field:", field1[0][row][col])
#                print("corresponding output:", output[counter])
            if(mask[row][col] > 0):
                for k in range(depth-1):
                    if (col <= 45):
                        newFieldOnMask[k][counter] =  field1[k][row][col]
                    elif (col > 45):
                        newFieldOnMask[k][counter] =  field2[k][row][col]
                counter += 1

    diff = (output - newFieldOnMask)
    print("difference between field and output:",diff>0.1)

    for k in range(len(field1)):
        plt.subplot(211)
        plt.plot(output[k], 'r--', linewidth=2.5, label="ob output")
        plt.plot(newFieldOnMask[k], 'b', label="original field")
        plt.title("OB Output VS. Original Field")
        plt.legend()

        plt.subplot(212)
        plt.plot(diff[k])
        plt.title("Difference between Output and Original field:")
        plt.show()



def test_ob_outputs3D_L1(output_dir, mask_dir, ob_mask, ob_output, fname):

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)
#    field = np.fromfile(str(output_dir / field), dtype='>f4')

    fname1 = fname+".0000036000.001.001"
    fname2 = fname + ".0000036000.002.001"

    field1 = mitgcm.rdmds(str(output_dir / fname1), itrs=-1)#.reshape(40,90)
    field2 = mitgcm.rdmds(str(output_dir / fname2), itrs=-1)#.reshape(40,90)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field1[0], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.subplot(212)
    plt.imshow(field2[0], origin='lower')
    plt.colorbar()
    plt.show()

    plt.title("mask")
    plt.imshow(mask, origin='lower')
    plt.colorbar()
    plt.show()

    print(mask.shape)

    output = np.fromfile(str(output_dir / ob_output), dtype='>f8').reshape(15,90)[0]

    newFieldOnMask = np.zeros(len(output))
    diff = np.zeros(len(output))

    print(field1[0][2][46])
    print(len(mask[0])-1)

    counter = 0

    for row in range(len(mask)-1):
        for col in range(len(mask[0])-1):
#                print("mask >= 1:", mask[row][col], "field:", field1[0][row][col])
#                print("corresponding output:", output[counter])
            if(mask[row][col] > 0):
                if (col <= 45):
                    newFieldOnMask[counter] =  field1[0][row][col]
                elif (col > 45):
                    newFieldOnMask[counter] =  field2[0][row][col]
                counter += 1

    print("len:", len(newFieldOnMask))
    diff = (output - newFieldOnMask)
    print("difference between field and output:",diff>0.1)

    plt.subplot(211)
    plt.plot(output, 'r--', linewidth=2.5, label="ob output")
    plt.plot(newFieldOnMask, 'b', label="original field")
    plt.title("OB Output VS. Original Field")
    plt.legend()

    plt.subplot(212)
    plt.plot(diff)
    plt.title("Difference between Output and Original field:")
    plt.show()





def new_test_ob_outputs3D_L1(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec):

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)

    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=36001)#np.NaN)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field[0][0], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.show()

    if (filePrec == 32):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
        depth = int(len(ob_out)/90)
        output = ob_out.reshape(depth,90)
        depth = int(depth/2)
    elif (filePrec == 64):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
        depth = int(len(ob_out)/90)
        output = ob_out.reshape(depth,90)

    print(output.shape)
    print("depth:",depth)

    newFieldOnMask = [[0 for x in range(90)] for y in range(depth)]
    diff = [[0 for x in range(90)] for y in range(depth)]
    print("shape of newFieldOnMask:",len(newFieldOnMask[0]))

    counter = 0

    for row in range(len(mask)):
        for col in range(len(mask[0])):
            if(mask[row][col] > 0):
                print("counter:",counter,"field:",field[fieldNum][0][row][col])
                for k in range(depth):
                    if (field[fieldNum][k][row][col] == 0):
                        newFieldOnMask[k][counter] = None
                    else:
                        newFieldOnMask[k][counter] =  field[fieldNum][k][row][col]
                counter += 1

    for i in range(90):
        for k in range(depth):
            if(output[k][i] == 0):
                output[k][i] = None

    print("output shape",output.shape)
    print("output first element:",output[0][0])
    print("newFieldOnMask last element",newFieldOnMask[0][-1] )
    print("newFieldOnMask:",newFieldOnMask[0])

    for k in range(depth):
        for i in range(len(newFieldOnMask)):
            if (newFieldOnMask[k][i] == None):
                diff[k][i] = None
            else:
                diff[k][i] = abs(output[k][i] - newFieldOnMask[k][i])
        print("abs difference between field and output:",diff[k])

        plt.subplot(211)
        plt.plot(output[k], 'r--', linewidth=2.5, label="ob output")
        plt.plot(newFieldOnMask[k], 'b', label="original field")
        plt.title("OB Output VS. Original Field")
        plt.legend(loc="center right")

        plt.subplot(212)
        plt.plot(diff[k])
        plt.title("Abs Difference between Output and Original field:")
        plt.show()



def test_ob_outputs2D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum):

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)


    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=36001)#np.NaN)

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

    ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
    depth = int(len(ob_out)/90)
    output = ob_out.reshape(depth,90)
    print(output.shape)

    newFieldOnMask = np.zeros(len(output[0]))
    diff = np.zeros(len(output[0]))

    counter = 0

    for row in range(len(mask)-1):
        for col in range(len(mask[0])-1):
            if(mask[row][col] > 0):
                if (field[row][col] == 0):
                    newFieldOnMask[counter] = None
                else:
#                newFieldOnMask[counter] =  field[fieldNum][row][col]
                    newFieldOnMask[counter] =  field[row][col]

                counter += 1

    for i in range(90):
            if(output[i] == 0):
                output[i] = None

    print("len:", len(newFieldOnMask))
    diff = (output[0] - newFieldOnMask)
    print("difference between field and output:",diff>0.1)


    plt.subplot(211)
    plt.plot(output[0], 'r--', linewidth=2.5, label="ob output")
    plt.plot(newFieldOnMask, 'b', label="original field")
    plt.title("OB Output VS. Original Field")
    plt.legend()

    plt.subplot(212)
    plt.plot(diff)
    plt.title("Difference between Output and Original field:")
    plt.show()



if __name__ == "__main__":



    fld_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/diags')
    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run')
    mask_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/input')
#    field = []
#    arr = np.zeros(45)
#    newfield=[]


#    test_ob_outputs3D(output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_01_THETA   _00036001.bin", 'T', 1)
#    test_ob_outputs2D(fld_dir, output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_01_THETA   _00036001.bin", 'THETA')
###PARAMS: test_ob_outputs2D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname)
#    test_ob_outputs2D(fld_dir, output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_01_ETAN    _00036001.bin", 'ETAN', 0)

### PARAMS: new_test_ob_outputs3D_L1(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum
    new_test_ob_outputs3D_L1(fld_dir, output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_01_SALT    _00036001.bin", 'SALT', 1, 64)
#    new_test_ob_outputs3D_L1(fld_dir, output_dir, mask_dir, "domain_flt32_mask1.bin", "MASK_01_THETA   _00036001.bin", 'THETA', 0, 64)

#    field = mitgcm.rdmds(str(output_dir / 'THETA'), itrs=-1,fill_value=0.0)
#    print(np.fromfile(str(output_dir / 'THETA.001.001.data'), dtype='>f4').shape)
#    field2 = (np.fromfile(str(output_dir / 'THETA.data'), dtype='>f4'))
#    print(field2.shape)
#    newfield = (np.append(arr,field2)).reshape(40,90)
#    print(newfield)
#    lines = open(str(output_dir / 'THETA   _on_mask1_global.    122400.meta'))
#    print(lines)
    #print(field.shape)

#    for line in lines:
#        print(line)




#    plt.close(1)
#    plt.figure(num=1,clear=True, figsize=(7,6))
#    plt.subplot(211)
#    plt.imshow(field, origin='lower')#, vmin=-2, vmax=16)
#    plt.colorbar()
#    plt.title('ETAH on Open Boundary Points')

#    plt.show()
#
