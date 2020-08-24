import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/home/mitgcm/Work/MITgcm/utils/python/MITgcmutils/')
import MITgcmutils as mitgcm


def test_ob_outputs3D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter, depth, nTimeLevels, tLevel):
    '''
    Some params:
    ob_mask :: filename of file containing open boundary mask
    ob_output :: filename of file containing open boundary output
    fname :: field name
    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
    filePrec :: file precision
    myIter :: iter number at which file was outputted
    depth :: Depth of field
    nTimeLevels :: number of time levels included in output binary file
    tLevel :: time level to compare
    '''

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(16,20)
    full_field = np.zeros([16,20])

    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "diagsTSUVW"), returnmeta=True, itrs=myIter)#np.NaN)
    print("field shape",field.shape)

#    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=np.NaN, fill_value=-9999)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field[0][0], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)
    plt.show()


    if (filePrec == 32):
        if (nTimeLevels == 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
            print("ob_out",ob_out)
            print(ob_out)
            num_obPnts = int((len(ob_out)/depth))
            print(num_obPnts)
            output = ob_out.reshape(depth,num_obPnts)
        elif (nTimeLevels > 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
            num_obPnts = int(len(ob_out)/(depth*nTimeLevels))
            output = ob_out.reshape(nTimeLevels,depth,num_obPnts)

    elif (filePrec == 64):
        if (nTimeLevels == 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
            print(ob_out)
            print(ob_out.shape)

            num_obPnts = int(len(ob_out)/depth)
            output = ob_out.reshape(depth,num_obPnts)
        elif (nTimeLevels > 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
            print(ob_out.shape)

            num_obPnts = int(len(ob_out)/(depth*nTimeLevels))
            output = ob_out.reshape(nTimeLevels,depth,num_obPnts)

    print("output.shape",output.shape)
    print("depth:",depth)
    print("num_obPnts",num_obPnts)
    print(output)


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
          if (nTimeLevels == 1):
              if(output[k][i] == 0):
                  output[k][i] = None
          elif (nTimeLevels > 1):
              if(output[tLevel][k][i] == 0):
                  output[tLevel][k][i] = None

    print("output shape",output.shape)
    print("output:",output)
    print("newFieldOnMask:",newFieldOnMask)

    for k in range(depth):
        for i in range(len(newFieldOnMask[0])):
            if (newFieldOnMask[k][i] == None):
                diff[k][i] = None
            else:
                if (nTimeLevels == 1):
                    print(output[k][i], newFieldOnMask[k][i])
                    diff[k][i] = abs(output[k][i] - newFieldOnMask[k][i])
                elif (nTimeLevels > 1):
                    print(output[tLevel][k][i], newFieldOnMask[k][i])
                    diff[k][i] = abs(output[tLevel][k][i] - newFieldOnMask[k][i])
        print("abs difference between field and output:",diff[k])

        plt.figure(num=1,clear=True, figsize=(7,6))
        plt.subplot(211)
        plt.imshow(full_field, origin='lower')#, vmin=-2, vmax=16)
        plt.colorbar()
        plt.title("field values on OB points")
        plt.show()

        plt.subplot(211)
        if (nTimeLevels == 1):
            plt.plot(output[k], 'r--', linewidth=2.5, label="ob output")
        elif (nTimeLevels > 1):
            plt.plot(output[tLevel][k], 'r--', linewidth=2.5, label="ob output")
        plt.plot(newFieldOnMask[k], 'b', label="original field")
        plt.title("OB Output VS. Original Field For Each Depth")
        plt.legend(loc="center right")

        plt.subplot(212)
        plt.plot(diff[k])
        plt.title("Abs Difference between Output and Original field:")
        plt.show()


def test_ob_outputs3D_allTime(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, depth, nx, ny):
    '''
    Some params:
    ob_mask :: filename of file containing open boundary mask
    ob_output :: filename of file containing open boundary output
    fname :: field name
    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
    filePrec :: file precision
    myIter :: iter number at which file was outputted
    depth :: Depth of field
    nTimeLevels :: number of time levels included in output binary file
    tLevel :: time level to compare
    '''

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(ny,nx)
    full_field = np.zeros([ny,nx])

#    mask_test = np.fromfile(str(output_dir / "m1.bin"), dtype='>f8').reshape(ny,nx)

    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "obDiag"), returnmeta=True, itrs=np.NaN, fill_value=-9999)
    print("field shape",field.shape)

    nTimeLevels = len(field)

    plt.figure(num=6,clear=True, figsize=(7,6))
    plt.imshow(mask, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title(fname)


    if (filePrec == 32):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
        print(ob_out.shape)
        num_obPnts = int(len(ob_out)/(depth*nTimeLevels))
        output = ob_out.reshape(nTimeLevels,depth,num_obPnts)

    elif (filePrec == 64):
        ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
        print(ob_out.shape)
        num_obPnts = int(len(ob_out)/(depth*nTimeLevels))
        output = ob_out.reshape(nTimeLevels,depth,num_obPnts)

    print("output.shape",output.shape)
    print("depth:",depth)
    print("num_obPnts",num_obPnts)


    newFieldOnMask = np.zeros([nTimeLevels, num_obPnts])
    diff = [[0.0 for x in range(num_obPnts)] for y in range(nTimeLevels)]

    counter = 0
    # Depth set to k = 1
    k = 1


    for row in range(len(mask)):
        for col in range(len(mask[0])):
            if(mask[row][col] > 0):
                for t in range(nTimeLevels):
                    if (field[t][fieldNum][k][row][col] == 0):
                        newFieldOnMask[t][counter] = None
                    else:
                        newFieldOnMask[t][counter] =  field[t][fieldNum][k][row][col]
                        full_field[row][col] = field[0][fieldNum][k][row][col]

                counter += 1

#    output_copy = [[0.0 for x in range(num_obPnts)] for y in range(nTimeLevels)]
    output_copy = np.zeros([nTimeLevels, num_obPnts])
    for i in range(num_obPnts):
        for t in range(nTimeLevels):
          if(output[t][k][i] == 0):
              output_copy[t][i] = None
          else:
              output_copy[t][i] = output[t][k][i]

    print("output_copy shape",output_copy.shape)
    print("output_cpy:",output_copy)
    print("newFieldOnMask",newFieldOnMask)

#    print("output:",output)
#    print("newFieldOnMask:",newFieldOnMask[0])
    for t in range(nTimeLevels):

        for i in range(len(newFieldOnMask[0])):
            if (newFieldOnMask[t][i] == None):
                diff[t][i] = None
            else:
                print(output_copy[t][i], newFieldOnMask[t][i])
                diff[t][i] = abs(output_copy[t][i] - newFieldOnMask[t][i])
#            print("abs difference between field and output:",diff)






    plt.figure(num=5,clear=True, figsize=(7,6))
    plt.imshow(full_field, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title("field values on OB points")
#    plt.show()

#    plt.figure(num=1,clear=True, figsize=(7,6))
#    plt.imshow(mask_test, origin='lower')#, vmin=-2, vmax=16)
#    plt.colorbar()
#    plt.title("mask_test")

    plt.figure(num=2,clear=True, figsize=(7,6))
    plt.imshow(output_copy)#, vmin=-2, vmax=16)
    plt.title(f'{fname} : k {k} : MASK 3 output')
    plt.colorbar()
#    plt.show()

    plt.figure(num=3,clear=True, figsize=(7,6))
    plt.imshow(newFieldOnMask)
    plt.title(f'{fname} : k {k} : MASK 3 actual')
    plt.colorbar()
#    plt.show()


    plt.figure(num=4,clear=True, figsize=(7,6))
    plt.imshow(diff, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title("Difference")
    plt.show()


def test_ob_outputs2D(fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter, nTimeLevels, tLevel):
    '''
    Some params:
    ob_mask :: filename of file containing open boundary mask
    ob_output :: filename of file containing open boundary output
    fname :: field name
    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
    filePrec :: file precision
    myIter :: iter number at which file was outputted
    nTimeLevels :: number of time levels included in output binary file
    tLevel :: time level to compare


    IMPORTANT NOTE: the array field may or may not have a third dimension if there are multiple fields outputted from the diagnostics package.
    Change accordingly:
    If yes:  field[fieldNum][y][x]
    Otherwise: field[y][x]
    ALL instances of field with extra fieldNum dimension are commented out for easy change: Line 172, 209, 210
    '''

    mask = np.fromfile(str(mask_dir / ob_mask), dtype='>f4').reshape(40,90)
    full_field = np.zeros([40,90])


    print(list(output_dir.glob("*")))
    field, itrs, meta = mitgcm.rdmds(str(fld_dir / "diagsTSUVW"), returnmeta=True, itrs=myIter)#np.NaN)

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
        if (nTimeLevels == 1):
            output = np.fromfile(str(output_dir / ob_output), dtype='>f4')
            num_obPnts = len(output)
        elif (nTimeLevels > 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f4')
            num_obPnts = int(len(ob_out)/(nTimeLevels))
            output = ob_out.reshape(nTimeLevels, num_obPnts)

    elif (filePrec == 64):
        if (nTimeLevels == 1):
            output = np.fromfile(str(output_dir / ob_output), dtype='>f8')
            num_obPnts = len(output)
        elif (nTimeLevels > 1):
            ob_out = np.fromfile(str(output_dir / ob_output), dtype='>f8')
            num_obPnts = int(len(ob_out)/(nTimeLevels))
            output = ob_out.reshape(nTimeLevels, num_obPnts)


    print(output.shape)
    newFieldOnMask = np.zeros(num_obPnts)
    diff = np.zeros(num_obPnts)

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

    for i in range(num_obPnts):
      if (nTimeLevels == 1):
          if(output[i] == 0):
              output[i] = None
      elif (nTimeLevels > 1):
          if(output[tLevel][i] == 0):
              output[tLevel][i] = None


    print("len:", len(newFieldOnMask))
    if (nTimeLevels == 1):
        diff = abs(output - newFieldOnMask)
    elif (nTimeLevels > 1):
        diff = abs(output[tLevel] - newFieldOnMask)

    print("abs difference between field and output:",diff)
    print("output",output[tLevel])
    print("newFieldOnMask",newFieldOnMask)

    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(full_field, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title("field values on OB points")
    plt.show()

    plt.subplot(211)
    if (nTimeLevels == 1):
        plt.plot(output, 'r--', linewidth=2.5, label="ob output")
    elif (nTimeLevels > 1):
        plt.plot(output[tLevel], 'r--', linewidth=2.5, label="ob output")
    plt.plot(newFieldOnMask, 'b', label="original field")
    plt.title("OB Output VS. Original Field")
    plt.legend()

    plt.subplot(212)
    plt.plot(diff)
    plt.title("Abs Difference between Output and Original field:")
    plt.show()



if __name__ == "__main__":

    import struct

    fld_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/diags')
    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run')
    mask_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/input')
#    fld_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/run_ob')
#    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/run_ob')
#    mask_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/input_ob')
    '''
    NOTE:
    1. When changing data.diagnostics_ob and rerunning model, make sure to empty the run directory first before rerunning.
    2. When changing data.diagnostics, make sure to empty diags directory in the run directory before rerunning the model.
    '''
#============================================================================================================================
    '''
    Plotting and comparing 2D field outputs:
    PARAMS: fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter,  nTimeLevels, tLevel
    Same descriptions as parameters for 3D
    '''

# Example for multiple time levels:
#    test_ob_outputs2D(fld_dir, output_dir, mask_dir, "flt32_mask3.bin", "MASK_03_THETA_00036030.bin", 'THETA', 0, 64, 36030,31, 29)
# Example for one time level: Always set nTimeLevels to 1
#    test_ob_outputs2D(fld_dir, output_dir, mask_dir, "flt32_mask1.bin", "MASK_01_ETAN_00036007.bin", 'ETAN', 0, 32, 36007,1, 0)

#============================================================================================================================

    '''
    # Plotting and comparing 3D field outputs:
    # PARAMS: fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, myIter, num_obPnts
    #    ob_mask :: filename of file containing open boundary mask
    #    ob_output :: filename of file containing open boundary output
    #    fname :: field name
    #    fieldNum :: Index to extract field from field array outputted from the diagnostics package (in data.diagnostics)
    #    filePrec :: file precision
    #    myIter :: iter number at which file was outputted
    #    depth :: Depth of field
    #    nTimeLevels :: number of time levels included in output binary file
    #    tLevel :: time level to compare (index from 0 through n time levels)
    '''
# Example for one time level: Always set nTimeLevels to 1
#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask3.bin", "MASK_03_THETA.bin", 'THETA', 0, 64, 36030, 15, 1, 0)
# Example for multiple time levels:
#    test_ob_outputs3D(fld_dir, output_dir, mask_dir, "flt32_mask1.bin", "MASK_01_THETA_00000002.bin", 'THETA', 0, 32, 2, 1, 1, 0)
#    test_ob_outputs3D_allTime(fld_dir, output_dir, mask_dir, "flt32_mask3.bin", "MASK_03_THETA.bin", 'THETA', 0, 32, 23, 20,16)

#fld_dir, output_dir, mask_dir, ob_mask, ob_output, fname, fieldNum, filePrec, depth, nx, ny
    test_ob_outputs3D_allTime(fld_dir, output_dir, mask_dir, "flt32_mask3.bin", "MASK_03_THETA.bin", 'THETA', 0, 32, 15, 90,40)
