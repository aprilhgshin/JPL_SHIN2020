# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/home/mitgcm/Work/MITgcm/utils/python/MITgcmutils/')
import MITgcmutils as mitgcm


#%%
# point this routine at a directory of model output and it will
# determine the time steps and the i's and j's of the tiles
# you have to give it the name of one of the output fields, 
# the default of 'Eta' since 'Eta' is likely output.
def find_time_steps_and_tiles(output_dir, basename = 'Eta'):
    
    time_steps = []
    tile_xs = []  # numbers of the COL blocks (mitgcm x)
    tile_ys = []  # numbers of the ROW blocks (mitgcm y)

    files = np.sort(list( output_dir.glob(basename + '.000*meta')))

    for file in files:
        print(file)
        tmp = (file.stem).split('.')
        
        time_step = int(tmp[1])
        tile_x = int(tmp[2])
        tile_y = int(tmp[3])
        
        print(time_step, tile_y, tile_x)
        if time_step not in time_steps:
            time_steps.append(time_step)
        if tile_y not in tile_ys:
            tile_ys.append(tile_y)
        if tile_x not in tile_xs:
            tile_xs.append(tile_x)
            
    print('time steps')
    pprint(time_steps)
    
    print('tile is')
    pprint(tile_ys)
    
    print('tile js')
    pprint(tile_xs)
    
    return time_steps, tile_xs, tile_ys
#%%

# loads a single 2D or 3D field given the directory of the output
# the base field name, the time step, the list of tiles in i and j,
# the domain size, the tile size, the dimensions of the field (2 or 3)
# and the filetype of the output (default float32)
    
# it returns a 2D or 3D numpy array of the same dimensions as the 
# 2D or 3D domain_size, depending whether the field is 2D or 3D
    
# basename should be something like
# Eta.000000001 (variable with time)
# or
# Depth (field with no time)
    
def load_single_field(output_dir, basename,\
                      tile_xs, tile_ys, domain_size, tile_size, dims,\
                      filetype='>f'):
        
    plt.figure(10, clear=True)
    
    if dims == 2:
        field = np.zeros((domain_size[1:]))
    elif dims == 3:
        field = np.zeros((domain_size))
        
    num_cols = int( domain_size[2]/tile_size[1])
    num_rows = int(domain_size[1]/tile_size[2])
    print(num_cols, num_rows)
    
    for tile_x in tile_xs:
       for tile_y in tile_ys:

            if len(tile_ys) == 1:            
                cur_row = int((tile_x-1)/num_cols)
                cur_col =np.mod(tile_x-1, num_cols)
            else:
                cur_row = tile_y-1
                cur_col = tile_x-1
            
            start_x = cur_col*tile_size[1]
            end_x   = (cur_col+1)*tile_size[1]
            
            start_y = cur_row*tile_size[2]
            end_y   = (cur_row+1)*tile_size[2]
            
            print(start_x, end_x, start_y, end_y)
            print(cur_row, cur_col)
            fname = create_filename(basename, tile_x, tile_y)
            
            f = open(output_dir / fname, 'rb')
            dt = np.dtype(filetype)
            print(fname)
            fileContent = np.fromfile(f, dtype=dt)
            
            
            if dims == 2:
                fileContent = np.reshape(fileContent, tile_size[1:])
                field[start_y:end_y, start_x:end_x] = fileContent
                plt.subplot(4,9,tile_x)
                plt.imshow(fileContent);plt.colorbar()
            elif dims==3:
                fileContent = np.reshape(fileContent, tile_size)
                field[:, start_y:end_y, start_x:end_x] = fileContent
                
                plt.subplot(4,9,tile_x)
                plt.imshow(fileContent[0,:])
        
    return field, fileContent
#%%      


# makes a filename based on a field basename, and the 
# tile i and j value
def create_filename(basename, tile_x, tile_y):
    fname = basename + \
        '.' + str(tile_x).zfill(3) + '.' + str(tile_y).zfill(3) + \
            '.data'
    return fname


#%%

if __name__ == "__main__":
       
    
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/tutorial_global_oce_latlon/run_mpi/output_1yr')
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/global_ocean.90x40x15/run_mpi/')
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/global_ocean.90x40x15/run_mpi_nPx_36/')

    output_filetype = '>f'

    # find the time steps and the list of i's and j's of the tiles
    
    time_steps, tile_xs, tile_ys =  \
        find_time_steps_and_tiles(output_dir, basename = 'Eta')
    
    
    # DOMAIN DIMENSIONS nx, ny, nz    
    nz = 15 # depths
    ny = 40 # rows    
    nx = 90 # cols

    print('domain size')
    domain_size = [nz, ny, nx]
    pprint(domain_size)
    
    # DETERMINE TILE SIZE based on DOMAIN SIZE and tile_yS, tile_xS
    
    tile_size = [15, 10, 10]
    print('tile size')
    pprint(tile_size)

    ## LOAD A SINGLE INSTANCE OF THE TEMPERATURE FIELD
    basename = 'Eta.' + str(time_steps[0]).zfill(10)
    #basename = 'Depth'
    # temperature is 3D
    field_dims = 2
    
    # load the firsttime step 
    time_step_to_load = time_steps[0]

    # load the field, the result is 3D
    field,fc = load_single_field(output_dir, basename,\
                                tile_xs, tile_ys, domain_size, tile_size, field_dims)

    #def rdmds(fnamearg,itrs=-1,machineformat='b',rec=None,fill_value=0,
    #      returnmeta=False,astype=float,region=None,lev=(),
    #      usememmap=False,mm=False,squeeze=True,verbose=False):
    #%% 
        
    field = []
    basename = 'Depth'
    field, its, meta = mitgcm.rdmds(str(output_dir / basename), returnmeta=True)
    
    print(field.shape)
    
    plt.close(1)
    plt.figure(num=1,clear=True, figsize=(7,2.85))
    plt.imshow(field, origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    #plt.title('T @ k=0 and t= ' + str(time_step_to_load))


    #%%
    ## LOAD ALL TIME STEPS OF THIS FIELD
    
    nt = len(time_steps)
    
    # create empty array of dimension time, depth, rows, cols
    T_all = np.zeros(([nt] + domain_size));
    
    # loop through all times, load file, add to T_all
    for ti, t in enumerate(time_steps):
        print(ti,t)
        basename = 'Eta.' + str(t).zfill(10)

        T_all[ti,:],fc = load_single_field(output_dir, basename,\
                          tile_xs, tile_ys, domain_size, tile_size, field_dims)
        


    ## plot the last time step at 3 different depth levels
    plt.close(2)
    plt.figure(num=2,clear=True, figsize=(8,10))
    
    plt.clf()
    plt.subplot(311)
    plt.imshow(T_all[ti,0,:], vmin=-2, vmax=35, origin='lower')
    plt.colorbar()
    plt.title('T @ k=0 and t= ' + str(t))
    
    plt.subplot(312)
    plt.imshow(T_all[ti,4,:], vmin=-2, vmax=15, origin='lower')
    plt.colorbar()
    plt.title('T @ k=4 and t= ' + str(t))
    
    plt.subplot(313)
    plt.imshow(T_all[ti,8,:], vmin=-2, vmax=5, origin='lower')
    plt.colorbar()
    plt.title('T @ k=8 and t= ' + str(t))
   
    plt.show()
                        