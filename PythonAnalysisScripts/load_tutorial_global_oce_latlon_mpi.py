# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint


#%%
# point this routine at a directory of model output and it will
# determine the time steps and the i's and j's of the tiles
# you have to give it the name of one of the output fields, 
# the default of 'Eta' since 'Eta' is likely output.
def find_time_steps_and_tile_ijs(output_dir, basename = 'Eta'):
    
    time_steps = []
    tile_is = []  # numbers of the ROW blocks (mitgcm y)
    tile_js = []  # numbers of the COL blocks (mitgcm x)

    files = np.sort(list( output_dir.glob(basename + '.000*meta')))
    
    for file in files:
        tmp = (file.stem).split('.')
        
        time_step = int(tmp[1])
        tile_j = int(tmp[2])
        tile_i = int(tmp[3])
        
        print(time_step, tile_i, tile_j)
        if time_step not in time_steps:
            time_steps.append(time_step)
        if tile_i not in tile_is:
            tile_is.append(tile_i)
        if tile_j not in tile_js:
            tile_js.append(tile_j)
            
    print('time steps')
    pprint(time_steps)
    
    print('tile is')
    pprint(tile_is)
    
    print('tile js')
    pprint(tile_js)
    
    return time_steps, tile_is, tile_js
#%%

# loads a single 2D or 3D field given the directory of the output
# the base field name, the time step, the list of tiles in i and j,
# the domain size, the tile size, the dimensions of the field (2 or 3)
# and the filetype of the output (default float32)
    
# it returns a 2D or 3D numpy array of the same dimensions as the 
# 2D or 3D domain_size, depending whether the field is 2D or 3D
    
def load_single_field(output_dir, basename, time_step,\
                      tile_is, tile_js, domain_size, tile_size, dims,\
                          filetype='>f'):
        
    if dims == 2:
        field = np.zeros((domain_size[1:]))
    elif dims == 3:
        field = np.zeros((domain_size))
        
    for tile_i in tile_is:
        for tile_j in tile_js:
            
            start_i = tile_size[1]*(tile_i -1)
            end_i = tile_size[1]*(tile_i )
            
            start_j = tile_size[2]*(tile_j-1)
            end_j = tile_size[2]*(tile_j)
            
            print(tile_i, tile_j)
            print(start_i, end_i, start_j, end_j)
                
            fname = create_filename_from_time_tiles(basename, time_step,\
                                                    tile_i, tile_j)
            print(fname)
            
            f = open(output_dir / fname, 'rb')
            dt = np.dtype(filetype)
    
            fileContent = np.fromfile(f, dtype=dt)
            if dims == 2:
                fileContent = np.reshape(fileContent, tile_size[1:])
                field[start_i:end_i, start_j:end_j] = fileContent
                
            elif dims==3:
                fileContent = np.reshape(fileContent, tile_size)
                field[:, start_i:end_i, start_j:end_j] = fileContent

    return field
#%%      


# makes a filename based on a field basename, the time step, and the 
# tile i and j value
def create_filename_from_time_tiles(basename, time_step, tile_i, tile_j):
    fname = basename + '.' + str(time_step).zfill(10) + \
        '.' + str(tile_j).zfill(3) + '.' + str(tile_i).zfill(3) + \
            '.data'
    return fname


#%%

if __name__ == "__main__":
       
    
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/tutorial_global_oce_latlon/run_mpi/output_1yr')
    output_filetype = '>f'

    # find the time steps and the list of i's and j's of the tiles
    
    time_steps, tile_is, tile_js =  \
        find_time_steps_and_tile_ijs(output_dir, basename = 'Eta')
    
    
    # DOMAIN DIMENSIONS nx, ny, nz    
    nk = 15 # depths
    nj = 90 # cols
    ni = 40 # rows
    
    print('domain size')
    domain_size = [nk, ni, nj]
    pprint(domain_size)
    
    # DETERMINE TILE SIZE based on DOMAIN SIZE and TILE_IS, TILE_JS
    tile_size = [nk, int(ni/len(tile_is)), int(nj/len(tile_js))]
    
    print('tile size')
    pprint(tile_size)


    ## LOAD A SINGLE INSTANCE OF THE TEMPERATURE FIELD
    basename = 'T'
    
    # temperature is 3D
    field_dims = 3
    
    # load the firsttime step 
    time_step_to_load = time_steps[0]

    # load the field, the result is 3D
    T_field = load_single_field(output_dir, basename, time_step_to_load,\
                              tile_is, tile_js, domain_size, tile_size, field_dims)

    
    print(T_field.shape)
    
    plt.close(1)
    plt.figure(num=1,clear=True, figsize=(7,2.85))
    plt.imshow(T_field[5,:], origin='lower', vmin=-2, vmax=16)
    plt.colorbar()
    plt.title('T @ k=0 and t= ' + str(time_step_to_load))


    ## LOAD ALL TIME STEPS OF THIS FIELD
    
    nt = len(time_steps)
    
    # create empty array of dimension time, depth, rows, cols
    T_all = np.zeros(([nt] + domain_size));
    
    # loop through all times, load file, add to T_all
    for ti, t in enumerate(time_steps):
        print(ti,t)
        T_all[ti,:] = load_single_field(output_dir, basename, t,\
                          tile_is, tile_js, domain_size, tile_size, field_dims)
        


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
                        