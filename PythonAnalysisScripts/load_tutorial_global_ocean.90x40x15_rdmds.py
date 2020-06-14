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



# makes a filename based on a field basename, and the 
# tile i and j value
def create_filename(basename, tile_x, tile_y):
    fname = basename + \
        '.' + str(tile_x).zfill(3) + '.' + str(tile_y).zfill(3) + \
            '.data'
    return fname


#%%

if __name__ == "__main__":
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/global_ocean.90x40x15/run_mpi/')
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/global_ocean.90x40x15/run_mpi_nPx_36/')
    output_dir = Path('/home/mitgcm/Work/MITgcm/verification/global_ocean.90x40x15/run_mpi_nPx_35_blanklist/')

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

    # mitgcm python util 'rdmds'
    
    #def rdmds(fnamearg,itrs=-1,machineformat='b',rec=None,fill_value=0,
    #      returnmeta=False,astype=float,region=None,lev=(),
    #      usememmap=False,mm=False,squeeze=True,verbose=False):
        
    field = []
    basename = 'T'
     # load the field, the result is 3D

    field, its, meta = mitgcm.rdmds(str(output_dir / basename), \
                                    returnmeta=True,\
                                    itrs=time_steps[0], fill_value=-9999)
    
    print(field.shape)
    
    plt.close(1)
    plt.figure(num=1,clear=True, figsize=(7,6))
    plt.subplot(211)
    plt.imshow(field[0,:], origin='lower')#, vmin=-2, vmax=16)
    plt.colorbar()
    plt.title('T @ k=0 and t= ' + str(time_steps[0]))

    plt.subplot(212)
    plt.imshow(field[0,:], origin='lower', vmin=-2, vmax=16)
    plt.colorbar()
    plt.title('T @ k=0 and t= ' + str(time_steps[0]))


    #%%
    ## LOAD ALL TIME STEPS OF THIS FIELD
    # the flag 'itrs=np.NaN' tells the program to load all time steps
    
    T_all, its, meta = mitgcm.rdmds(str(output_dir / basename),\
                                    returnmeta=True, 
                                    itrs=np.NaN)
    
 
    ## plot the last time step at 3 different depth levels
    plt.close(2)
    plt.figure(num=2,clear=True, figsize=(8,10))
    
    ti = -1
    plt.clf()
    plt.subplot(311)
    plt.imshow(T_all[ti,0,:], vmin=-2, vmax=35, origin='lower')
    plt.colorbar()
    plt.title('T @ k=0 and t= ' + str(time_steps[ti]))
    
    plt.subplot(312)
    plt.imshow(T_all[ti,4,:], vmin=-2, vmax=15, origin='lower')
    plt.colorbar()
    plt.title('T @ k=4 ')
    
    plt.subplot(313)
    plt.imshow(T_all[ti,8,:], vmin=-2, vmax=5, origin='lower')
    plt.colorbar()
    plt.title('T @ k=8 ')
   
    plt.show()
                        