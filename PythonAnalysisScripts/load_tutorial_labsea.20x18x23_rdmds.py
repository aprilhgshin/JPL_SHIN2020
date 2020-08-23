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

    files = np.sort(list( output_dir.glob(basename + '*meta')))

    for file in files:
        print(file)
        tmp = (file.stem).split('.')
        
        time_step = int(tmp[1])
#        tile_x = int(tmp[2])
#        tile_y = int(tmp[3])
        
#        print(time_step, tile_y, tile_x)
        if time_step not in time_steps:
            time_steps.append(time_step)
        
            
    print('time steps')
    pprint(time_steps)
    
    
    return time_steps
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
    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/run_ob/')

    output_filetype = '>f'


    # find the time steps and the list of i's and j's of the tiles
    
    time_steps =  \
        find_time_steps_and_tiles(output_dir, basename = 'diagsSI')
    
    print(time_steps)
    
    # DOMAIN DIMENSIONS nx, ny, nz    
    nz = 23 # depths
    ny = 16 # rows    
    nx = 20 # cols

    print('domain size')
    domain_size = [nz, ny, nx]
    pprint(domain_size)
    
    bathy_fname = output_dir / 'bathy.labsea1979'
    bathy = np.fromfile(bathy_fname, dtype='>f4').reshape(ny,nx)
    land_mask = np.where(bathy == 0, np.nan, 1)
    plt.figure(1000);
    plt.subplot(121);plt.imshow(bathy,origin='lower')
    plt.subplot(122);plt.imshow(land_mask, origin='lower')

    #%%    
    # DETERMINE TILE SIZE based on DOMAIN SIZE and tile_yS, tile_xS
    
    tile_size = [23, 10, 8]
    print('tile size')
    pprint(tile_size)

    ## LOAD A SINGLE INSTANCE OF THE TEMPERATURE FIELD

    # mitgcm python util 'rdmds'
    
    #def rdmds(fnamearg,itrs=-1,machineformat='b',rec=None,fill_value=0,
    #      returnmeta=False,astype=float,region=None,lev=(),
    #      usememmap=False,mm=False,squeeze=True,verbose=False):
        

    basename = 'diagsSI'
     # load the field, the result is 3D

    diagsSI,diagsSI_its, diagsSI_meta = mitgcm.rdmds(str(output_dir / basename), \
                                    returnmeta=True,\
                                    itrs=np.NaN, fill_value=-9999)
    
    print(diagsSI.shape)
    
    plt.close(1)
    plt.figure(num=1,clear=True, figsize=(20,6))
    
    for i in range(5):
        plt.subplot(2,3,i+1)
        plt.imshow(diagsSI[0,i,:], origin='lower')#, vmin=-2, vmax=16)
        plt.colorbar()
        plt.title(diagsSI_meta['fldlist'][i])
                  
    #%%
    
    basename = 'diagsTSUVW'
     # load the field, the result is 3D

    diagsTSUVW, diagsTSUVW_its, diagsTSUVW_meta = mitgcm.rdmds(str(output_dir / basename), \
                                    returnmeta=True,\
                                    itrs=np.NaN, fill_value=-9999)
    
    print(diagsTSUVW.shape)
    
    plt.close(2)
    plt.figure(num=2,clear=True, figsize=(20,6))
    
    for i in range(5):
        plt.subplot(2,3,i+1)
        plt.imshow(diagsTSUVW[0,i,0,:], origin='lower')#, vmin=-2, vmax=16)
        plt.colorbar()
        plt.title(diagsTSUVW_meta['fldlist'][i])       


    #%%

    mask_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/lab_sea/run_ob/')
    
    mask_files = np.sort(list(mask_dir.glob('flt32*mask*.bin')))
    
    masks = []
    masks_length = []
    for mf in mask_files:
        print(mf)
        masks.append(np.fromfile(mf, dtype='>f4').reshape(ny,nx).astype(np.int))
        masks_length.append(np.max(masks[-1]))

    plt.close(3)
    plt.figure(num=3,figsize=(12,12))
    
    plt.close(4)
    plt.figure(num=4,figsize=(12,12))
    
    for i in range(len(mask_files)):
 
        tmp = np.where(masks[i] > 0, 1, 0)
        tmp = np.where(land_mask == 1, tmp + .5, tmp)
       
        plt.figure(3)
        plt.subplot(2,2,i+1);
        plt.imshow(tmp,origin='lower')
        plt.title('MASK ' + str(i) + ' : length = ' + str(masks_length[i]))
    
        plt.figure(4)
        plt.subplot(2,2,i+1);
        plt.imshow(masks[i],origin='lower')
        plt.title('MASK ' + str(i) + ' : length = ' + str(masks_length[i]))
    
    #%%
    
    nt = diagsSI.shape[0]
    nf = diagsSI.shape[1]
    nm = len(mask_files)
    
    SI_obs = []
    for m in range(nm):
        SI_obs.append(np.zeros((nt, nf, masks_length[m])))
        
    
    for m in range(nm):
        for t in range(nt):
            for f in range(nf):
                tmp_field = diagsSI[t,f,:]
                for ob_point in range(masks_length[m]):
                    SI_obs[m][t,f,ob_point] = \
                            tmp_field[np.where(masks[m]==ob_point+1)]
        
    
    plt.close(4)
    plt.figure(num=4,figsize=(15,8))
    for m in range(nm):
        for f in range(nf):
            plt.subplot(nm,nf,m*nf + f+1);
            plt.imshow(SI_obs[m][:,f,:]);plt.axis('auto')
            plt.colorbar()
            plt.title(f'{diagsSI_meta["fldlist"][f]} : MASK {m}')
   
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
            
            
    #%%
    TSUVW_obs = []
    for i in range(nm):
        TSUVW_obs.append(np.zeros((nt, nf, nz, masks_length[i])))
        
    
    for m in range(nm): # loop throught masks
        for t in range(nt): # loop through time
            for f in range(nf): # loop through fields
                for k in range(nz): # loop through vertical levels
                    tmp_field = diagsTSUVW[t,f,k,:]
                    # loop through points on the mask [1...n]
                    for ob_point in range(masks_length[m]):
                        TSUVW_obs[m][t,f,k,ob_point] = \
                            tmp_field[np.where(masks[m]==ob_point+1)]
            
    
    k1 = 1
    k2= 5
    
    plt.close(5)
    plt.figure(num=5,figsize=(15,8))
    for m in range(nm):
        for f in range(nf):
            plt.subplot(nm,nf,m*nf + f+1);
            plt.imshow(TSUVW_obs[m][:,f,k1,:]);plt.axis('auto')
            plt.colorbar()
            plt.title(f'{diagsTSUVW_meta["fldlist"][f]} : k {k1} : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
       
    plt.close(6)
    plt.figure(num=6,figsize=(15,8))
    for m in range(nm):
        for f in range(nf):
            plt.subplot(nm,nf,m*nf + f+1);
            plt.imshow(TSUVW_obs[m][:,f,k2,:]);plt.axis('auto')
            plt.colorbar()
            plt.title(f'{diagsTSUVW_meta["fldlist"][f]} : k {k2} : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
    
             
    #%% LOAD MASK_NN_FIELD.bin files
 
    fields_2D = ['AREA','ETAN','HEFF','HSNOW','UICE','VICE']
    fields_3D = ['SALT','THETA','UVEL','VVEL','WVEL']
   
    mask_files = []          
    ob = dict()
    obf = dict()
    
    field_names = []
    for m in range(nm):
        print('\nMASK ' + str(m))
        mask_length = masks_length[m]
        print('-- length: ' + str(mask_length))

        tmp_str = 'MASK_' + str.zfill(str(m+1),2) + '*'
        print(tmp_str)
        print('\n')
        mask_files = np.sort(list(output_dir.glob(tmp_str)))
        
        tmp = None
        
        for file in mask_files:
            print('\nloading ', file.name)
            field_name = file.name[8:].split('.')[0]
            tmp = np.fromfile(file, dtype='>f4')
            
            print(len(tmp))
            if field_name not in field_names:
                field_names.append(field_name)
                ob[field_name] = dict()
               
              
            if field_name in fields_2D:
               
                print(f'{field_name} is 2D')
                nt = int(len(tmp)/mask_length)
                tmp = np.reshape(tmp, [nt, mask_length])
               
            elif field_name in fields_3D:
                print(f'{field_name} is 3D')
                nt = int(len(tmp)/mask_length/nz)
                tmp = np.reshape(tmp, [nt, nz, mask_length])

            else:
                print(f'{field_name} not found')
                
            ob[field_name][m] = tmp 
            obf[file.stem] = tmp
    #%%
            
    # plot all of the AREA open boundary files
    plt.figure(10,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['AREA'][m]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'AREA : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
   
    # plot all of the AREA open boundary files
    plt.figure(11,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['HEFF'][m]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'HEFF : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
   
     # plot all of the AREA open boundary files
    plt.figure(12,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['HSNOW'][m]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'HSNOW : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
   
    #%%
        
    # PLOT T S U V AT DIFFERENT LEVELS
   
    # pick z level
    k1=1
    
    # pick z level
    k2=5
    
    # plot all of the THETA open boundary files @ k1
    plt.figure(20,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['THETA'][m][:,k1,:]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'THETA : k {k1} : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
    
    # plot all of the THETA open boundary files @ k2
    plt.figure(21,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['THETA'][m][:,k2,:]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'THETA : k {k2} :  MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)   


    #SALT
    plt.figure(30,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['SALT'][m][:,k1,:]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'SALT : k {k1} : MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)       
  
    # pick z level
    k2=5
    
    # plot all of the THETA open boundary files
    plt.figure(31,clear=True, figsize=(3,10));
    for m in range(nm):
        plt.subplot(4,1,m+1);
        plt.imshow(ob['SALT'][m][:,k2,:]);
        plt.axis('auto');plt.colorbar()
        plt.title(f'SALT : k {k2} :  MASK {m}')
    plt.subplots_adjust(wspace=.3, hspace=.3,top=.95, bottom=.05)     
  