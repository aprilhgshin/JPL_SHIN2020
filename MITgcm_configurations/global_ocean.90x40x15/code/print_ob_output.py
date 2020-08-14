import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pprint import pprint

import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.append('/home/mitgcm/Work/MITgcm/utils/python/MITgcmutils/')
import MITgcmutils as mitgcm

if __name__ == "__main__":



    output_dir = Path('/home/mitgcm/Work/JPL_SHIN2020/MITgcm_configurations/global_ocean.90x40x15/run/')
    field = []
    arr = np.zeros(45)
    newfield=[]

#    field = mitgcm.rdmds(str(output_dir / 'THETA'), itrs=-1,fill_value=0.0)
#    print(np.fromfile(str(output_dir / 'THETA.001.001.data'), dtype='>f4').shape)
    field2 = (np.fromfile(str(output_dir / 'THETA.data'), dtype='>f4'))
    print(field2.shape)
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
