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

    field, its = mitgcm.rdmds(str(output_dir / 'THETA'), itrs=-1,rec=1)

    print(np.fromfile(str(output_dir / 'THETA.data'), dtype='>f4').shape)
#    lines = open(str(output_dir / 'THETA   _on_mask1_global.    122400.meta'))
#    print(lines)
    #print(field.shape)

#    for line in lines:
#        print(line)







#    plt.close(1)
#    plt.figure(num=1,clear=True, figsize=(7,6))
#    plt.subplot(211)
#    plt.imshow(field[0,:], origin='lower')#, vmin=-2, vmax=16)
#    plt.colorbar()
#    plt.title('THETA on OB @ k=0')

#    plt.show()
