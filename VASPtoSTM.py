'''

VASP2STM.py

Based on STM-2DScan.py by ShuangLeung (sleung1924@gmail.com) 

Written by Siri, July 2024.

'''

import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from ase.calculators.vasp import VaspChargeDensity

class charge_density():

    def __init__(self, filename='CHG'):
        self.filename = filename

        # Read CHG file
        vasp_charge = VaspChargeDensity(filename)
        self.density = vasp_charge.chg[-1]
        self.atoms = vasp_charge.atoms[-1]
        del vasp_charge
        
        # Read size of XYZ grids.
        self.ngridpoints = np.array(self.density.shape)

        # Read scaling factor and unit cell of the crystal structure.
        self.unit_cell = self.atoms.get_cell()
        self.cell_lengths = np.sqrt(np.dot(self.unit_cell,self.unit_cell.transpose()).diagonal())

    def make_cut(self, 
                 tip_height,
                 scan_mode='2D_slice',
                 cut_direction='c',
                 supercell=[1,1]):
        self.tip_height = tip_height

        # Choose cut_direction, rearranges coordinate order
        if cut_direction == 'c':
            pass
        elif cut_direction == 'b':
            pass

        # Find index of cut height slice
        cut_height_idx = int(round(tip_height/self.cell_lengths[2]*self.density.shape[2]))
        
        # Find cut density
        if scan_mode == '2D_slice':
            self.cut_density = self.density[:,:,cut_height_idx]
        elif scan_mode == 'constant_current':
            self.cut_density = self.density[:,:,cut_height_idx]
    
    def plot(self,save_fig=False):
        c = plt.matshow(self.cut_density, cmap='RdPu_r')
        plt.colorbar(c)
        plt.xticks(())
        plt.yticks(())

        if save_fig == True:
            plt.savefig(self.filename+'_'+str(round(self.tip_height,3))+'.pdf', dpi=300, bbox_inches='tight')