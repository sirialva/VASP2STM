'''

VASP2STM.py

Based on STM-2DScan.py by ShuangLeung (sleung1924@gmail.com) 

Written by Siri, July 2024.

'''

import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
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

    def swap_columns(self, array, column1, column2):
        if array.ndim == 1:
            temp = array[column1].copy()
            array[column1]  = array[column2]
            array[column2] =  temp
        elif array.ndim == 3:
            array = np.swapaxes(array,column1,column2)
        return array

    def make_cut(self, 
                 tip_height,
                 scan_mode='2D_slice',
                 cut_direction='c'):
        self.tip_height = tip_height

        # Choose cut_direction, rearranges coordinate order
        if cut_direction == 'c':
            pass
        elif cut_direction == 'b':
            self.unit_cell = self.swap_columns(self.unit_cell, 1, 2)
            self.cell_lengths = self.swap_columns(self.cell_lengths, 1, 2)
            self.density = self.swap_columns(self.density, 1, 2)
            print('NOTICE: Choosing cut_direction = b swaps b and c axis in data.')

        # Find index of cut height slice
        cut_height_idx = int(round(tip_height/self.cell_lengths[2]*self.density.shape[2]))
        
        # Find cut density
        if scan_mode == '2D_slice':
            self.cut_density = self.density[:,:,cut_height_idx]
        elif scan_mode == 'constant_current':
            self.cut_density = self.density[:,:,cut_height_idx]

    def make_supercell(self,size=1):
        # Tile cut_density to a supercell of given size, 
        # either give singular integer (gives a square) or [rows,columns] for not square
        self.size = size
        if isinstance(self.size, int):
            self.size = [self.size,self.size]
        self.cut_density = np.tile(self.cut_density, reps=self.size)

    def custom_cmap(self, colorlist=['xkcd:dark violet', 'xkcd:deep pink', 'xkcd:baby pink', 'xkcd:pale']):
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colorlist)
        return cmap
    
    def plot(self,save_fig=False):
        c = plt.matshow(self.cut_density, cmap=self.custom_cmap())
        plt.colorbar(c)
        plt.xticks(())
        plt.yticks(())

        if save_fig == True:
            plt.savefig(self.filename+'_'+str(round(self.tip_height,3))+'.pdf', dpi=300, bbox_inches='tight')