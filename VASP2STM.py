'''

VASP2STM.py

Based on STM-2DScan.py by ShuangLeung (sleung1924@gmail.com)             

'''

import sys, os
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ase.calculators.vasp import VaspChargeDensity

class VASP2STM():

    def _init_(self, filename='CHG'):
        vasp_charge = VaspChargeDensity(filename)
        self.density = vasp_charge.chg[-1]
        self.atoms = vasp_charge.atoms[-1]
        del vasp_charge
        
        # Read size of XYZ grids.
        self.ngridpoints = np.array(self.density.shape)

        # Read scaling factor and unit cell of the crystal structure.
        self.unit_cell = self.atoms.get_cell()
        self.cell_lengths = np.sqrt(np.dot(self.unit_cell,self.unit_cell.transpose()).diagonal())

    def make_cut(self,self.scan_mode='constant_current',tip_height):
        '''
        Select the STM scan mode to obtain 2D slice for visualization:
        1: Constant-height mode;[Default option],
        2: Constant-current mode;,
		3: 2D-slice at a specified height.

        Select the way of setting the height of tip position to make the 2D slice. 
        [1] Specifiying one height. [Default]
        [2] Specifiying a range of heights. Makes multiple images.
         
        '''
        heights = []

        scan_mode_options=['constant_height','constant_current','2Dslice']

        

    def plot_cut():