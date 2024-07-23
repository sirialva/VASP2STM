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

    def cut_plane_idx(self):
        # Changes indices of x,y (plane) and z (out-of-plane) depending on desired cut direction.
        if self.cut_direction == 'c':
            x = 0
            y = 1
            z = 2
        elif self.cut_direction == 'b':
            x = 0
            y = 2
            z = 1
        elif self.cut_direction == 'a':
            x = 1
            y = 2
            z = 0
        else:
            print("ERROR: Cut direction is not specified. Default is cut perpendicular to c-axis.")
            x = 0
            y = 1
            z = 2
        return x,y,z

    def make_cut(self,tip_height,
                 self.scan_mode='constant_current',
                 self.cut_direction='c',
                 supercell=[1,1]):
        '''
        cut_direction_options = ['a', 'b', 'c']
        scan_mode_options=['constant_height','constant_current','2Dslice']

        NOT IMPLEMENTED BELOW:
        Select the STM scan mode to obtain 2D slice for visualization:
        1: Constant-height mode;[Default option],
        2: Constant-current mode;,
		3: 2D-slice at a specified height.

        Select the way of setting the height of tip position to make the 2D slice. 
        [1] Specifiying one height. [Default]
        [2] Specifiying a range of heights. Makes multiple images.
         
        '''

        self.heights = []
        self.heights.append(tip_height)

        x,y,z = cut_plane_idx(self.cut_direction)

        # The size of grids along x/y (in-plane) vectors.
        supercell_xngridpoints = (ngridpoints[x]-1)*supercell[x]+1
        supercell_yngridpoints = (ngridpoints[y]-1)*supercell[y]+1

        # Make arrays of x and y values.
        supercell_xarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        supercell_yarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)

        # Make arrays of supercell_density2D/I/H with the same dimensions as x/y arrays.
        supercell_density2D = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        I = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        H = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)

        #Find projection of y vector onto a vector.
        ytox = np.dot(unit_cell[x],unit_cell[y].T)/cell_lengths[x]
        #Find component of y vector perpendicular to a vector.
        ynormal = np.cross(unit_cell[x],unit_cell[y].T)/cell_lengths[x]
        ynormal = np.sqrt(np.dot(ynormal,ynormal.T))

        ##################################

        for h in heights:
            # Contant height mode.
            if scan_mode_option == "1":
                n1 = h/cell_lengths[2]*ngridpoints[2]
                dn1 = n1-np.floor(n1)
                n1 = int(n1)%ngridpoints[2]
                ldos = (1-dn1)*density[:,:,n1]+dn1*density[:,:,(n1+1)%ngridpoints[2]]

            # Contant current mode.
            elif self.scan_mode == 'constant_current':
                n2 = h/self.cell_lengths[z]*self.ngridpoints[z]
                dn2 = n2-np.floor(n2)
                n2 = int(n2)%self.ngridpoints[z]
                # Get the averaged current.
                averaged_current = ((1-dn2)*self.density[:,:,n2].mean()+dn2*self.density[:,:,(n2+1)%self.ngridpoints[2]].mean())
                c1 = self.density[:,:,n2]
                c2 = self.density[:,:,(n2+1)%self.ngridpoints[z]]
            
            # 2D-slice at a specified height.
            else:
                plane_index = int(round(h/cell_lengths[2]*ngridpoints[2]))%ngridpoints[2]
                density2D = density[:,:,plane_index]
                
            for i in range(supercell_xngridpoints):
                for j in range(supercell_yngridpoints):
                    supercell_xarray[i][j] = float(i)/float(supercell_xngridpoints)*cell_lengths[0]*m+float(j)/float(supercell_yngridpoints)*ytox*n
                    supercell_yarray[i][j] = float(j)/float(supercell_yngridpoints)*ynormal*n
                    mi = i%(ngridpoints[0]-1)
                    nj = j%(ngridpoints[1]-1)

                    if scan_mode_option == "1":
                        I[i][j] = ldos[mi][nj]

                    elif scan_mode_option == "2":
                        if c2[mi][nj]-c1[mi][nj] == 0:
                            H[i][j] = n2*cell_lengths[2]/ngridpoints[2]

                        else:
                            H[i][j] = (n2+(averaged_current-c1[mi][nj])/(c2[mi][nj]-c1[mi][nj]))*cell_lengths[2]/ngridpoints[2]
                        
                    else:
                        supercell_density2D[i][j] = density2D[mi][nj]

    def plot_cut():