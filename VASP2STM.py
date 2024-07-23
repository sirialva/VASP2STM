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

class charge_density:

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

    def make_cut(self,tip_height,
                 scan_mode='constant_current',
                 cut_direction='c',
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
        self.scan_mode = scan_mode
        self.cut_direction = cut_direction

        self.heights = []
        self.heights.append(tip_height)

        m = supercell[0]
        n = supercell[1]

        # The size of grids along x/y (in-plane) vectors.
        supercell_xngridpoints = (self.ngridpoints[0]-1)*supercell[0]+1
        supercell_yngridpoints = (self.ngridpoints[1]-1)*supercell[1]+1

        # Make arrays of x and y values.
        self.supercell_xarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        self.supercell_yarray = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)

        # Make arrays of supercell_density2D/I/H with the same dimensions as x/y arrays.
        self.supercell_density2D = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        self.I = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)
        self.H = np.zeros((supercell_xngridpoints,supercell_yngridpoints),np.float64)

        #Find projection of y vector onto a vector.
        ytox = np.dot(self.unit_cell[0],self.unit_cell[1].T)/self.cell_lengths[0]
        #Find component of y vector perpendicular to a vector.
        ynormal = np.cross(self.unit_cell[0],self.unit_cell[1].T)/self.cell_lengths[0]
        ynormal = np.sqrt(np.dot(ynormal,ynormal.T))

        ##################################

        for h in self.heights:
            # Contant height mode.
            if self.scan_mode == 'constant_height':
                n1 = h/self.cell_lengths[2]*self.ngridpoints[2]
                dn1 = n1-np.floor(n1)
                n1 = int(n1)%self.ngridpoints[2]
                ldos = (1-dn1)*self.density[:,:,n1]+dn1*self.density[:,:,(n1+1)%self.ngridpoints[2]]

            # Contant current mode.
            elif self.scan_mode == 'constant_current':
                n2 = h/self.cell_lengths[2]*self.ngridpoints[2]
                dn2 = n2-np.floor(n2)
                n2 = int(n2)%self.ngridpoints[2]
                # Get the averaged current.
                averaged_current = ((1-dn2)*self.density[:,:,n2].mean()+dn2*self.density[:,:,(n2+1)%self.ngridpoints[2]].mean())
                c1 = self.density[:,:,n2]
                c2 = self.density[:,:,(n2+1)%self.ngridpoints[2]]
            
            # 2D-slice at a specified height.
            elif self.scan_mode == '2Dslice':
                plane_index = int(round(h/self.cell_lengths[2]*self.ngridpoints[2]))%self.ngridpoints[2]
                density2D = self.density[:,:,plane_index]
                
            for i in range(supercell_xngridpoints):
                for j in range(supercell_yngridpoints):
                    self.supercell_xarray[i][j] = float(i)/float(supercell_xngridpoints)*self.cell_lengths[0]*m+float(j)/float(supercell_yngridpoints)*ytox*n
                    self.supercell_yarray[i][j] = float(j)/float(supercell_yngridpoints)*ynormal*n
                    mi = i%(self.ngridpoints[0]-1)
                    nj = j%(self.ngridpoints[1]-1)

                    if self.scan_mode == 'constant_height':
                        self.I[i][j] = ldos[mi][nj]

                    elif self.scan_mode == 'constant_current':
                        if c2[mi][nj]-c1[mi][nj] == 0:
                            self.H[i][j] = n2*self.cell_lengths[2]/self.ngridpoints[2]

                        else:
                            self.H[i][j] = (n2+(averaged_current-c1[mi][nj])/(c2[mi][nj]-c1[mi][nj]))*self.cell_lengths[2]/self.ngridpoints[2]
                        
                    elif self.scan_mode == '2Dslice':
                        self.supercell_density2D[i][j] = density2D[mi][nj]

    def plot(self):
        # Plot the 2D contour in matplotlib. 
	    # The newly generated images will be named after the letter "H/C/S"(constant-height mode/constant-current mode/2D-slice)
        # and height value (of tip position).
        if self.scan_mode == 'constant_height':
            P = self.I
            mode_label = "H" 
        elif self.scan_mode == 'constant_current':
            P = self.H
            mode_label = "C"
        elif self.scan_mode == '2Dslice':
            P = self.supercell_density2D
            mode_label = "S"

        plt.figure()
        plt.rcParams['figure.max_open_warning'] = 50
        plt.axis('square')
        plt.axis('off')
        plt.xticks(())
        plt.yticks(())
        cm = plt.colormaps.get_cmap('bone')
        plt.contourf(self.supercell_xarray,self.supercell_yarray,P, 900, cmap=cm)
        plt.colorbar()
        
        #plt.savefig(mode_label+'_'+str(round(h,3))+'.png', dpi=300, bbox_inches='tight')