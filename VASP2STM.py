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

# read CHGCAR
vasp_charge = VaspChargeDensity('CHGCAR')
density = vasp_charge.chg[-1]
atoms = vasp_charge.atoms[-1]
del vasp_charge