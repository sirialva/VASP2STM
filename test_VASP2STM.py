## Use of VASP2STM script

import matplotlib.pyplot as plt
import VASPtoSTM

data = VASPtoSTM.charge_density(filename='CHG_Sr2RuO4')

data.make_cut(tip_height=23)

data.plot()
plt.show()