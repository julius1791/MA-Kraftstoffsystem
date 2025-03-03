# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
import h2_properties as h2
csfont = {'family': "serif", "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=False)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

save_dir = os.path.join(os.getcwd(), "diagrams")

h1 = h2.calc_H2_enthalpy(280, )