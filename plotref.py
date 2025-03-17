# -*- coding: utf-8 -*-
import json
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# import h2_properties as h2
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.cm as cm

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

systemnames = ["Referenz", "Pumpe", "Verdampfer", "Vormischung"]
save_dir = os.path.join(os.getcwd(), "diagrams")

# # Referenz
# Referenz = {
#     "$P_{HPFP}$": 1713.758302,
#     "$P_{LPFP}$": 663.9537217,
#     "$\dot{Q}_{FOHE}$": 122e3
# }
# data = {
#     # Pumpe
#     "Pumpe": {
#         "$P_{HPFP}$": 17104.04348,
#         "$P_{RV}$": 137864.765,
#         "$\dot{Q}_{FOHE}$": 159e3,
#         "$\dot{Q}_{PHC}$": 170232.0819
#     },
    
#     # Verdampfer
#     "Verdampfer": {
#         "$P_{HPFC}$": 38613.0293,
#         "$P_{RV}$": 112593.5955,
#         "$\dot{Q}_{FOHE}$": 159e3,
#         "$\dot{Q}_{PHC}$": 174105.822
#     },
    
#     # Vormischung
#     "Vormischung": {
#         "$P_{HPFC}$": 43131.6977,
#         "$P_{RV}$": 124787.0186,
#         "$\dot{Q}_{FOHE}$": 159e3,
#         "$\dot{Q}_{PHC}$": 156788.6055
#     }
# }

colors = ["darkred", "orangered", "navy", "cyan"]
label_list = ["$P_{HPFP/C}$", "$P_{RV/LPFP}$", "$\dot{Q}_{FOHE}$", "$\dot{Q}_{PHC}$"]
data_list = np.array([
    [1713.758302, 663.9537217, 127e3, 0],
    [17104.04348, 137864.765, 159e3, 170232.0819],
    [38613.0293, 112593.5955, 159e3, 174105.822],
    [43131.6977, 124787.0186, 159e3, 156788.6055],
])/1e3


fig, ax = plt.subplots(1, 2, width_ratios=[1, 5])
fig.subplots_adjust(wspace=0.3)
bottom = np.zeros(len(data_list))
for i in range(len(label_list)):
    dps = []
    for j in range(len(data_list)):
        dps.append(data_list[j][i])
    ax[1].bar(systemnames, dps, bottom=bottom, label=label_list[i], color=colors[i], width = 0.5)
    bottom += dps

ax[1].set_ylabel("Leistung [kW]")
ax[1].set_ylim([0, 500])

bottom = 0
for i in range(len(label_list)):
    dps = data_list[0][i]
    ax[0].bar(systemnames[0], dps, bottom=bottom, label=label_list[i], color=colors[i], width=0.5)
    bottom += dps
ax[0].set_ylabel("Leistung [kW]")
ax[0].set_ylim([0, 9.9])
ax[0].legend(bbox_to_anchor=[0.5, 1.4], loc="upper center", ncols=2, columnspacing=0.5, fontsize=12)
ax[0].set_position([0, 0.05, 0.1, 0.6])
ax[1].set_position([0.3, 0.05, 0.65, 0.9])


fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 9/2.54)


fig.savefig(os.path.join(save_dir, 'refcomp.png'), dpi=600, bbox_inches="tight")
plt.show()
    