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

# Referenz
Referenz = {
    "$P_{HPFP}$": 1713.758302,
    "$P_{LPFP}$": 663.9537217,
    "$\dot{Q}_{FOHE}$": 122e3
}
data = {
    # Pumpe
    "Pumpe": {
        "$P_{HPFP}$": 17104.04348,
        "$P_{RV}$": 137864.765,
        "$\dot{Q}_{FOHE}$": 159e3,
        "$\dot{Q}_{PHC}$": 170232.0819
    },
    
    # Verdampfer
    "Verdampfer": {
        "$P_{HPFC}$": 38613.0293,
        "$P_{RV}$": 112593.5955,
        "$\dot{Q}_{FOHE}$": 159e3,
        "$\dot{Q}_{PHC}$": 174105.822
    },
    
    # Vormischung
    "Vormischung": {
        "$P_{HPFC}$": 43131.6977,
        "$P_{RV}$": 124787.0186,
        "$\dot{Q}_{FOHE}$": 159e3,
        "$\dot{Q}_{PHC}$": 156788.6055
    }
}

colors = ["orangered", "red", "orange", "darkred", "navy", "cyan", "dodgerblue"]
label_list = ["$P_{HPFP}$", "$P_{LPFP}$", "$P_{HPFC}$", "$P_{RV}$", "$\dot{Q}_{FOHE}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{IDG}$"]
data_list = np.array([
    [1713.758302, 663.9537217, 0, 0, 122e3, 0, 5e3],
    [17104.04348, 0, 0, 137864.765, 159e3, 170232.0819, 0],
    [0, 0, 38613.0293, 112593.5955, 159e3, 174105.822, 0],
    [0, 0, 43131.6977, 124787.0186, 159e3, 156788.6055, 0],
])/1e3


fig, ax = plt.subplots()
bottom = np.zeros(len(data_list))
for i in range(len(label_list)):
    dps = []
    for j in range(len(data_list)):
        dps.append(data_list[j][i])
    ax.bar(systemnames, dps, bottom=bottom, label=label_list[i], color=colors[i])
    bottom += dps
plt.legend()
plt.ylabel("Leistung [kW]")

fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 9/2.54)


fig.savefig(os.path.join(save_dir, 'refcomp.png'), dpi=600, bbox_inches="tight")
plt.show()
    