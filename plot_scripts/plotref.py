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
import pathlib
import csv

cwd = os.getcwd()
os.chdir(os.path.dirname(cwd))

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
plt.rcParams['hatch.color'] = 'gray'
plt.rcParams['hatch.linewidth'] = 4
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

systemnames = ["Referenz", "Pumpe", "Verdampfer", "Vormischung"]
save_dir = os.path.join(os.getcwd(), "diagrams")
folder = os.path.join(os.getcwd(), "single_results")

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

edgecols = ["none", "none", "none", "none", "red", "none"]
hatches = ["", "", "", "", "//", ""]
colors = ["darkred", "yellow", "orangered", "navy", "navy", "cyan"]
label_list = ["$P_\mathrm{HP}$", "$P_\mathrm{LP}$", "$P_\mathrm{RV}$", "$\dot{Q}_\mathrm{FOHE}$", "Überschuss", "$\dot{Q}_\mathrm{PHC}$"]
# # tw=250
# data_list = np.array([
#     [1713.758305, 663.9556996, 127e3, 0],
#     [17415.53696, 140419.1372, 159e3, 176184.3309],
#     [39300.28627, 114551.1699, 159e3, 180081.0622],
#     [43978.56931, 127221.6915, 159e3, 162949.7701],
# ])/1e3

# tw=160
data_list = np.array([
    [1713.758305, 663.9556996, 127e3, 0],
    [17347.81637, 36680.74905, 159e3, 278073.3287],
    [39155.20789, 25198.42996, 159e3, 267759.2037],
    [43819.02388, 27917.47608, 159e3, 260612.5962],
])/1e3


# find files in target directory
files = [f for f in os.listdir(folder) if os.path.isfile(
    os.path.join(pathlib.Path().resolve(), folder, f))]
data = dict()
for file in files:
    # only import .csv files
    if file.split(".")[1] == "csv":
        filename = os.path.join(folder, file)
        with open(filename, newline="") as f:
            data_f = csv.reader(f)
            _ = next(data_f)
            param_row = next(data_f)
            _ = next(data_f)
            P_row = next(data_f)
            _ = next(data_f)
            Q_row = next(data_f)
            _ = next(data_f)
            qm_row = next(data_f)
        func = param_row[0]
        data.update({param_row[0]: {
            "P_hp": float(P_row[0]), "P_2": float(P_row[1]), 
            "Q_fohe": float(Q_row[0])-float(Q_row[1]), 
            "Q_phc": float(Q_row[1]), "qm_phc": float(qm_row[1])*1e3}
            })


ref = data["reference"]
ref2 = data["reference2"]
pump = data["pump"]
after = data["after"]
dual = data["dual"]

data_list = np.array([
    [ref["P_hp"], ref["P_2"], 0, ref2["Q_fohe"], ref["Q_fohe"]-ref2["Q_fohe"], 0],
    [pump["P_hp"], 0, pump["P_2"], pump["Q_fohe"], 0, pump["Q_phc"]],
    [after["P_hp"], 0, after["P_2"], after["Q_fohe"], 0, after["Q_phc"]],
    [dual["P_hp"], 0, dual["P_2"], dual["Q_fohe"], 0, dual["Q_phc"]],
])/1e3

qm_phc_list = [0, pump["qm_phc"], after["qm_phc"], dual["qm_phc"]]

fig, ax = plt.subplots(1, 2, width_ratios=[1, 5])
ax2 = ax[1].twinx()
fig.subplots_adjust(wspace=0.3)
bottom = np.zeros(len(data_list))
xvals = np.arange(0, len(data_list), 1)-np.full(len(data_list), 0.05)
xvals[0] = 0
for i in range(len(label_list)):
    dps = []
    for j in range(len(data_list)):
        dps.append(data_list[j][i])

    ax[1].bar(xvals, dps, bottom=bottom, label=label_list[i], fc=colors[i], hatch=hatches[i], width = 0.5)
    bottom += dps
ax2.bar(np.arange(0, len(data_list), 1)+np.full(len(data_list), 0.25), qm_phc_list, label="$\dot{m}_\mathrm{PHC}$", fc = "green", width = 0.1)
ax[1].set_ylabel("Leistung [kW]")
ax2.set_ylabel("PHC Wasserstoffbedarf $\dot{m}_\mathrm{PHC}$ [g/s]")
ax[1].set_ylim([0, 500])
ax[1].set_xticks(np.arange(0, len(data_list), 1), systemnames)

bottom = 0
ax[0].bar(0, 0, bottom=0, label="$\dot{m}_\mathrm{PHC}$", fc = "green", width = 0)
for i in range(len(label_list)):
    dps = data_list[0][i]
    ax[0].bar(systemnames[0], dps, bottom=bottom, label=label_list[i], fc=colors[i], hatch=hatches[i], width=0.4)
    bottom += dps

ax[0].set_ylabel("Leistung [kW]")
ax[0].set_ylim([0, 3.6])
ax[0].legend(bbox_to_anchor=[0.5, 1.55], loc="upper center", ncols=2, columnspacing=0.5, fontsize=12)
ax[0].set_position([0, 0.05, 0.1, 0.6])
ax[1].set_position([0.3, 0.05, 0.65, 0.9])


fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 9/2.54)


fig.savefig(os.path.join(save_dir, 'refcomp.pdf'), dpi=600, bbox_inches="tight")
plt.show()


# H_h2 = 119.96 #MJ/kg
# H_jeta = 43.1 # MJ/kg

# energy = np.array([0.311338*H_jeta, 0.11155807*H_h2, 0.111560557*H_h2, 0.111614224*H_h2]) - np.array([0.311338*H_jeta, 0.311338*H_jeta, 0.311338*H_jeta, 0.311338*H_jeta])




# fig, ax = plt.subplots(1, 2, width_ratios=[1, 5])
# ax[1].bar(systemnames, energy, bottom=[0,0,0,0], width=0.5, color="navy")
# ax[1].set_ylabel("Energieverbrauch Differenz $\Delta \dot{H}_\mathrm{ref}$ [MW]")
# ax[0].bar("Referenz", 0.311338*H_jeta, bottom=0, width=0.5, color="orangered")
# ax[0].set_ylabel("Energieverbrauch $\dot{H}$ [MW]")
# fig.set_size_inches(16/2.54, 9/2.54)
# ax[0].set_position([0, 0.05, 0.1, 0.9])
# ax[1].set_position([0.3, 0.05, 0.65, 0.9])
# fig.savefig(os.path.join(save_dir, 'refenergy.pdf'), dpi=600, bbox_inches="tight")
# plt.show()

