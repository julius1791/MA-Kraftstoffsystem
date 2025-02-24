# -*- coding: utf-8 -*-
import json
import os
import matplotlib.pyplot as plt
import numpy as np

folder = os.path.join(os.getcwd(), "results")

subfolders = ["pump_test", "after_test", "pre_test", "dual_test"]

data = dict()
for subfolder in subfolders:
    file_path = os.path.join(folder, subfolder + ".json")
    name = subfolder[:-5]
    with open(file_path) as jsonfile:
        data.update({name: json.load(jsonfile)})

for key in data:
    data_i = data[key]
    P_mfp = np.array(data_i["P_mfp"])
    t_wu = np.array(data_i["t_wu"])
    t_bk = np.array(data_i["t_bk"])
    if "P_r" in data_i:
        if np.array(data_i["P_r"]).size == P_mfp.size:
            P_mfp += np.array(data_i["P_r"])
    idx = t_wu == 260
    t_bk1 = t_bk[idx]
    P_mfp1 = P_mfp[idx]
    t_wu1 = t_wu[idx]
    plt.plot(t_bk1, P_mfp1/1e3, label=key)
plt.legend()
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Leistungsbedarf [kW]")
plt.title("Wärmetauschereintrittstemperatur = 260 [K]")
plt.show()

data_i = data["pump"]
P_mfp = np.array(data_i["P_mfp"])
t_wu = np.array(data_i["t_wu"])
t_bk = np.array(data_i["t_bk"])
qm_r = np.array(data_i["qm_r"])
if "P_r" in data_i:
    if np.array(data_i["P_r"]).size == P_mfp.size:
        P_mfp += np.array(data_i["P_r"])
t_hx_list = range(100, 300, 40)
for t_hx in t_hx_list:
    
    idx = t_wu == t_hx
    t_bk1 = t_bk[idx]
    P_mfp1 = P_mfp[idx]
    t_wu1 = t_wu[idx]
    plt.plot(t_bk1, P_mfp1/1e3, label=t_hx)
plt.legend()
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Leistungsbedarf [kW]")
plt.title("pump für unterschiedliche HX Temperaturen [K]")
plt.show()

t_wu_u = np.unique(t_wu)
t_bk_u = np.unique(t_bk)
P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
i=0

for t_bki in t_bk_u:
    j=0
    for t_wui in t_wu_u:
        
        idx = t_bk == t_bki
        t_wu2 = t_wu[idx]
        P_mfp2 = P_mfp[idx]
        qm_r2 = qm_r[idx]
        idx2 = t_wu2 == t_wui
        try:
            P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
            qm_r3[j, i] = qm_r2[idx2][0]
        except:
            qm_r3[j, i] = np.nan
        j+=1
    i+=1


cs = plt.contour(t_bk_u, t_wu_u, P_mfp3/1e3, levels=[*range(0, 20, 5), *range(20, 100, 10), *range(100, 200, 20)])
plt.clabel(cs, levels=[*range(0, 20, 5), *range(20, 100, 10), *range(100, 200, 20)])
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Wärmetauschereintrittstemperatur [K]")
plt.title("Pumpenleistung [kW]")
plt.show()

cs = plt.contour(t_bk_u, t_wu_u, qm_r3/0.1, levels=[*np.linspace(0, 0.06, 6)/0.1, *np.linspace(0.08, 0.3, 6)/0.1, *np.linspace(0.4, 1, 4)/0.1])
plt.clabel(cs)
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Wärmetauschereintrittstemperatur [K]")
plt.title("Rezirkulationsmassenstrom/Brennkammermassenstrom [-]")
plt.show()