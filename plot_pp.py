# -*- coding: utf-8 -*-
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import h2_properties as h2

folder = os.path.join(os.getcwd(), "results2")

subfolders = ["pump", "after", "pre", "dual"]

data = dict()
for subfolder in subfolders:
    file_path = os.path.join(folder, subfolder + ".json")
    name = subfolder
    with open(file_path) as jsonfile:
        data.update({name: json.load(jsonfile)})

t_list = [100, 160, 220, 280]
for t in t_list:
    for key in data:
        data_i = data[key]
        P_mfp = np.array(data_i["P_mfp"])
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        if "P_r" in data_i:
            if np.array(data_i["P_r"]).size == P_mfp.size:
                P_mfp += np.array(data_i["P_r"])
        idx = t_wu == t
        t_bk1 = t_bk[idx]
        P_mfp1 = P_mfp[idx]
        t_wu1 = t_wu[idx]
        plt.plot(t_bk1, P_mfp1/1e3, label=key)
    plt.legend()
    plt.xlabel("Brennkammereintrittstemperatur [K]")
    plt.ylabel("Leistungsbedarf [kW]")
    plt.title("Wärmetauschereintrittstemperatur = " + str(t) + " [K]")
    plt.show()



data_i = data["pump"]
P_mfp = np.array(data_i["P_mfp"])
P_fp = P_mfp

t_wu = np.array(data_i["t_wu"])
t_bk = np.array(data_i["t_bk"])
qm_r = np.array(data_i["qm_r"])
Q = np.array(data_i["Q"])
if "P_r" in data_i:
    if np.array(data_i["P_r"]).size == P_mfp.size:
        P_r = np.array(data_i["P_r"])
        P_mfp += np.array(data_i["P_r"])
t_hx_list = [100, 160, 220, 280]
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
    idx = t_bk == t_bki
    t_wu2 = t_wu[idx]
    P_mfp2 = P_mfp[idx]
    qm_r2 = qm_r[idx]
    for t_wui in t_wu_u:
        idx2 = t_wu2 == t_wui
        try:
            P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
            qm_r3[j, i] = qm_r2[idx2][0]
        except:
            qm_r3[j, i] = np.nan
            P_mfp3[j, i] = np.nan
        j+=1
    i+=1


cs = plt.contour(t_bk_u, t_wu_u, P_mfp3/1e3, levels=[*np.linspace(0, 6.5, 14), *np.linspace(7, 11, 5), *np.linspace(12, 20, 5), 25, *np.linspace(30, 80, 6)])
plt.clabel(cs)
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Wärmetauschereintrittstemperatur [K]")
plt.title("Pumpenleistung [kW]")
plt.show()

cs = plt.contour(t_bk_u, t_wu_u, qm_r3/0.1, levels=[*np.linspace(0.0, 0.05, 6)/0.1, *np.linspace(0.06, 0.12, 4)/0.1, *np.linspace(0.15, 0.3, 4)/0.1, *np.linspace(0.4, 1, 4)/0.1])
plt.clabel(cs)
plt.xlabel("Brennkammereintrittstemperatur [K]")
plt.ylabel("Wärmetauschereintrittstemperatur [K]")
plt.title("Rezirkulationsmassenstrom/Brennkammermassenstrom [-]")
plt.show()


for subfolder in subfolders:
    data_i = data[subfolder]
    P_mfp = np.array(data_i["P_mfp"])
    P_fp = P_mfp
    
    t_wu = np.array(data_i["t_wu"])
    t_bk = np.array(data_i["t_bk"])
    qm_r = np.array(data_i["qm_r"])
    qm_cb = np.array(data_i["qm_cb"])
    qm_phc = np.array(data_i["qm_phc"])
    Q = np.array(data_i["Q"])
    if "P_r" in data_i:
        if np.array(data_i["P_r"]).size == P_mfp.size:
            P_r = np.array(data_i["P_r"])
            P_mfp += np.array(data_i["P_r"])
            
    id_200 = t_wu == 180
    t_bk_200 = t_bk[id_200]
    P_fp_200 = P_fp[id_200]/1e3
    if subfolder != "pre":
        P_r = P_r[id_200]/1e3
        P_fp_200 -= P_r
    else:
        P_r = P_fp_200 * 0
    Q_200 = Q[id_200]/1e3
    qm_200 = qm_cb[id_200] + qm_phc[id_200]
    Q_hx = Q_200.clip(0, 200)
    Q_phc = Q_200 - Q_hx
    Q_excess = Q_200 - 200
    Q_excess = Q_excess.clip(-200, 0)
    sp = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_fp_200, P_r, baseline="zero")
    sp = plt.stackplot(t_bk_200, Q_excess, baseline="zero")
    plt.xlabel("Brennkammereintrittstemperatur [K]")
    plt.ylabel("Leistung [kW]")
    plt.title(subfolder)
    plt.xlim([200,600])
    plt.ylim([-10, 600])
    p = 1.7e6
    h = list()
    for t in t_bk_200:
        qmi = qm_200[t_bk_200 == t]
        Q_hxi = Q_hx[t_bk_200 == t]/1e3
        Q_phci = Q_phc[t_bk_200 == t]/1e3
        P_fpi = P_fp_200[t_bk_200 == t]/1e3
        P_ri = P_r[t_bk_200 == t]/1e3
        hi = (h2.calc_H2_enthalpy(t, p) - h2.calc_H2_enthalpy(22, 4.2e5))/1e3 * qmi
        #print((hi-Q_phci - Q_hxi-P_fpi-P_ri)/hi)
        h.append(hi)
        
    plt.plot(t_bk_200, h)
    plt.show()


