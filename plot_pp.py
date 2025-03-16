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

tbk_lims = [150,400]
twu_lims = [100,280]

folder = os.path.join(os.getcwd(), "results")
save_dir = os.path.join(os.getcwd(), "diagrams")

subfolders = ["dual", "after", "pump"]
systemnames = ["Vormischung", "Verdampfer", "Pumpe"]

###############################################################################
############################# Daten Sammeln ###################################
###############################################################################

data = dict()
for subfolder in subfolders:
    file_path = os.path.join(folder, subfolder + ".json")
    name = subfolder
    with open(file_path) as jsonfile:
        data.update({name: json.load(jsonfile)})
stylelist = ["-", ":", "--", "-."]
colors = ["red", "cyan", "lawngreen", "navy", "orangered", "green"]

###############################################################################
######################## Leistungsbedarf Systemvergleich ######################
###############################################################################

t_list = [100, 160, 220, 280]
for t in t_list:
    for key, style, name, color in zip(data, stylelist, systemnames, colors):
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
        plt.plot(t_bk1, P_mfp1/1e3, label=name, color=color, linestyle=style)
    plt.legend(title="Kraftstoffsystem", loc="best")
    plt.xlabel("$T_{BK}$ [K]")
    plt.ylabel("$P_{mech}$ [kW]")
    plt.xlim([max(140,t+20), tbk_lims[1]])
    plt.ylim([0,300])
    plt.title("$T_{W}$ = " + str(t) + " K")
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    fig.savefig(os.path.join(save_dir, 'arch_' + str(t) + '.png'), dpi=600, bbox_inches="tight")
    plt.show()
    

###############################################################################
################# Leistungsbedarf Pumpe T_W Vergleich #########################
###############################################################################


data_i = data[subfolders[2]]
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
t_hx_list = reversed([100, 160, 220, 280])
for t_hx, style, color in zip(t_hx_list, stylelist, colors):
    
    idx = t_wu == t_hx
    t_bk1 = t_bk[idx]
    P_mfp1 = P_mfp[idx]
    t_wu1 = t_wu[idx]
    plt.plot(t_bk1, P_mfp1/1e3, label=str(t_hx) + " K", linestyle=style, color=color)
plt.legend(title="$T_{W}$ [K]")
plt.xlabel("$T_{BK}$ [K]")
plt.ylabel("$P_{mech}$ [kW]")
plt.title("Architektur: Pumpe")
plt.xlim(tbk_lims)
plt.ylim([0,300])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'pump.png'), dpi=600)
plt.show()


###############################################################################
##################### Leistungsbedarf T_W Vergleich ###########################
###############################################################################

fig, ax = plt.subplots()
for arch, color in zip(subfolders, colors):
    data_i = data[arch]
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
    t_hx_list = [280, 220, 160, 100]
    for t_hx, style in zip(t_hx_list, stylelist):
        
        idx = t_wu == t_hx
        t_bk1 = t_bk[idx]
        P_mfp1 = P_mfp[idx]
        t_wu1 = t_wu[idx]
        plt.plot(t_bk1, P_mfp1/1e3, linestyle=style, color=color)
tit = Rectangle([0,0], 1, 0.1, fc="none", ec="none")
lvm = Rectangle([0,0], 1, 1, fc=colors[0])
lvd = Rectangle([0,0], 1, 1, fc=colors[1])
lp = Rectangle([0,0], 1, 1, fc=colors[2])
l100 = Line2D([0], [0], color="black", linestyle=stylelist[0])
l160 = Line2D([0], [0], color="black", linestyle=stylelist[1])
l220 = Line2D([0], [0], color="black", linestyle=stylelist[2])
l280 = Line2D([0], [0], color="black", linestyle=stylelist[3])
leg_handles = [tit, tit, l100, lvm, l160, lvd, l220, lp, l280]
leg_labels = ["$T_{W}$ [K]", "Architektur", t_hx_list[0], subfolders[0], t_hx_list[1], subfolders[1], t_hx_list[2], subfolders[2], t_hx_list[3]]
plt.xlabel("$T_{BK}$ [K]")
plt.ylabel("$P_{mech}$ [kW]")
plt.xlim(tbk_lims)
plt.ylim([0,300])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 12/2.54)
plt.subplots_adjust(top=0.7)
plt.legend(ncols=5, handles=leg_handles, labels=leg_labels, bbox_to_anchor=(0.5,1.15), loc="upper center")
fig.savefig(os.path.join(save_dir, 'pmech.png'), dpi=600)
plt.show()


    

###############################################################################
########################### Energie Stackplots ################################
###############################################################################


for subfolder,name in zip(subfolders, systemnames):
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
            
    id_200 = t_wu == 120
    t_bk_200 = t_bk[id_200]
    P_fp_200 = P_fp[id_200]/1e3
    if subfolder != "pre":
        P_r = P_r[id_200]/1e3
        P_fp_200 -= P_r
    else:
        P_r = P_fp_200 * 0
    Q_200 = Q[id_200]/1e3
    qm_200 = qm_cb[id_200] + qm_phc[id_200]
    Q_hx = Q_200.clip(0, 159)
    Q_phc = Q_200 - Q_hx
    Q_excess = Q_200 - 159
    Q_excess = Q_excess.clip(-159, 0)
    Q_hx2 = [159 for i in Q_200]
    
    
    spQ = plt.stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
    #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
    spP = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
    p1 = Rectangle((0, 0), 1, 1, fc="navy")
    p2 = Rectangle((0, 0), 1, 1, fc="cyan")
    p3 = Rectangle((0, 0), 1, 1, fc="darkred")
    p4 = Rectangle((0, 0), 1, 1, fc="orangered")
    if subfolder != "pump":
        plt.legend([p3, p4, p2, p1], ["$P_{HPFC}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)
    else:
        plt.legend([p3, p4, p2, p1], ["$P_{HPFP}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)

    plt.xlabel("$T_{BK}$ [K]")
    plt.ylabel("Leistung [kW]")
    plt.title("Architektur: " + name)
    plt.xlim([150, 320])
    plt.ylim([0, 500])
    p = 1.7e6
    h = list()
    # for t in t_bk_200:
    #     qmi = qm_200[t_bk_200 == t]
    #     Q_hxi = Q_hx[t_bk_200 == t]/1e3
    #     Q_phci = Q_phc[t_bk_200 == t]/1e3
    #     P_fpi = P_fp_200[t_bk_200 == t]/1e3
    #     P_ri = P_r[t_bk_200 == t]/1e3
    #     hi = (h2.calc_H2_enthalpy(t, p) - h2.calc_H2_enthalpy(22, 4.2e5))/1e3 * qmi
    #     #print((hi-Q_phci - Q_hxi-P_fpi-P_ri)/hi)
    #     h.append(hi)
        
    # plt.plot(t_bk_200, h)
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    fig.savefig(os.path.join(save_dir, 'stackplot_' + subfolder + '.png'), dpi=600, bbox_inches="tight")
    plt.show()
    
    
###############################################################################
########################### Energie Stackplots ################################
###############################################################################


for subfolder,name in zip(subfolders, systemnames):
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
            
    id_200 = t_wu == t_bk - 40
    t_bk_200 = t_bk[id_200]
    P_fp_200 = P_fp[id_200]/1e3
    if subfolder != "pre":
        P_r = P_r[id_200]/1e3
        P_fp_200 -= P_r
    else:
        P_r = P_fp_200 * 0
    Q_200 = Q[id_200]/1e3
    qm_200 = qm_cb[id_200] + qm_phc[id_200]
    Q_hx = Q_200.clip(0, 159)
    Q_phc = Q_200 - Q_hx
    Q_excess = Q_200 - 159
    Q_excess = Q_excess.clip(-159, 0)
    Q_hx2 = [159 for i in Q_200]
    
    
    spQ = plt.stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
    #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
    spP = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
    p1 = Rectangle((0, 0), 1, 1, fc="navy")
    p2 = Rectangle((0, 0), 1, 1, fc="cyan")
    p3 = Rectangle((0, 0), 1, 1, fc="darkred")
    p4 = Rectangle((0, 0), 1, 1, fc="orangered")
    if subfolder != "pump":
        plt.legend([p3, p4, p2, p1], ["$P_{HPFC}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)
    else:
        plt.legend([p3, p4, p2, p1], ["$P_{HPFP}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)

    plt.xlabel("$T_{BK}$ [K]")
    plt.ylabel("Leistung [kW]")
    plt.title("Architektur: " + name)
    plt.xlim([150, 320])
    plt.ylim([0, 500])
    p = 1.7e6
    h = list()
    # for t in t_bk_200:
    #     qmi = qm_200[t_bk_200 == t]
    #     Q_hxi = Q_hx[t_bk_200 == t]/1e3
    #     Q_phci = Q_phc[t_bk_200 == t]/1e3
    #     P_fpi = P_fp_200[t_bk_200 == t]/1e3
    #     P_ri = P_r[t_bk_200 == t]/1e3
    #     hi = (h2.calc_H2_enthalpy(t, p) - h2.calc_H2_enthalpy(22, 4.2e5))/1e3 * qmi
    #     #print((hi-Q_phci - Q_hxi-P_fpi-P_ri)/hi)
    #     h.append(hi)
        
    # plt.plot(t_bk_200, h)
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    fig.savefig(os.path.join(save_dir, 'stackplot_' + subfolder + '.png'), dpi=600, bbox_inches="tight")
    plt.show()
    
###############################################################################
######################## Leistungssplit vs Rezirkulation ######################
###############################################################################




t = 300

for key, style, name, color in zip(data, stylelist, systemnames, colors):
    data_i = data[key]
    P_mfp = np.array(data_i["P_mfp"])
    t_bk = np.array(data_i["t_bk"])
    t_wu = np.array(data_i["t_wu"])
    qm_r = np.array(data_i["qm_r"])
    if "P_r" in data_i:
        if np.array(data_i["P_r"]).size == P_mfp.size:
            P_mfp += np.array(data_i["P_r"])
            RV_ratio = np.array(data_i["P_r"])/P_mfp
    idx = t_bk == t
    t_wu=t_wu[idx]
    qm_r = qm_r[idx]
    RV_ratio = RV_ratio[idx]
    plt.plot(t-t_wu, RV_ratio, color=color, label=name)
plt.legend()
plt.xlabel("Temperaturdifferenz $T_{BK}-T_{W}$ [K]")
plt.ylabel("$P_{RV}/P_{mech}$ [kW]")
plt.xlim([20,200])
plt.ylim([0,1])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'rezirkulation.png'), dpi=600, bbox_inches="tight")
plt.show()




