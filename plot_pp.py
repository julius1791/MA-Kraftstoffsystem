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
    plt.xlim([max(140,t+20), 600])
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
plt.xlim([140, 600])
plt.ylim([0,250])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'pump.png'), dpi=600)
plt.show()


###############################################################################
##################### Leistungsbedarf T_W Vergleich ###########################
###############################################################################

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
plt.xlim([140, 600])
plt.ylim([0,250])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
    
plt.legend(ncols=5, handles=leg_handles, labels=leg_labels)
fig.savefig(os.path.join(save_dir, 'pmech.png'), dpi=600)
plt.show()

###############################################################################
################# Leistungsbedarf Pumpe Contour plot ##########################
###############################################################################

fig, ax = plt.subplots(3,1, sharex=True)
fig.subplots_adjust(wspace=0.1)
norm = mpl.colors.Normalize(30, 250)
for arch, plot, name in zip(subfolders, [0,1,2], systemnames):
    
    t_wu_u = np.unique(t_wu)
    t_bk_u = np.unique(t_bk)
    P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    i=0
    
    for t_bki in t_bk_u:
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
    lev = [*np.linspace(28, 250, 400)]
    cs = ax[plot].contourf(t_bk_u, t_wu_u, P_mfp3/1e3, levels = lev, cmap="Reds", norm=norm)
    ax[plot].set_title(name)
    if plot == 2:
        levs = [32, 40, 50, 60, 80, 120, 170]
    elif plot == 1:
        levs = [50, 60, 70, 80, 100, 120, 170]
    else:
        levs = [50, 60, 70, 80, 100, 170]
    
    cs = ax[plot].contour(t_bk_u, t_wu_u, P_mfp3/1e3, levels = levs, colors='black')
    plt.clabel(cs, fontsize=10)

cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.1)
cb.set_label("$P_{mech}$ [kW]", rotation=0, labelpad=-50, loc="center")
for sp in ax[0:1]:
    sp.tick_params(labelbottom=False)
clevs = [30, 40, 50, 60, 80, 100, 120, 160, 200, 250]
ax[1].set_ylabel("$T_{W}$ [K]")
ax[2].set_xlabel("$T_{BK}$ [K]")

fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 26/2.54)
plt.xlim([200,600])
cb.set_ticks(clevs)
cb.set_ticklabels(clevs)
plt.ticklabel_format(style="plain")
fig.savefig(os.path.join(save_dir, 'powercontour.png'), dpi=600, bbox_inches="tight")
plt.show()

###############################################################################
######################### Wärmebedarf Contour plot ############################
###############################################################################

fig, ax = plt.subplots(3,1, sharex=True)
fig.subplots_adjust(wspace=0.1)
norm = mpl.colors.Normalize(120, 1000)
for arch, plot, name in zip(subfolders, [0,1,2], systemnames):
    
    t_wu_u = np.unique(t_wu)
    t_bk_u = np.unique(t_bk)
    Q3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    i=0
    
    for t_bki in t_bk_u:
        data_i = data[arch]
        P_mfp = np.array(data_i["P_mfp"])
        P_fp = P_mfp
        Q = np.array(data_i["Q"])
    
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        qm_r = np.array(data_i["qm_r"])
        Q = np.array(data_i["Q"])
        
        j=0
        idx = t_bk == t_bki
        t_wu2 = t_wu[idx]
        Q2 = Q[idx]
        qm_r2 = qm_r[idx]
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                Q3[j, i] = np.extract(True, Q2[idx2])
                qm_r3[j, i] = qm_r2[idx2][0]
            except:
                qm_r3[j, i] = np.nan
                Q3[j, i] = np.nan
            j+=1
        i+=1
    lev = [*np.linspace(150, 1500, 600)]
    cs = ax[plot].contourf(t_bk_u, t_wu_u, Q3/1e3, levels = lev, cmap="Reds", norm=norm)
    ax[plot].set_title(name)

    levs = [200, 300, 400, 500, 600, 700, 800]

    
    cs = ax[plot].contour(t_bk_u, t_wu_u, Q3/1e3, levels = levs, colors='black')
    plt.clabel(cs, fontsize=10)

cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.1)
cb.set_label("$Q_{gsmt}$ [kW]", rotation=0, labelpad=-50, loc="center")
for sp in ax[0:1]:
    sp.tick_params(labelbottom=False)
clevs = [120, 200, 300, 400, 500, 600, 700, 800, 1000]
ax[1].set_ylabel("$T_{W}$ [K]")
ax[2].set_xlabel("$T_{BK}$ [K]")

fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 26/2.54)
plt.xlim([200,600])
cb.set_ticks(clevs)
cb.set_ticklabels(clevs)
plt.ticklabel_format(style="plain")
fig.savefig(os.path.join(save_dir, 'heatcontour.png'), dpi=600, bbox_inches="tight")
plt.show()

###############################################################################
################### Massenstromverhältnis Pumpe ###############################
###############################################################################


lev = [*np.linspace(0.1, 10, 800)]
cs = plt.contourf(t_bk_u, t_wu_u, qm_r3/0.1, levels=lev, cmap="Reds", norm=mpl.colors.LogNorm(0.08, 10))
cb = plt.colorbar(label="$\dot{m}_r/\dot{m}_{BK}$ [-]")
clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3, 5, 10]
cb.set_ticks(clevs)
cb.set_ticklabels(clevs)
plt.xlabel("$T_{BK}$ [K]")
plt.ylabel("$T_{W}$ [K]")
clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3]
cs = plt.contour(t_bk_u, t_wu_u, qm_r3/0.1, levels = clevs, colors='black')
plt.clabel(cs, fontsize=8)
plt.title("Architektur: " + arch)
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)

plt.xlim([200,600])
fig.savefig(os.path.join(save_dir, 'pump_massflow.png'), dpi=600)
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
            
    id_200 = t_wu == 140
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
    Q_hx2 = [200 for i in Q_200]
    
    
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
    plt.xlim([160,400])
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
            
    id_200 = t_wu == 140
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
    plt.xlim([160,400])
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




