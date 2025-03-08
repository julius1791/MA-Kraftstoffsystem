# -*- coding: utf-8 -*-
import json
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# import h2_properties as h2
from matplotlib.patches import Rectangle

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

folder = os.path.join(os.getcwd(), "results2")
save_dir = os.path.join(os.getcwd(), "diagrams")

subfolders = ["pre", "dual", "after", "pump"]

data = dict()
for subfolder in subfolders:
    file_path = os.path.join(folder, subfolder + ".json")
    name = subfolder
    with open(file_path) as jsonfile:
        data.update({name: json.load(jsonfile)})
stylelist = ["-", ":", "--", "-."]

t_list = [100, 160, 220, 280]
for t in t_list:
    for key, style in zip(data, stylelist):
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
        plt.plot(t_bk1, P_mfp1/1e3, label=key, color="black", linestyle=style)
    plt.legend(title="Architecture", loc="upper right")
    plt.xlabel("$T_{BK,ein}$ [K]")
    plt.ylabel("$P_{mech}$ [kW]")
    plt.title("$T_{HX,ein}$ = " + str(t) + " K")
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    fig.savefig(os.path.join(save_dir, 'arch_' + str(t) + '.png'), dpi=600, bbox_inches="tight")
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
t_hx_list = reversed([100, 160, 220, 280])
for t_hx, style in zip(t_hx_list, stylelist):
    
    idx = t_wu == t_hx
    t_bk1 = t_bk[idx]
    P_mfp1 = P_mfp[idx]
    t_wu1 = t_wu[idx]
    plt.plot(t_bk1, P_mfp1/1e3, label=str(t_hx) + " K", color="black", linestyle=style)
plt.legend(title="$T_{HX,ein}$ [K]")
plt.xlabel("$T_{BK,ein}$ [K]")
plt.ylabel("$P_{mech}$ [kW]")
plt.title("Architektur: pump [K]")
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'pump.png'), dpi=600)
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
lev = [*np.linspace(4, 75, 400)]
cs = plt.contourf(t_bk_u, t_wu_u, P_mfp3/1e3, levels = lev, cmap="Greys", norm=mpl.colors.LogNorm(2, 75))
cb = plt.colorbar(label="$P_{mech}$ [kW]")
clevs = [4, 5, 6, 8, 10, 12, 15, 20, 30, 50, 75]
cb.set_ticks(clevs)
cb.set_ticklabels(clevs)
plt.xlabel("$T_{BK,ein}$ [K]")
plt.ylabel("$T_{HX,ein}$ [K]")
levs = [4, 5, 6, 8, 10, 12, 15, 20]
cs = plt.contour(t_bk_u, t_wu_u, P_mfp3/1e3, levels = levs, colors='black')
plt.clabel(cs, fontsize=8)
plt.title("Leistungsbedarf [Architektur: pump]")
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'pump_power.png'), dpi=600)
plt.show()

lev = [*np.linspace(0.1, 10, 800)]
cs = plt.contourf(t_bk_u, t_wu_u, qm_r3/0.1, levels=lev, cmap="Greys", norm=mpl.colors.LogNorm(0.08, 10))
cb = plt.colorbar(label="$\dfrac{\dot{m}_r}{\dot{m}_{bk}}$ [-]")
clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3, 5, 10]
cb.set_ticks(clevs)
cb.set_ticklabels(clevs)
plt.xlabel("$T_{BK,ein}$ [K]")
plt.ylabel("$T_{HX,ein}$ [K]")
clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3]
cs = plt.contour(t_bk_u, t_wu_u, qm_r3/0.1, levels = clevs, colors='black')
plt.clabel(cs, fontsize=8)
plt.title("Massenstromverhältnis [Architektur: pump]")
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'pump_massflow.png'), dpi=600)
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
    
    spP = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['white', 'white', 'black', 'grey'], baseline="zero")
    spQ = plt.stackplot(t_bk_200, Q_hx, Q_phc, colors=['white', 'white'], baseline="zero", hatch=['//', "+"], edgecolors=["black", "black"])

    p1 = Rectangle((0, 0), 1, 1, fc="white", hatch="//", ec="black")
    p2 = Rectangle((0, 0), 1, 1, fc="white", hatch="+", ec="black")
    p3 = Rectangle((0, 0), 1, 1, fc="grey")
    p4 = Rectangle((0, 0), 1, 1, fc="black")
    if subfolder != "pre":
        plt.legend([p3, p4, p2, p1], ["$P_{M}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{HX}$"], loc="lower right", framealpha=1)
    else:
        plt.legend([p3, p2, p1], ["$P_{M}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{HX}$"], loc="lower right", framealpha=1)

    plt.xlabel("$T_{BK,ein}$ [K]")
    plt.ylabel("Leistung [kW]")
    plt.title("Architektur: " + subfolder)
    plt.xlim([200,700])
    plt.ylim([0, 800])
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


