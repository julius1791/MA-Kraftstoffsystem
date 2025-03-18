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

subfolders = ["pump"]
systemnames = ["Pumpe"]


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


data_i = data[subfolders[0]]
P_mfp = np.array(data_i["P_mfp"])
P_fp = P_mfp

t_wu = np.array(data_i["t_wu"])
t_bk = np.array(data_i["t_bk"])


###############################################################################
################# Kraftstoffverbrauch Contour plot ##########################
###############################################################################


norm = mpl.colors.Normalize(0.111, 0.1125)
for arch, plot, name in zip(subfolders, [0], systemnames):
    
    t_wu_u = np.unique(t_wu)
    t_bk_u = np.unique(t_bk)
    P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    qm3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    i=0
    
    for t_bki in t_bk_u:
        data_i = data[arch]
        P_mfp = np.array(data_i["P_mfp"])
        P_fp = P_mfp
    
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        qm_r = np.array(data_i["qm_r"])
        qm = np.array(data_i["qm_cb"])+np.array(data_i["qm_phc"])+np.array(data_i["qm_pot"])
        Q = np.array(data_i["Q"])
        if "P_r" in data_i:
            if np.array(data_i["P_r"]).size == P_mfp.size:
                P_r = np.array(data_i["P_r"])
                P_mfp += np.array(data_i["P_r"])
        P_mfp=P_mfp/(P_mfp+Q)
        j=0
        idx = t_bk == t_bki
        qm2 = qm[idx]
        t_wu2 = t_wu[idx]
        P_mfp2 = P_mfp[idx]
        qm_r2 = qm_r[idx]
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
                qm_r3[j, i] = qm_r2[idx2][0]
                qm3[j,i] = qm2[idx2]
            except:
                qm_r3[j, i] = np.nan
                P_mfp3[j, i] = np.nan
                qm3[j,i] = np.nan
            j+=1
        i+=1
    lev = [*np.linspace(0.111, 0.1125, 400)]
    cs = plt.contourf(t_bk_u, t_wu_u, qm3, levels = lev, cmap="Reds", norm=norm)
    levs = np.linspace(0.111, 0.1125, 21)
    
    cs = plt.contour(t_bk_u, t_wu_u, qm3, levels = levs, colors='black')
    plt.clabel(cs, fontsize=10)

    #cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.08)
    #cb.set_label("$P_{mech}$ [kW]", rotation=0, labelpad=-50, loc="left")
    #clevs = [30, 40, 50, 60, 80, 100, 120, 160, 200, 250]
    plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_{W}$ [K]", fontsize=12)
    plt.xlabel("Brennkammer-Eintrittstemperatur $T_{BK}$ [K]")
    
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    plt.xlim([140, 800])
    #cb.set_ticks(clevs)
    #cb.set_ticklabels(clevs)
    plt.ticklabel_format(style="plain")
    fig.savefig(os.path.join(save_dir, name+'massflowcontour.png'), dpi=600, bbox_inches="tight")
    plt.show()


