import json
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# import h2_properties as h2
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
import matplotlib.cm as cm

parent_folder = os.path.dirname(os.getcwd())

folder = os.path.join(parent_folder, "results")
save_dir = os.path.join(parent_folder, "diagrams")

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

tbk_lims = [150,500]
twu_lims = [100,280]



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


data_i = data[subfolders[0]]
P_mfp = np.array(data_i["P_mfp"])
P_fp = P_mfp

t_wu = np.array(data_i["t_wu"])
t_bk = np.array(data_i["t_bk"])


# Define circle positions and corresponding letters
circle_positions = [(300, 160), (180, 160), (495, 105), (300, 275)]
letters = ['A', 'B', 'C', 'D']



###############################################################################
################# Leistungsbedarf Pumpe Contour plot ##########################
###############################################################################


norm = mpl.colors.Normalize(30, 300)
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
    lev = [*np.linspace(28, 300, 400)]
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.contourf(t_bk_u, t_wu_u, P_mfp3/1e3, levels = lev, cmap="Reds", norm=norm, zorder=0)
    if plot == 2:
        levs = [32, 40, 50, 60, 80, 120, 170, 250]
    elif plot == 1:
        levs = [50, 60, 70, 80, 100, 120, 170, 250]
    else:
        levs = [50, 60, 70, 80, 100, 120, 170, 250]
    
    cs = ax.contour(t_bk_u, t_wu_u, P_mfp3/1e3, levels = levs, colors='black', zorder=5)
    plt.clabel(cs, fontsize=10)

    #cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.08)
    #cb.set_label("$P_{mech}$ [kW]", rotation=0, labelpad=-50, loc="left")
    #clevs = [30, 40, 50, 60, 80, 100, 120, 160, 200, 250]
    plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]", fontsize=10)
    plt.xlabel("Brennkammer-Eintrittstemperatur $T_\mathrm{BK}$ [K]")
    plt.title("Leistungsbedarf $P$ [kW]")
    
    
    
    
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    plt.xlim(tbk_lims)
    
    # Plot circles with letters
    for (cx, cy), letter in zip(circle_positions, letters):
        circle = plt.Circle((cx, cy), 6, fc='white', ec='black', lw=0.5, zorder=10)
        ax.add_patch(circle)
        ax.text(cx, cy, letter, ha='center', va='center', fontsize=10, fontweight='bold', zorder=15)
    #cb.set_ticks(clevs)
    #cb.set_ticklabels(clevs)
    plt.ticklabel_format(style="plain")
    fig.savefig(os.path.join(save_dir, name+'powercontour.pdf'), dpi=600, bbox_inches="tight")
    plt.show()






###############################################################################
######################### Wärmebedarf Contour plot ############################
###############################################################################

# fig, ax = plt.subplots(3,1, sharex=True)
# fig.subplots_adjust(wspace=0.1)
# norm = mpl.colors.Normalize(120, 700)
# for arch, plot, name in zip(subfolders, [0,1,2], systemnames):
    
#     t_wu_u = np.unique(t_wu)
#     t_bk_u = np.unique(t_bk)
#     Q3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     i=0
    
#     for t_bki in t_bk_u:
#         data_i = data[arch]
#         P_mfp = np.array(data_i["P_mfp"])
#         P_fp = P_mfp
#         Q = np.array(data_i["Q"])
    
#         t_wu = np.array(data_i["t_wu"])
#         t_bk = np.array(data_i["t_bk"])
#         qm_r = np.array(data_i["qm_r"])
#         Q = np.array(data_i["Q"])
        
#         j=0
#         idx = t_bk == t_bki
#         t_wu2 = t_wu[idx]
#         Q2 = Q[idx]
#         qm_r2 = qm_r[idx]
#         for t_wui in t_wu_u:
#             idx2 = t_wu2 == t_wui
#             try:
#                 Q3[j, i] = np.extract(True, Q2[idx2])
#                 qm_r3[j, i] = qm_r2[idx2][0]
#             except:
#                 qm_r3[j, i] = np.nan
#                 Q3[j, i] = np.nan
#             j+=1
#         i+=1
#     lev = [*np.linspace(150, 700, 600)]
#     cs = ax[plot].contourf(t_bk_u, t_wu_u, Q3/1e3, levels = lev, cmap="Reds", norm=norm)
#     ax[plot].set_title(name)

#     levs = [200, 300, 400, 500, 600]

    
#     cs = ax[plot].contour(t_bk_u, t_wu_u, Q3/1e3, levels = levs, colors='black')
#     plt.clabel(cs, fontsize=10)

# cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.1)
# cb.set_label("$Q_{gsmt}$ [kW]", rotation=0, labelpad=-50, loc="left")
# for sp in ax[0:1]:
#     sp.tick_params(labelbottom=False)
# clevs = [120, 200, 300, 400, 500, 600, 700]
# ax[1].set_ylabel("$T_{W}$ [K]")
# ax[2].set_xlabel("$T_{BK}$ [K]")

# fig = mpl.pyplot.gcf()
# fig.set_size_inches(16/2.54, 26/2.54)
# plt.xlim(tbk_lims)
# cb.set_ticks(clevs)
# cb.set_ticklabels(clevs)
# plt.ticklabel_format(style="plain")
# fig.savefig(os.path.join(save_dir, 'heatcontour.png'), dpi=600, bbox_inches="tight")
# plt.show()

###############################################################################
################### Massenstromverhältnis Pumpe ###############################
###############################################################################


# lev = [*np.linspace(0.1, 10, 800)]
# cs = plt.contourf(t_bk_u, t_wu_u, qm_r3/0.1, levels=lev, cmap="Reds", norm=mpl.colors.LogNorm(0.08, 10))
# cb = plt.colorbar(label="$\dot{m}_r/\dot{m}_{BK}$ [-]")
# clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3, 5, 10]
# cb.set_ticks(clevs)
# cb.set_ticklabels(clevs)
# plt.xlabel("$T_{BK}$ [K]")
# plt.ylabel("$T_{W}$ [K]")
# clevs = [0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2, 3]
# cs = plt.contour(t_bk_u, t_wu_u, qm_r3/0.1, levels = clevs, colors='black')
# plt.clabel(cs, fontsize=8)
# plt.title("Architektur: " + arch)
# fig = mpl.pyplot.gcf()
# fig.set_size_inches(16/2.54, 10.5/2.54)

# plt.xlim(tbk_lims)
# fig.savefig(os.path.join(save_dir, 'pump_massflow.png'), dpi=600)
# plt.show()


###############################################################################
######################### Wärmebedarf Contour plot ############################
###############################################################################


norm = mpl.colors.Normalize(120, 900)
for arch, name in zip(subfolders, systemnames):
    
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
    lev = [*np.linspace(150, 900, 600)]
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.contourf(t_bk_u, t_wu_u, Q3/1e3, levels = lev, cmap="Reds", norm=norm, zorder=0)

    levs = [149, 200, 300, 400, 500, 600, 700, 800]

    
    cs = ax.contour(t_bk_u, t_wu_u, Q3/1e3, levels = levs, colors='black', zorder=5)
    plt.clabel(cs, fontsize=10)

    # cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.1)
    # cb.set_label("$Q_{gsmt}$ [kW]", rotation=0, labelpad=-50, loc="left")
    # for sp in ax[0:1]:
    #     sp.tick_params(labelbottom=False)
    # clevs = [120, 200, 300, 400, 500, 600, 700]
    plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]", fontsize=10)
    plt.xlabel("Brennkammer-Eintrittstemperatur $T_\mathrm{BK}$ [K]")
    plt.title("Wärmebedarf $\dot{Q}$ [kW]")
    
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    plt.xlim(tbk_lims)
    # Plot circles with letters
    for (cx, cy), letter in zip(circle_positions, letters):
        circle = plt.Circle((cx, cy), 6, fc='white', ec='black', lw=0.5, zorder=10)
        ax.add_patch(circle)
        ax.text(cx, cy, letter, ha='center', va='center', fontsize=10, fontweight='bold', zorder=15)
    # cb.set_ticks(clevs)
    # cb.set_ticklabels(clevs)
    plt.ticklabel_format(style="plain")
    fig.savefig(os.path.join(save_dir, name+'heatcontour.pdf'), dpi=600, bbox_inches="tight")
    plt.show()
    
    
###############################################################################
##################### Wasserstoffbedarf Contour plot ##########################
###############################################################################


norm = mpl.colors.Normalize(0, 6)
for arch, name in zip(subfolders, systemnames):
    
    t_wu_u = np.unique(t_wu)
    t_bk_u = np.unique(t_bk)
    Q3 = np.zeros([len(t_wu_u), len(t_bk_u)])
    qm_phc3 = np.zeros([len(t_wu_u), len(t_bk_u)])
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
        qm_phc = np.array(data_i["qm_phc"])
        
        j=0
        idx = t_bk == t_bki
        t_wu2 = t_wu[idx]
        Q2 = Q[idx]
        qm_phc2 = qm_phc[idx]
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                Q3[j, i] = np.extract(True, Q2[idx2])
                qm_phc3[j, i] = qm_phc2[idx2][0]
            except:
                qm_phc3[j, i] = np.nan
                Q3[j, i] = np.nan
            j+=1
        i+=1
    lev = [*np.linspace(0, 6, 600)]
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    cs = plt.contourf(t_bk_u, t_wu_u, qm_phc3*1e3, levels = lev, cmap="Reds", norm=norm, zorder = 0)

    levs = [*np.linspace(0, 6, 7)]

    
    cs = plt.contour(t_bk_u, t_wu_u, qm_phc3*1e3, levels = levs, colors='black', zorder = 5)
    plt.clabel(cs, fontsize=10)

    # cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.1)
    # cb.set_label("$Q_{gsmt}$ [kW]", rotation=0, labelpad=-50, loc="left")
    # for sp in ax[0:1]:
    #     sp.tick_params(labelbottom=False)
    # clevs = [120, 200, 300, 400, 500, 600, 700]
    plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]", fontsize=10)
    plt.xlabel("Brennkammer-Eintrittstemperatur $T_\mathrm{BK}$ [K]")
    plt.title("Wasserstoffbedarf $\dot{m}_\mathrm{PHC}$ [g/s]")
    
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 10.5/2.54)
    plt.xlim(tbk_lims)
    # Plot circles with letters
    for (cx, cy), letter in zip(circle_positions, letters):
        circle = plt.Circle((cx, cy), 6, fc='white', ec='black', lw=0.5, zorder=10)
        ax.add_patch(circle)
        ax.text(cx, cy, letter, ha='center', va='center', fontsize=10, fontweight='bold', zorder=15)
    # cb.set_ticks(clevs)
    # cb.set_ticklabels(clevs)
    plt.ticklabel_format(style="plain")
    fig.savefig(os.path.join(save_dir, name+'heatcontour.pdf'), dpi=600, bbox_inches="tight")
    plt.show()
    

    
###############################################################################
################# Leistungsbedarf Pumpe Contour plot ##########################
###############################################################################


# norm = mpl.colors.Normalize(0, 0.6)
# for arch, plot, name in zip(subfolders, [0,1,2], systemnames):
    
#     t_wu_u = np.unique(t_wu)
#     t_bk_u = np.unique(t_bk)
#     P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     i=0
    
#     for t_bki in t_bk_u:
#         data_i = data[arch]
#         P_mfp = np.array(data_i["P_mfp"])
#         P_fp = P_mfp
    
#         t_wu = np.array(data_i["t_wu"])
#         t_bk = np.array(data_i["t_bk"])
#         qm_r = np.array(data_i["qm_r"])
#         Q = np.array(data_i["Q"])
#         if "P_r" in data_i:
#             if np.array(data_i["P_r"]).size == P_mfp.size:
#                 P_r = np.array(data_i["P_r"])
#                 P_mfp += np.array(data_i["P_r"])
#         P_mfp=P_mfp/(P_mfp+Q)
#         j=0
#         idx = t_bk == t_bki
#         t_wu2 = t_wu[idx]
#         P_mfp2 = P_mfp[idx]
#         qm_r2 = qm_r[idx]
#         for t_wui in t_wu_u:
#             idx2 = t_wu2 == t_wui
#             try:
#                 P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
#                 qm_r3[j, i] = qm_r2[idx2][0]
#             except:
#                 qm_r3[j, i] = np.nan
#                 P_mfp3[j, i] = np.nan
#             j+=1
#         i+=1
#     lev = [*np.linspace(0, 1, 400)]
#     cs = plt.contourf(t_bk_u, t_wu_u, P_mfp3, levels = lev, cmap="Reds", norm=norm)
#     levs = np.linspace(0, 1, 11)
    
#     cs = plt.contour(t_bk_u, t_wu_u, P_mfp3, levels = levs, colors='black')
#     plt.clabel(cs, fontsize=10)

#     #cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.08)
#     #cb.set_label("$P_{mech}$ [kW]", rotation=0, labelpad=-50, loc="left")
#     #clevs = [30, 40, 50, 60, 80, 100, 120, 160, 200, 250]
#     plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_{W}$ [K]", fontsize=12)
#     plt.xlabel("Brennkammer-Eintrittstemperatur $T_{BK}$ [K]")
    
#     fig = mpl.pyplot.gcf()
#     fig.set_size_inches(16/2.54, 10.5/2.54)
#     plt.xlim(tbk_lims)
#     #cb.set_ticks(clevs)
#     #cb.set_ticklabels(clevs)
#     plt.ticklabel_format(style="plain")
#     fig.savefig(os.path.join(save_dir, name+'powerratiocontour.pdf'), dpi=600, bbox_inches="tight")
#     plt.show()
    
###############################################################################
################# Kraftstoffverbrauch Contour plot ##########################
###############################################################################


# norm = mpl.colors.Normalize(0.1113, 0.1127)
# for arch, plot, name in zip(subfolders, [0,1,2], systemnames):
    
#     t_wu_u = np.unique(t_wu)
#     t_bk_u = np.unique(t_bk)
#     P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     qm3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     i=0
    
#     for t_bki in t_bk_u:
#         data_i = data[arch]
#         P_mfp = np.array(data_i["P_mfp"])
#         P_fp = P_mfp
    
#         t_wu = np.array(data_i["t_wu"])
#         t_bk = np.array(data_i["t_bk"])
#         qm_r = np.array(data_i["qm_r"])
#         qm = np.array(data_i["qm_cb"])+np.array(data_i["qm_phc"])+np.array(data_i["qm_pot"])
#         Q = np.array(data_i["Q"])
#         if "P_r" in data_i:
#             if np.array(data_i["P_r"]).size == P_mfp.size:
#                 P_r = np.array(data_i["P_r"])
#                 P_mfp += np.array(data_i["P_r"])
#         P_mfp=P_mfp/(P_mfp+Q)
#         j=0
#         idx = t_bk == t_bki
#         qm2 = qm[idx]
#         t_wu2 = t_wu[idx]
#         P_mfp2 = P_mfp[idx]
#         qm_r2 = qm_r[idx]
#         for t_wui in t_wu_u:
#             idx2 = t_wu2 == t_wui
#             try:
#                 P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
#                 qm_r3[j, i] = qm_r2[idx2][0]
#                 qm3[j,i] = qm2[idx2]
#             except:
#                 qm_r3[j, i] = np.nan
#                 P_mfp3[j, i] = np.nan
#                 qm3[j,i] = np.nan
#             j+=1
#         i+=1
#     lev = [*np.linspace(0.1113, 0.1127, 400)]
#     cs = plt.contourf(t_bk_u, t_wu_u, qm3, levels = lev, cmap="Reds", norm=norm)
#     levs = np.linspace(0.1113, 0.1127, 11)
    
#     cs = plt.contour(t_bk_u, t_wu_u, qm3, levels = levs, colors='black')
#     plt.clabel(cs, fontsize=10)

#     #cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Reds"), orientation="horizontal", location="bottom", ax=ax, pad=0.08)
#     #cb.set_label("$P_{mech}$ [kW]", rotation=0, labelpad=-50, loc="left")
#     #clevs = [30, 40, 50, 60, 80, 100, 120, 160, 200, 250]
#     plt.ylabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]", fontsize=12)
#     plt.xlabel("Brennkammer-Eintrittstemperatur $T_\mathrm{BK}$ [K]")
    
#     fig = mpl.pyplot.gcf()
#     fig.set_size_inches(16/2.54, 10.5/2.54)
#     plt.xlim(tbk_lims)
#     #cb.set_ticks(clevs)
#     #cb.set_ticklabels(clevs)
#     plt.ticklabel_format(style="plain")
#     fig.savefig(os.path.join(save_dir, name+'massflowcontour.pdf'), dpi=600, bbox_inches="tight")
#     plt.show()



