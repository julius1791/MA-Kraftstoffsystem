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

cwd = os.getcwd()
os.chdir(os.path.dirname(cwd))

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

tbk_lims = [150,500]
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
colors = ["orangered", "navy", "green", "cyan"]

###############################################################################
######################## Leistungsbedarf Systemvergleich ######################
###############################################################################

# t_list = [100, 160, 220, 280]
# for t in t_list:
#     for key, style, name, color in zip(data, stylelist, systemnames, colors):
#         data_i = data[key]
#         P_mfp = np.array(data_i["P_mfp"])
#         t_wu = np.array(data_i["t_wu"])
#         t_bk = np.array(data_i["t_bk"])
#         if "P_r" in data_i:
#             if np.array(data_i["P_r"]).size == P_mfp.size:
#                 P_mfp += np.array(data_i["P_r"])
#         idx = t_wu == t
#         t_bk1 = t_bk[idx]
#         P_mfp1 = P_mfp[idx]
#         t_wu1 = t_wu[idx]
#         plt.plot(t_bk1, P_mfp1/1e3, label=name, color=color, linestyle=style)
#     plt.legend(title="Kraftstoffsystem", loc="best")
#     plt.xlabel("$T_{BK}$ [K]")
#     plt.ylabel("$P_{mech}$ [kW]")
#     plt.xlim([max(140,t+20), tbk_lims[1]])
#     plt.ylim([0,300])
#     plt.title("$T_{W}$ = " + str(t) + " K")
#     fig = mpl.pyplot.gcf()
#     fig.set_size_inches(16/2.54, 10.5/2.54)
#     fig.savefig(os.path.join(save_dir, 'arch_' + str(t) + '.png'), dpi=600, bbox_inches="tight")
#     plt.show()
    

###############################################################################
######################## Power split - all - tbk ######################
###############################################################################

t = 160
for key, style, name, color in zip(data, stylelist, systemnames, colors):
    data_i = data[key]
    P_mfp = np.array(data_i["P_mfp"])
    t_wu = np.array(data_i["t_wu"])
    t_bk = np.array(data_i["t_bk"])
    if "P_r" in data_i:
        if np.array(data_i["P_r"]).size == P_mfp.size:
            P_r = np.array(data_i["P_r"])
    idx = t_wu == t
    t_bk1 = t_bk[idx]
    P_mfp1 = P_mfp[idx]
    P_r1= P_r[idx]
    t_wu1 = t_wu[idx]
    if key == "pump":
        plt.plot(t_bk1, P_mfp1/1e3, label=name+" HPFP", color=color, linestyle=stylelist[0])
    else:
        plt.plot(t_bk1, P_mfp1/1e3, label=name+" HPFC", color=color, linestyle=stylelist[0])
    plt.plot(t_bk1, P_r1/1e3, label=name+" RV", color=color, linestyle=stylelist[1])
plt.legend(loc="best")
plt.xlabel("Brennkammer-Eintrittstemperatur $T_{\mathrm{BK}}$ [K]")
plt.ylabel("Leistung $P_{i}$ [kW]")
plt.title("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}= 160$ K")
plt.xlim([max(140,t+20), tbk_lims[1]])
plt.ylim([0,120])
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 10.5/2.54)
fig.savefig(os.path.join(save_dir, 'tbkcomp.pdf'), dpi=600, bbox_inches="tight")
plt.show()
    


###############################################################################
######################## Powersplit - single - TW ######################
###############################################################################
# for tbk in [300, 400]:
#     fig, ax = plt.subplots()
#     for key, name, color in zip(data, systemnames, colors):
#         data_i = data[key]
#         P_mfp = np.array(data_i["P_mfp"])
#         t_wu = np.array(data_i["t_wu"])
#         t_bk = np.array(data_i["t_bk"])
#         if "P_r" in data_i:
#             if np.array(data_i["P_r"]).size == P_mfp.size:
#                 P_r = np.array(data_i["P_r"])
#         idx = t_bk == tbk
#         t_wu1 = t_wu[idx]
#         P_mfp1 = P_mfp[idx]
#         t_bk1 = t_bk[idx]
#         P_r1 = P_r[idx]
#         if name == "Pumpe":
#             plt.plot(t_wu1, P_mfp1/1e3, label="HPFP", color=colors[0], linestyle=stylelist[0])
#         else:
#             plt.plot(t_wu1, P_mfp1/1e3, label="HPFC", color=colors[0], linestyle=stylelist[0])
#         plt.plot(t_wu1, P_r1/1e3, label="RV", color=colors[1], linestyle=stylelist[0])
        
#         plt.legend(loc="best")
#         plt.xlabel("Wärmeübertrager-Eintrittstemperatur $T_{W}$ [K]")
#         plt.ylabel("Leistung $P_{i}$ [kW]")
#         plt.xlim(twu_lims)
#         plt.ylim([0,200])
#         fig = mpl.pyplot.gcf()
#         fig.set_size_inches(12/2.54, 8.5/2.54)
#         fig.savefig(os.path.join(save_dir, str(tbk)+name+'_powersplit.png'), dpi=600, bbox_inches="tight")
#         plt.show()


###############################################################################
######################## Powersplit - single - 2TW ######################
###############################################################################


for key, name in zip(data, systemnames):
    fig, ax = plt.subplots(1, 2, sharey=True)
    for tbk, i in zip([300, 400], [0, 1]):
        data_i = data[key]
        P_mfp = np.array(data_i["P_mfp"])
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        if "P_r" in data_i:
            if np.array(data_i["P_r"]).size == P_mfp.size:
                P_r = np.array(data_i["P_r"])
        idx = t_bk == tbk
        t_wu1 = t_wu[idx]
        P_mfp1 = P_mfp[idx]
        t_bk1 = t_bk[idx]
        P_r1 = P_r[idx]
        if name == "Pumpe":
            ax[i].plot(t_wu1, P_mfp1/1e3, label="HPFP", color=colors[0], linestyle=stylelist[0])
        else:
            ax[i].plot(t_wu1, P_mfp1/1e3, label="HPFC", color=colors[0], linestyle=stylelist[0])
        ax[i].plot(t_wu1, P_r1/1e3, label="RV", color=colors[1], linestyle=stylelist[0])
        
        plt.legend(loc="best")
        ax[i].set_title("$T_\mathrm{BK}="+str(tbk)+"$ K")
        
        ax[i].set_xlim(twu_lims)
        ax[i].set_ylim([0,200])
        
    fig.supxlabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]", y=0.1, fontsize=12)
    ax[0].set_ylabel("Leistung $P_{i}$ [kW]")
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 8.5/2.54)
    plt.tight_layout()
    fig.savefig(os.path.join(save_dir, name+'_powersplit.pdf'), dpi=600, bbox_inches="tight")
    plt.show()

###############################################################################
######################## Powersplit - all - TW ######################
###############################################################################
for tbk in [300, 400]:
    fig, ax = plt.subplots()
    for key, name, color in zip(data, systemnames, colors):
        data_i = data[key]
        P_mfp = np.array(data_i["P_mfp"])
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        if "P_r" in data_i:
            if np.array(data_i["P_r"]).size == P_mfp.size:
                P_r = np.array(data_i["P_r"])
        idx = t_bk == tbk
        t_wu1 = t_wu[idx]
        P_mfp1 = P_mfp[idx]
        t_bk1 = t_bk[idx]
        P_r1 = P_r[idx]
        if name == "Pumpe":
            plt.plot(t_wu1, P_mfp1/1e3, label=name+" HPFP", color=color, linestyle=stylelist[0])
        else:
            plt.plot(t_wu1, P_mfp1/1e3, label=name+" HPFC", color=color, linestyle=stylelist[0])
        plt.plot(t_wu1, P_r1/1e3, label=name+" RV", color=color, linestyle=stylelist[1])
        
    plt.legend(loc="best")
    plt.xlabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]")
    plt.ylabel("Leistung $P_{i}$ [kW]")
    plt.title("Brennkammer-Eintrittstemperatur $T_\mathrm{BK}="+str(tbk)+"$ K")
    plt.xlim(twu_lims)
    plt.ylim([0,200])
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(16/2.54, 8.5/2.54)
    fig.savefig(os.path.join(save_dir, str(tbk)+'summary_powersplit.pdf'), dpi=600, bbox_inches="tight")
    plt.show()
    
###############################################################################
######################## Leistung gesamt - all - TW ######################
###############################################################################

fig, axs = plt.subplots(1,2, sharey=True)
fig.subplots_adjust(wspace=0.12)
for tbk, j in zip([300, 400], [0,1]):
    for key, name, i in zip(data, systemnames, [0,1,2]):
        data_i = data[key]
        P_mfp = np.array(data_i["P_mfp"])
        t_wu = np.array(data_i["t_wu"])
        t_bk = np.array(data_i["t_bk"])
        Q = np.array(data_i["Q"])
        if "P_r" in data_i:
            if np.array(data_i["P_r"]).size == P_mfp.size:
                P_r = np.array(data_i["P_r"])
        idx = t_bk == tbk
        t_wu1 = t_wu[idx]
        P_mfp1 = P_mfp[idx]
        t_bk1 = t_bk[idx]
        P_r1 = P_r[idx]
        Q1 = Q[idx]
        
        if name == "Pumpe":
            axs[j].plot(t_wu1, (P_mfp1+P_r1)/1e3, label=name, color=colors[i], linestyle=stylelist[i])
        else:
            axs[j].plot(t_wu1, (P_mfp1+P_r1)/1e3, label=name, color=colors[i], linestyle=stylelist[i])
        #plt.plot(t_wu1, Q1/1e3, label="$Q$", color=colors[i], linestyle=stylelist[1])
        axs[j].set_title("$T_\mathrm{BK}=" + str(tbk)+ "$ K")
        plt.ylim([0,200])
handles, labels = plt.gca().get_legend_handles_labels()
axs[1].legend(handles, labels, loc="upper right")
fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 7.5/2.54)
#fig.subplots_adjust(wspace=0.1, hspace=0.1)
axs[0].set_ylabel("Leistung $P_\mathrm{gsmt}$ [kW]")
fig.add_subplot(111, frameon=False)
axs[0].set_xlim([100,250])
axs[1].set_xlim([100,250])
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.xlabel("Wärmeübertrager-Eintrittstemperatur $T_\mathrm{W}$ [K]")
#axs[0].legend(ncols=3, bbox_to_anchor=(0.5,1.2), loc="upper center")
fig.savefig(os.path.join(save_dir, 'summary_power.pdf'), dpi=600, bbox_inches="tight")
plt.show()

###############################################################################
################# Leistungsbedarf Pumpe T_W Vergleich #########################
###############################################################################


# data_i = data[subfolders[2]]
# P_mfp = np.array(data_i["P_mfp"])
# P_fp = P_mfp

# t_wu = np.array(data_i["t_wu"])
# t_bk = np.array(data_i["t_bk"])
# qm_r = np.array(data_i["qm_r"])
# Q = np.array(data_i["Q"])
# if "P_r" in data_i:
#     if np.array(data_i["P_r"]).size == P_mfp.size:
#         P_r = np.array(data_i["P_r"])
#         P_mfp += np.array(data_i["P_r"])
# t_hx_list = reversed([100, 160, 220, 280])
# for t_hx, style, color in zip(t_hx_list, stylelist, colors):
    
#     idx = t_wu == t_hx
#     t_bk1 = t_bk[idx]
#     P_mfp1 = P_mfp[idx]
#     t_wu1 = t_wu[idx]
#     plt.plot(t_bk1, P_mfp1/1e3, label=str(t_hx) + " K", linestyle=style, color=color)
# plt.legend(title="$T_{W}$ [K]")
# plt.xlabel("$T_{BK}$ [K]")
# plt.ylabel("$P_{mech}$ [kW]")
# plt.title("Architektur: Pumpe")
# plt.xlim(tbk_lims)
# plt.ylim([0,300])
# fig = mpl.pyplot.gcf()
# fig.set_size_inches(16/2.54, 10.5/2.54)
# fig.savefig(os.path.join(save_dir, 'pump.png'), dpi=600)
# plt.show()


###############################################################################
##################### Leistungsbedarf T_W Vergleich ###########################
###############################################################################

# fig, ax = plt.subplots()
# for arch, color in zip(subfolders, colors):
#     data_i = data[arch]
#     P_mfp = np.array(data_i["P_mfp"])
#     P_fp = P_mfp
    
#     t_wu = np.array(data_i["t_wu"])
#     t_bk = np.array(data_i["t_bk"])
#     qm_r = np.array(data_i["qm_r"])
#     Q = np.array(data_i["Q"])
#     if "P_r" in data_i:
#         if np.array(data_i["P_r"]).size == P_mfp.size:
#             P_r = np.array(data_i["P_r"])
#             P_mfp += np.array(data_i["P_r"])
#     t_hx_list = [280, 220, 160, 100]
#     for t_hx, style in zip(t_hx_list, stylelist):
        
#         idx = t_wu == t_hx
#         t_bk1 = t_bk[idx]
#         P_mfp1 = P_mfp[idx]
#         t_wu1 = t_wu[idx]
#         plt.plot(t_bk1, P_mfp1/1e3, linestyle=style, color=color)
# tit = Rectangle([0,0], 1, 0.1, fc="none", ec="none")
# lvm = Rectangle([0,0], 1, 1, fc=colors[0])
# lvd = Rectangle([0,0], 1, 1, fc=colors[1])
# lp = Rectangle([0,0], 1, 1, fc=colors[2])
# l100 = Line2D([0], [0], color="black", linestyle=stylelist[0])
# l160 = Line2D([0], [0], color="black", linestyle=stylelist[1])
# l220 = Line2D([0], [0], color="black", linestyle=stylelist[2])
# l280 = Line2D([0], [0], color="black", linestyle=stylelist[3])
# leg_handles = [tit, tit, l100, lvm, l160, lvd, l220, lp, l280]
# leg_labels = ["$T_{W}$ [K]", "Architektur", t_hx_list[0], subfolders[0], t_hx_list[1], subfolders[1], t_hx_list[2], subfolders[2], t_hx_list[3]]
# plt.xlabel("$T_{BK}$ [K]")
# plt.ylabel("$P_{mech}$ [kW]")
# plt.xlim(tbk_lims)
# plt.ylim([0,300])
# fig = mpl.pyplot.gcf()
# fig.set_size_inches(16/2.54, 12/2.54)
# plt.subplots_adjust(top=0.7)
# plt.legend(ncols=5, handles=leg_handles, labels=leg_labels, bbox_to_anchor=(0.5,1.15), loc="upper center")
# fig.savefig(os.path.join(save_dir, 'pmech.png'), dpi=600)
# plt.show()


    

###############################################################################
########################### Energie Stackplots ################################
###############################################################################


# for subfolder,name in zip(subfolders, systemnames):
#     t_wu_u = np.unique(t_wu)
#     t_bk_u = np.unique(t_bk)
#     P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     qm_r3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     i=0

#     for t_bki in t_bk_u:
#         j=0
#         idx = t_bk == t_bki
#         t_wu2 = t_wu[idx]
#         P_mfp2 = P_mfp[idx]
#         for t_wui in t_wu_u:
#             idx2 = t_wu2 == t_wui
#             try:
#                 P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
    
#             except:
#                 P_mfp3[j, i] = np.nan
#             j+=1
#         i+=1
#     data_i = data[subfolder]
#     P_mfp = np.array(data_i["P_mfp"])
#     P_fp = P_mfp
    
#     t_wu = np.array(data_i["t_wu"])
#     t_bk = np.array(data_i["t_bk"])
#     qm_r = np.array(data_i["qm_r"])
#     qm_cb = np.array(data_i["qm_cb"])
#     qm_phc = np.array(data_i["qm_phc"])
#     Q = np.array(data_i["Q"])
#     if "P_r" in data_i:
#         if np.array(data_i["P_r"]).size == P_mfp.size:
#             P_r = np.array(data_i["P_r"])
#             P_mfp += np.array(data_i["P_r"])
            
#     id_200 = t_wu == 120
#     t_bk_200 = t_bk[id_200]
#     P_fp_200 = P_fp[id_200]/1e3
#     if subfolder != "pre":
#         P_r = P_r[id_200]/1e3
#         P_fp_200 -= P_r
#     else:
#         P_r = P_fp_200 * 0
#     Q_200 = Q[id_200]/1e3
#     qm_200 = qm_cb[id_200] + qm_phc[id_200]
#     Q_hx = Q_200.clip(0, 159)
#     Q_phc = Q_200 - Q_hx
#     Q_excess = Q_200 - 159
#     Q_excess = Q_excess.clip(-159, 0)
#     Q_hx2 = [159 for i in Q_200]
    
    
#     spQ = plt.stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
#     #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
#     spP = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
#     p1 = Rectangle((0, 0), 1, 1, fc="navy")
#     p2 = Rectangle((0, 0), 1, 1, fc="cyan")
#     p3 = Rectangle((0, 0), 1, 1, fc="darkred")
#     p4 = Rectangle((0, 0), 1, 1, fc="orangered")
#     if subfolder != "pump":
#         plt.legend([p3, p4, p2, p1], ["$P_{HPFC}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)
#     else:
#         plt.legend([p3, p4, p2, p1], ["$P_{HPFP}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)

#     plt.xlabel("Brennkammer-Eintrittstemperatur $T_{BK}$ [K]")
#     plt.ylabel("Leistung [kW]")
#     plt.title("Architektur: " + name)
#     plt.xlim([150, 320])
#     plt.ylim([0, 500])
#     p = 1.7e6
#     h = list()
#     # for t in t_bk_200:
#     #     qmi = qm_200[t_bk_200 == t]
#     #     Q_hxi = Q_hx[t_bk_200 == t]/1e3
#     #     Q_phci = Q_phc[t_bk_200 == t]/1e3
#     #     P_fpi = P_fp_200[t_bk_200 == t]/1e3
#     #     P_ri = P_r[t_bk_200 == t]/1e3
#     #     hi = (h2.calc_H2_enthalpy(t, p) - h2.calc_H2_enthalpy(22, 4.2e5))/1e3 * qmi
#     #     #print((hi-Q_phci - Q_hxi-P_fpi-P_ri)/hi)
#     #     h.append(hi)
        
#     # plt.plot(t_bk_200, h)
#     fig = mpl.pyplot.gcf()
#     fig.set_size_inches(16/2.54, 10.5/2.54)
#     fig.savefig(os.path.join(save_dir, 'stackplot_' + subfolder + '.png'), dpi=600, bbox_inches="tight")
#     plt.show()
    
    
###############################################################################
########################### Energie Stackplots ################################
###############################################################################


# for subfolder,name in zip(subfolders, systemnames):
#     t_wu_u = np.unique(t_wu)
#     t_bk_u = np.unique(t_bk)
#     P_mfp3 = np.zeros([len(t_wu_u), len(t_bk_u)])
#     i=0

#     for t_bki in t_bk_u:
#         j=0
#         idx = t_bk == t_bki
#         t_wu2 = t_wu[idx]
#         P_mfp2 = P_mfp[idx]
#         for t_wui in t_wu_u:
#             idx2 = t_wu2 == t_wui
#             try:
#                 P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
#             except:
#                 P_mfp3[j, i] = np.nan
#             j+=1
#         i+=1
#     data_i = data[subfolder]
#     P_mfp = np.array(data_i["P_mfp"])
#     P_fp = P_mfp
    
#     t_wu = np.array(data_i["t_wu"])
#     t_bk = np.array(data_i["t_bk"])
#     qm_r = np.array(data_i["qm_r"])
#     qm_cb = np.array(data_i["qm_cb"])
#     qm_phc = np.array(data_i["qm_phc"])
#     Q = np.array(data_i["Q"])
#     if "P_r" in data_i:
#         if np.array(data_i["P_r"]).size == P_mfp.size:
#             P_r = np.array(data_i["P_r"])
#             P_mfp += np.array(data_i["P_r"])
            
#     id_200 = t_wu == t_bk - 40
#     t_bk_200 = t_bk[id_200]
#     P_fp_200 = P_fp[id_200]/1e3
#     if subfolder != "pre":
#         P_r = P_r[id_200]/1e3
#         P_fp_200 -= P_r
#     else:
#         P_r = P_fp_200 * 0
#     Q_200 = Q[id_200]/1e3
#     qm_200 = qm_cb[id_200] + qm_phc[id_200]
#     Q_hx = Q_200.clip(0, 159)
#     Q_phc = Q_200 - Q_hx
#     Q_excess = Q_200 - 159
#     Q_excess = Q_excess.clip(-159, 0)
#     Q_hx2 = [159 for i in Q_200]
    
    
#     spQ = plt.stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
#     #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
#     spP = plt.stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
#     p1 = Rectangle((0, 0), 1, 1, fc="navy")
#     p2 = Rectangle((0, 0), 1, 1, fc="cyan")
#     p3 = Rectangle((0, 0), 1, 1, fc="darkred")
#     p4 = Rectangle((0, 0), 1, 1, fc="orangered")
#     if subfolder != "pump":
#         plt.legend([p3, p4, p2, p1], ["$P_{HPFC}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)
#     else:
#         plt.legend([p3, p4, p2, p1], ["$P_{HPFP}$", "$P_{R}$", "$\dot{Q}_{PHC}$", "$\dot{Q}_{FOHE}$"], loc="lower right", framealpha=1)

#     plt.xlabel("Brennkammer-Eintrittstemperatur $T_{BK}$ [K]")
#     plt.ylabel("Leistung [kW]")
#     plt.title("Architektur: " + name)
#     plt.xlim([150, 320])
#     plt.ylim([0, 500])
#     p = 1.7e6
#     h = list()
#     # for t in t_bk_200:
#     #     qmi = qm_200[t_bk_200 == t]
#     #     Q_hxi = Q_hx[t_bk_200 == t]/1e3
#     #     Q_phci = Q_phc[t_bk_200 == t]/1e3
#     #     P_fpi = P_fp_200[t_bk_200 == t]/1e3
#     #     P_ri = P_r[t_bk_200 == t]/1e3
#     #     hi = (h2.calc_H2_enthalpy(t, p) - h2.calc_H2_enthalpy(22, 4.2e5))/1e3 * qmi
#     #     #print((hi-Q_phci - Q_hxi-P_fpi-P_ri)/hi)
#     #     h.append(hi)
        
#     # plt.plot(t_bk_200, h)
#     fig = mpl.pyplot.gcf()
#     fig.set_size_inches(16/2.54, 10.5/2.54)
#     fig.savefig(os.path.join(save_dir, 'stackplot_' + subfolder + '.png'), dpi=600, bbox_inches="tight")
#     plt.show()
    
###############################################################################
######################## Leistungssplit vs Rezirkulation ######################
###############################################################################




# t = 300
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# for key, name, color in zip(data, systemnames, colors):
#     data_i = data[key]
#     P_mfp = np.array(data_i["P_mfp"])
#     t_bk = np.array(data_i["t_bk"])
#     t_wu = np.array(data_i["t_wu"])
#     qm_r = np.array(data_i["qm_r"])
#     qm_cb = np.array(data_i["qm_cb"])
#     qm_r = qm_r/qm_cb
#     if "P_r" in data_i:
#         if np.array(data_i["P_r"]).size == P_mfp.size:
#             P_mfp += np.array(data_i["P_r"])
#             RV_ratio = np.array(data_i["P_r"])/P_mfp
#     if "qm_v" in data_i:
#         if np.array(data_i["qm_v"]).size == P_mfp.size:
#             qm_r += np.array(data_i["qm_v"])
#     idx = t_bk == t
#     t_wu=t_wu[idx]
#     qm_r = qm_r[idx]
#     RV_ratio = RV_ratio[idx]
#     ax1.plot(t-t_wu, RV_ratio, color=color, label="$P_{RV}/P_{mech}$ "+name, linestyle=stylelist[0])
#     ax2.plot(t-t_wu, qm_r, color=color, label="$\dot{m}_r/\dot{m}_{BK}$ "+name, linestyle=stylelist[1])
# h1, l1 = ax1.get_legend_handles_labels()
# h2, l2 = ax2.get_legend_handles_labels()
# legend = ax1.legend(h1,l1, loc="upper right", fontsize=8)
# legend2 = ax2.legend(h2, l2,loc="lower left", fontsize=8)



# ax1.set_xlabel("Temperaturdifferenz $T_{BK}-T_{W}$ [K]")
# ax1.set_ylabel("$P_{RV}/P_{mech}$ [-]")
# ax2.set_ylabel("$\dot{m}_r/\dot{m}_{BK}$")
# plt.xlim([20,200])
# ax1.set_ylim([0,1])
# ax2.set_ylim([0,5])
# fig = mpl.pyplot.gcf()
# fig.set_size_inches(16/2.54, 10.5/2.54)
# fig.savefig(os.path.join(save_dir, 'rezirkulation.png'), dpi=600, bbox_inches="tight")
# plt.show()

###############################################################################
######################## Stackplot Grid ######################
###############################################################################

csfont = {'family': "serif", "serif": ["lmr"], "size": 8}
plt.rc('font',**csfont)
systemnames = ["Vormisch.", "Verdampfer", "Pumpe"]


fig, axs = plt.subplots(3,3)


###############################################################################
# T_W = 160
p2 = 0
for subfolder,name, p1 in zip(subfolders, systemnames, [0, 1, 2]):
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
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
    
            except:
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
            
    id_200 = t_wu == 160
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
    
    
    spQ = axs[p1,p2].stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
    #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
    spP = axs[p1,p2].stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
    axs[2,0].set_xlabel("$T_\mathrm{BK}$ [K]")
    axs[1,0].set_ylabel("Leistung [kW]")
    axs[p1,p2].set_title(name+" $T_\mathrm{W}=160$ K")
    axs[p1,p2].set_xlim([180, 320])
    axs[p1,p2].set_ylim([0, 500])
    if p1 != 2:
        axs[p1,p2].set_xticklabels([])
    
###############################################################################
# T_W = T_BK - 40
p2 = 1
for subfolder,name, p1 in zip(subfolders, systemnames, [0, 1, 2]):
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
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
    
            except:
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
    
    
    spQ = axs[p1,p2].stackplot(t_bk_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
    #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
    spP = axs[p1,p2].stackplot(t_bk_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
    
    axs[2,1].set_xlabel("$T_\mathrm{BK}$ [K]")
    axs[p1,p2].set_title(name+" $T_\mathrm{W}=T_\mathrm{BK}-40$ K")
    axs[p1,p2].set_xlim([180, 320])
    axs[p1,p2].set_ylim([0, 500])
    axs[p1,p2].set_yticklabels([])
    if p1 != 2:
        axs[p1,p2].set_xticklabels([])
    
    
###############################################################################
# T_BK = 300
p2 = 2
for subfolder,name, p1 in zip(subfolders, systemnames, [0, 1, 2]):
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
        for t_wui in t_wu_u:
            idx2 = t_wu2 == t_wui
            try:
                P_mfp3[j, i] = np.extract(True, P_mfp2[idx2])
    
            except:
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
            
    id_200 = t_bk == 300
    t_wu_200 = t_wu[id_200]
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
    
    
    spQ = axs[p1,p2].stackplot(t_wu_200, Q_hx2 , Q_phc, colors=['navy', 'cyan'], baseline="zero", ec=["none"])
    #sp3 = plt.stackplot(t_bk_200, 200, colors=['white'], baseline="zero", hatch=['//'], edgecolors=["black", "black"])
    spP = axs[p1,p2].stackplot(t_wu_200, Q_hx, Q_phc, P_r, P_fp_200, colors=['none', 'none', 'orangered', 'darkred'], baseline="zero", ec=["none"])
    
    axs[2,2].set_xlabel("$T_\mathrm{W}$ [K]")
    axs[p1,p2].set_title(name+" $T_\mathrm{BK}=300$ K")
    axs[p1,p2].set_xlim([100, 280])
    axs[p1,p2].set_ylim([0, 500])
    axs[p1,p2].set_yticklabels([])
    if p1 != 2:
        axs[p1,p2].set_xticklabels([])
    
    
p1 = Rectangle((0, 0), 1, 1, fc="navy")
p2 = Rectangle((0, 0), 1, 1, fc="cyan")
p3 = Rectangle((0, 0), 1, 1, fc="darkred")
p4 = Rectangle((0, 0), 1, 1, fc="orangered")

plt.legend([p3, p4, p2, p1], ["$P_\mathrm{HPFC/P}$", "$P_\mathrm{RV}$", "$\dot{Q}_\mathrm{PHC}$", "$\dot{Q}_\mathrm{FOHE}$"], loc="lower left", framealpha=1)


fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 16/2.54)
fig.savefig(os.path.join(save_dir, 'stackplot_summary.pdf'), dpi=600, bbox_inches="tight")
plt.show()


