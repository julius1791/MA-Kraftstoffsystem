# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
import h2_properties as h2
import h2flow
csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

save_dir = os.path.join(os.getcwd(), "diagrams")
qm=0.100412959+0.042830625
t1 = 200
h1 = h2.calc_H2_enthalpy(t1, 1756736.842)
H1 = h1*(qm)/1e3
H2 = H1 + 200
t2 = h2flow.calc_tp(H2*1e3/qm, 1756736.842, 0, True)

H3 = H2 + 824765.2865/1e3-200
t3 = h2flow.calc_tp(H3*1e3/qm, 1668900, 0, True)


fig, ax = plt.subplots()
ax.plot([0, H2-H1, H3-H1], [t1, t2, t3], label=" ", color="black", linestyle="-")
#plt.plot([],[], label="Wärmequellen:", color="none")
ax.plot([0, H2-H1], [370, 400], label="Ölsystem", color="black", linestyle=":")
ax.plot([],[], label=" ", color="none")
ax.plot([H2-H1, H3-H1], [600, 1600], label="Parallele Verbennung", color="black", linestyle="--")
h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*2
handles = ph + h
labels = ["Wasserstoff:", "Wärmequellen:"] + l

plt.ylim([0, 1200])
plt.xlim([0, H3-H1])
plt.xlabel("Leistung $\dot{H}$ [kW]")
plt.ylabel("Temperatur $T$ [K]")
ax.annotate('', xy=(0, 80), xytext=(H2-H1, 80), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|'})
ax.annotate('$\dot{Q}_{FOHE}$', xy=((H2-H1)/2, 135), ha='center', va='center')
ax.annotate('', xy=(H2-H1, 80), xytext=(H3-H1, 80), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|'})
ax.annotate('$\dot{Q}_{PHC}$', xy=((H2+H3-2*H1)/2, 135), ha='center', va='center')

ax.arrow(H2-H1, t2, 0, 400-t2, length_includes_head=True, width=2, head_width=10, head_length=20, linestyle="--", fill=True, ec="none", fc="black")
ax.arrow(H2-H1, 400, 0, t2-400, length_includes_head=True, width=2, head_width=10, head_length=20, linestyle="--", fill=True, ec="none", fc="black")
ax.annotate("$\Delta T_{min}$", xytext=[H2-H1+10, 380], xy=[H2-H1, t2])
fig.set_size_inches(16/2.54, 10.5/2.54)
ax.legend(handles, labels, ncol=3, loc="lower center", bbox_to_anchor=(0.43, 1.05))
plt.tight_layout()

    

# t1, t2, t3, t4 = leg.get_texts()
# t1.set_size(9)
# t3.set_size(9)
# t4.set_size(9)
# t2._fontproperties = t3._fontproperties.copy()
# t2.set_size(12)







fig.savefig(os.path.join(save_dir, 'hx.png'), dpi=600, bbox_inches="tight")
plt.show()