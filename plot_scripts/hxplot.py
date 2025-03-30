# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os
cwd = os.getcwd()
os.chdir(os.path.dirname(cwd))
import h2flow

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

save_dir = os.path.join(os.getcwd(), "diagrams")
qm=0.12
t1 = 140
H1 = h2flow.h2.calc_H2_enthalpy(t1, 1889258.86)*qm/1e3
t2 = 200
H2 = h2flow.h2.calc_H2_enthalpy(t2, 1889258.86)*qm/1e3
t3 = 330
H3 = h2flow.h2.calc_H2_enthalpy(t3, 1889258.86)*qm/1e3
t4 = 400
H4 = h2flow.h2.calc_H2_enthalpy(t4, 1889258.86)*qm/1e3


fig, ax = plt.subplots()
ax.plot([0, H2-H1, H3-H1, H4-H1], [t1, t2, t3, t4], label=" ", color="black", linestyle="-")
ax.plot([0, H2-H1], [297, 336], label="Klimasystem", color="black", linestyle=":")
ax.plot([],[], label=" ", color="none")
ax.plot([H2-H1, H3-H1], [400, 440], label="Ölsystem", color="black", linestyle="-.")
ax.plot([],[], label=" ", color="none")
ax.plot([H3-H1, H4-H1], [370, 800], label="Parallele Verbennung", color="black", linestyle="--")
h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*2
handles = ph + h
labels = ["Wasserstoff", "Wärmequellen:"] + l


plt.ylim([0, 600])
plt.xlim([0, H3-H1+100])
plt.xlabel("Wärme $\dot{Q}$")
plt.ylabel("Temperatur $T$")
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
plt.yticks([])


ax.annotate('', xy=(0, 70), xytext=(H2-H1, 70), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|'})
ax.annotate('$\dot{Q}_{\mathrm{Klima}}$', xy=((H2-H1)/2, 40), ha='center', va='center')
ax.annotate('', xy=(H2-H1, 70), xytext=(H3-H1, 70), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|'})
ax.annotate('$\dot{Q}_{\mathrm{FOHE}}$', xy=((H2+H3-2*H1)/2, 40), ha='center', va='center')
ax.annotate('', xy=(H3-H1, 70), xytext=(H3-H1+100, 70), xycoords='data', textcoords='data',
            arrowprops={'arrowstyle': '|-|'})
ax.annotate('$\dot{Q}_{\mathrm{PHC}}$', xy=(H3-H1+50, 40), ha='center', va='center')

#ax.arrow(H2-H1, t2, 0, 400-t2, length_includes_head=True, width=2, head_width=10, head_length=20, linestyle="--", fill=True, ec="none", fc="black")
#ax.arrow(H2-H1, 400, 0, t2-400, length_includes_head=True, width=2, head_width=10, head_length=20, linestyle="--", fill=True, ec="none", fc="black")
#ax.annotate("$\Delta T_{min}$", xytext=[H2-H1+10, 380], xy=[H2-H1, t2])

ax.legend(handles, labels, ncols=4, loc="lower center", bbox_to_anchor=(0.5, 1.02), columnspacing=0.2)
#plt.tight_layout()

    

# t1, t2, t3, t4 = leg.get_texts()
# t1.set_size(9)
# t3.set_size(9)
# t4.set_size(9)
# t2._fontproperties = t3._fontproperties.copy()
# t2.set_size(12)






fig.set_size_inches(16/2.54, 8.5/2.54)
fig.savefig(os.path.join(save_dir, 'hx.pdf'), dpi=600, bbox_inches="tight")
plt.show()