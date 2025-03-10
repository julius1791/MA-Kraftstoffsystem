# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import numpy as np

csfont = {'family': "serif", "serif": ["lmr"], "size": 12}
plt.rc('font',**csfont)
plt.rc('text', usetex=True)
par = {'mathtext.default': 'regular'}
plt.rcParams.update(par)

save_dir = os.path.join(os.getcwd(), "diagrams")
csv_dir = os.path.join(os.getcwd(), "material properties", "h2spin.csv")

data = np.genfromtxt(csv_dir, delimiter=",")
T = data[0]
ortho = data[1]*100
para = data[2]*100

plt.plot(T, ortho, label="ortho-H$_2$", color="black", linestyle=":")
plt.plot(T, para, label="para-H$_2$", color="black", linestyle="--")
plt.vlines([20, 290], 0, 100, colors=["black"], linestyles=["-"])
plt.annotate("Siedetemperatur", xytext=[9, 30], xy=[0,0], rotation=90)
plt.annotate("[Bei 1 bar]", xytext=[22, 30], xy=[0,0], rotation=90)
plt.annotate("Raumtemperatur", xytext=[278, 30], xy=[0,0], rotation=90)

leg = plt.legend(loc="lower center")
plt.ylim([0, 100])
plt.xlim([0, 300])
plt.xlabel("Temperatur $T$ [K]")
plt.ylabel("Anteil in \%")

fig = mpl.pyplot.gcf()
fig.set_size_inches(16/2.54, 9/2.54)


fig.savefig(os.path.join(save_dir, 'spin.png'), dpi=600, bbox_inches="tight")
plt.show()