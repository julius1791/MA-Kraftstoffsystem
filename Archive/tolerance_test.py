# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 23:32:09 2025

@author: Julius
"""

import fuelsystem
import matplotlib.pyplot as plt

tolerance_list = range(-6, 4, 1)
ref_i = list()
pre_i = list()
aft_i = list()
dual_i = list()
pump_i = list()
ref_p = list()
pre_p = list()
aft_p = list()
dual_p = list()
pump_p = list()

t_bk = 270
t_wu = 200
eta_p = 0.88
p_cb = 1331391
qm_cb = 0.1
t0 = 22
p0 = 4.2e5

first = True
for tol in tolerance_list:
    delta = 10**tol
    print(delta)

    _ , _, _, P_hpfp, i = fuelsystem.reference(430, 1, 200000, 5500, p_cb, 0.88, 5e5, 0.83, 0.3, 250, 0.4e5, delta)
    ref_i.append(i)
    ref_p.append(P_hpfp)

    _, _, P_hpfp, _, i = fuelsystem.h2dual(t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0, delta)
    dual_i.append(i)
    dual_p.append(P_hpfp)

    _, P_hpfp, _, i = fuelsystem.h2pump(t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0, delta)
    pump_i.append(i)
    pump_p.append(P_hpfp)

    _, P_hpfp, _, i = fuelsystem.h2after(t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0, delta)
    aft_i.append(i)
    aft_p.append(P_hpfp)

    _, P_hpfp, _, i = fuelsystem.h2pre(t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0, delta)
    pre_i.append(i)
    pre_p.append(P_hpfp)

ref = ref_p[0]
dual = dual_p[0]
pump = pump_p[0]
aft = aft_p[0]
pre = pre_p[0]

ref_p = [x / ref for x in ref_p]
dual_p = [x / dual for x in dual_p]
pump_p = [x / pump for x in pump_p]
pre_p = [x / pre for x in pre_p]
aft_p = [x / aft for x in aft_p]

fig, ax = plt.subplots(2)

ax[0].set_title("Difference in Pump/Compressor Power (normalised)")
ax[0].set_xlabel("Tolerance (logarithmic)")
ax[0].set_ylabel("relative Pump/Compressor Power")
ax[0].plot(tolerance_list, ref_p, label="reference")
ax[0].plot(tolerance_list, dual_p, label = "dual")
ax[0].plot(tolerance_list, aft_p, label = "after")
ax[0].plot(tolerance_list, pre_p, label = "pre")
ax[0].plot(tolerance_list, pump_p, label = "pump")

ax[1].set_title("Number of Iterations required")
ax[1].set_xlabel("Tolerance (logarithmic)")
ax[1].set_ylabel("number of iterations")
ax[1].plot(tolerance_list, ref_p, label="reference")
ax[1].plot(tolerance_list, dual_i, label = "dual")
ax[1].plot(tolerance_list, aft_i, label = "after")
ax[1].plot(tolerance_list, pre_i, label = "pre")
ax[1].plot(tolerance_list, pump_i, label = "pump")
ax[0].legend(loc="center left")
fig.tight_layout()
fig.show()
