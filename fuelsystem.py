# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import math
import json
import fuelflow

# def interpolate(X, Y, x: float):
#     """
#     Linear interpolation along axis X
    
#     Parameters
#     ----------
#     X : array_like 
#         sorted list of Physical properties forming the abscissa of the interpolation
#     Y : array_like 
#         sorted list of Physical properties forming the ordinate of the interpolation
#     x : float 
#         Point of interest on the abscissa

#     Returns
#     -------
#     y : float
#         Point of interest on the ordinate

#     """
#     # find the nearest point above point of interest
#     x1 = X[X > x].min()
#     # find the nearest point below point of interest
#     x0 = X[X < x].max()
#     # get indeces of the known points
#     id_1 = np.where(X == x1)
#     id_0 = np.where(X == x0)
#     # get ordinates of the known points
#     y1 = Y[id_1]
#     y0 = Y[id_0]
#     # apply linear interpolation
#     y = y0 + (y1-y0)/(x1-x0)*(x-x0)
#     return y










# jeta_props = jeta_properties(300, 1e5)
# dh = jeta_properties(400, 1e5)[1]-jeta_properties(300, 1e5)[1]
# print(jeta_props)
# print(dh)
# jeta_props = jeta_properties(400, 1e5)
# print(jeta_props)
# tpp = import_tpp()

# h2_2_cp = tpp['h2_2bar']["Cp (J/g*K)"].to_numpy()
# h2_2_t = tpp['h2_2bar']["Temperature (K)"].to_numpy()
# h2_30_cp = tpp['h2_30bar']["Cp (J/g*K)"].to_numpy()
# h2_30_t = tpp['h2_30bar']["Temperature (K)"].to_numpy()

# diff_2 = list()
# diff_30 = list()
# cp_2 = list()
# cp_30 = list()
# for i in range(np.prod(h2_2_cp.shape)):
#     cp_2.append(H2_cp(h2_2_t[i], 2*10**5)/1000)
#     diff_2.append(abs(1 - 1/(h2_2_cp[i]/cp_2[i])))
# for i in range(np.prod(h2_30_cp.shape)):
#     cp_30.append(H2_cp(h2_30_t[i], 30*10**5)/1000)
#     diff_30.append(abs(1 - 1/(h2_30_cp[i]/cp_30[i])))
# fig, axs = plt.subplots(2)


# axs[0].plot(h2_2_t, h2_2_cp)
# axs[0].plot(h2_2_t, cp_2)
# axs[1].plot(h2_30_t, h2_30_cp)
# axs[1].plot(h2_30_t, cp_30)

# axs[0].set(ylim=(10, 14))
# axs[1].set(ylim=(10, 36))
# axs[0].title.set_text("2 bar")
# axs[1].title.set_text("30 bar")
# for i in range(2):
#     axs[i].set(xlim=(0, 200), xlabel="Temperature [K]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
#     axs[i].legend(["nist.gov", "stoffmodell"])
# fig.tight_layout()
# plt.savefig("h2stoffmodell_abweichungen.png", dpi=600)

# min_p = 10**5
# max_p = 2*10**6
# p_increment = 10**4
# cp_t30 = list()
# cp_t40 = list()
# cp_t300 = list()
# cp_t600 = list()
# p_list = list()

# for i in range(min_p, max_p, p_increment):
#     p_list.append(i)
#     cp_t30.append(H2_cp(30, i)/1000)
#     cp_t40.append(H2_cp(40, i)/1000)
#     cp_t300.append(H2_cp(300, i)/1000)
#     cp_t600.append(H2_cp(600, i)/1000)

# fig, ax = plt.subplots()
# ax.plot(p_list, cp_t30)
# ax.plot(p_list, cp_t40)
# ax.plot(p_list, cp_t300)
# ax.plot(p_list, cp_t600)

# ax.legend(["T = 30K", "T = 40K", "T = 300K", "T = 600K"])
# ax.set(xlabel="Pressure [Pa]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
# ax.set_xlim(xmin=min_p, xmax=max_p)
# plt.savefig("isothermen_cp", dpi=600)

# min_t = 10
# max_t = 200
# t_increment = 1
# cp_p1 = list()
# cp_p10 = list()
# cp_p15 = list()
# cp_p20 = list()
# cp_p30 = list()
# h_p1 = list()
# h_p10 = list()
# h_p15 = list()
# h_p20 = list()
# h_p30 = list()
# t_list = list()

# for i in range(min_t, max_t, t_increment):
#     i = i/2
#     t_list.append(i)
#     cp_p1.append(H2_cp(i, 10**5)/1000)
#     cp_p10.append(H2_cp(i, 10**6)/1000)
#     cp_p15.append(H2_cp(i, 1.5*10**6)/1000)
#     cp_p20.append(H2_cp(i, 2*10**6)/1000)
#     cp_p30.append(H2_cp(i, 3*10**6)/1000)
#     h_p1.append(H2_h(i, 10**5)/1000)
#     h_p10.append(H2_h(i, 10**6)/1000)
#     h_p15.append(H2_h(i, 1.5*10**6)/1000)
#     h_p20.append(H2_h(i, 2*10**6)/1000)
#     h_p30.append(H2_h(i, 3*10**6)/1000)

# fig, ax = plt.subplots()
# ax.plot(t_list, cp_p1)
# ax.plot(t_list, cp_p10)
# ax.plot(t_list, cp_p15)
# ax.plot(t_list, cp_p20)
# ax.plot(t_list, cp_p30)

# ax.legend(["p = 1 bar", "p = 10 bar", "p = 15 bar",  "p = 20 bar", "p = 30 bar"])
# ax.set(xlabel="Temperature [K]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
# ax.set_xlim(xmin=min_t, xmax=max_t/2)
# plt.savefig("isobaren_cp", dpi=600)

# fig, ax = plt.subplots()
# ax.plot(t_list, h_p1)
# ax.plot(t_list, h_p10)
# ax.plot(t_list, h_p15)
# ax.plot(t_list, h_p20)
# ax.plot(t_list, h_p30)

# ax.legend(["p = 1 bar", "p = 10 bar", "p = 15 bar",  "p = 20 bar", "p = 30 bar"])
# ax.set(xlabel="Temperature [K]", ylabel=r"$\ h [\frac{kJ}{kg}]$ ")
# ax.set_xlim(xmin=min_t, xmax=max_t/2)
# plt.savefig("isobaren_h", dpi=600)

# t_hot = H2_find_t_for_h(H2_h(80, 1e5), 1e5)
# print(t_hot)



# ff1 = Jeta_Flow(1, 300, 1e5)
# ff1.heat(50000)
# print(ff1.t)
# ff1.heat(-50000)
# print(ff1.t)

# print(jeta_properties(300, 1e5)[0])
# print(jeta_properties(300, 3*1e6)[1]-jeta_properties(300, 1e6)[1])
# work = ff1.pump(3e6, 0.8)
# print(ff1.p, ff1.t)
# print(work)




t_r1 = 420
qm_cb = 0.3
qm_r = 0.3
qm_t = 0.1
t0 = 250
p0 = 0.4e5
p_lpfp = 3e5
eta_lpfp = 0.83
Q_fohe = 200000
p_hpfp = 3e6
eta_hpfp = 0.88
Q_idg = 5500

t_r0 = 1000
i = 0
while abs(t_r0 - t_r1) > 1e-6:
    jetaflow = fuelflow.JetaFlow(qm_cb+qm_t, t0, p0)
    i+=1
    t_r0 = t_r1
    t_lpfp = jetaflow.pump_hydraulic(p_lpfp, eta_lpfp)
    t_mix = jetaflow.mix_flows(fuelflow.JetaFlow(qm_r, t_r0, jetaflow.p))
    t_fohe = jetaflow.heat_fixed_power(Q_fohe)
    t_hpfp = jetaflow.pump_hydraulic(p_hpfp, eta_hpfp)
    cb_ff = jetaflow.split_flows(qm_cb)
    t_idg = jetaflow.heat_fixed_power(Q_idg)
    to_tank_ff = jetaflow.split_flows(qm_t)
    t_r1 = jetaflow.t
    
print("\nReference \n")
print("T_HX", t_mix)
print("T_CB ", t_hpfp[1], " iteration ", i)
print("P_v ", t_hpfp[0])
print("Q ", Q_fohe+Q_idg)


t_r1 = 320
qm_cb = 0.1
qm_r = 0.017
t0 = 20
p0 = 2e5
eta_lpfp = 0.8
Q_gsmt = 331000
p_hpfp = 3e6
eta_hpfp = 0.88

t_r0 = 1000
i = 0

while abs(t_r0 - t_r1) > 1e-6:
    i+=1
    t_r0 = t_r1
    h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_r0, p_hpfp, 1))
    # print("p_mix ", h2flow.p)
    # print("T_mix ", t_mix)
    # print("T_sat ", sat_t(h2flow.p, True))
    t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
    # print("T_HPFP ", t_hpfp)
    t_fohe = h2flow.heat_fixed_power(Q_gsmt) 
    cb_ff = h2flow.split_flows(qm_cb)
    t_r1 = h2flow.t
print("\nCompressor Pre-Mix \n")
print("T_HX", t_hpfp[1])
print("T_CB ", t_fohe, " iteration ", i)
print("P_v ", t_hpfp[0])
print("Q ", Q_gsmt)

sumpq = t_hpfp[0] + Q_gsmt
#dh = qm_cb*(H2_h(t_fohe, p_hpfp, 1)-H2_h(t0, p0, 0))
# print(sumpq, dh)


t_r1 = 320
qm_cb = 0.1
qm_r = 0.017
qm_r2 = 0.45
t0 = 20
p0 = 2e5
eta_lpfp = 0.8
Q_gsmt = 331000
p_hpfp = 3e6
eta_hpfp = 0.88

t_r0 = 1000
i = 0

while abs(t_r0 - t_r1) > 1e-6:
    i+=1
    t_r0 = t_r1
    h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_r0, p_hpfp, 1))
    # print("p_mix ", h2flow.p)
    # print("T_mix ", t_mix)
    # print("T_sat ", sat_t(h2flow.p, True))
    t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
    t_mix2 = h2flow.mix_flows(fuelflow.H2Flow(qm_r2, t_r0, p_hpfp, 1))
    # print("T_HPFP ", t_hpfp)
    t_fohe = h2flow.heat_fixed_power(Q_gsmt) 
    cb_ff = h2flow.split_flows(qm_cb)
    t_r1 = h2flow.t
print("\nCompressor with Dual Mix \n")
print("T_HX", t_mix2)
print("T_CB ", t_fohe, " iteration ", i)
print("P_v ", t_hpfp[0])
print("Q ", Q_gsmt)

sumpq = t_hpfp[0] + Q_gsmt
# dh = qm_cb*(H2_h(t_fohe, p_hpfp, 1)-H2_h(t0, p0, 0))
# print(sumpq, dh)



t_r1 = 320
qm_cb = 0.1
qm_r = 0.45
t0 = 20
p0 = 2e5
eta_lpfp = 0.8
Q_gsmt = 380000
p_hpfp = 3e6
eta_hpfp = 0.88

t_r0 = 1000
i = 0

while abs(t_r0 - t_r1) > 1e-6:
    i+=1
    t_r0 = t_r1
    h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    Q_vap = h2flow.heat_to_saturation()
    t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
    t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_r0, p_hpfp, 1))
    t_fohe = h2flow.heat_fixed_power(Q_gsmt- Q_vap) 
    cb_ff = h2flow.split_flows(qm_cb)
    t_r1 = h2flow.t
print("\nCompressor After-Mix \n")
print("T_HX", t_mix)
print("T_CB ", t_fohe, " iteration ", i)
print("P_v ", t_hpfp[0])
print("Q ", Q_gsmt)

sumpq = t_hpfp[0] + Q_gsmt
#dh = qm_cb*(H2_h(t_fohe, p_hpfp, 1)-H2_h(t0, p0, 0))
# print(sumpq, dh)

t_r1 = 320
qm_cb = 0.1
qm_r = 0.5
t0 = 20
p0 = 2e5
eta_lpfp = 0.8
Q_gsmt = 422000
p_hpfp = 3e6
eta_hpfp = 0.8

t_r0 = 1000
i = 0

while abs(t_r0 - t_r1) > 1e-6:
    i+=1
    t_r0 = t_r1
    h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    t_hpfp = h2flow.pump_hydraulic(p_hpfp, eta_hpfp)
    # print(t_hpfp)
    t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_r0, p_hpfp, 1))
    # print("p_mix ", h2flow.p)
    # print("T_mix ", t_mix)
    # print("T_sat ", sat_t(h2flow.p, True))
    t_fohe = h2flow.heat_fixed_power(Q_gsmt)
    cb_ff = h2flow.split_flows(qm_cb)
    t_r1 = h2flow.t
print("\nPump \n")
print("T_HX", t_mix)
print("T_CB ", t_fohe, " iteration ", i)
print("P_p ", t_hpfp[0])
print("Q ", Q_gsmt)

sumpq = t_hpfp[0] + Q_gsmt
# dh = qm_cb*(H2_h(t_fohe, p_hpfp, 1)-H2_h(t0, p0, 0))
# print(sumpq, dh)