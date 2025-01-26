# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import fuelflow



def reference(t_cb, qm_hpfp, Q_fohe, Q_idg, p_hpfp, eta_hpfp, p_lpfp, eta_lpfp, qm_cb, t0, p0):
    t_r1 = 420
    i = 0
    qm_t = 0.1
    #h0 = fuelflow.JetaFlow(1, t0, p0).calc_h()
    #hcb = fuelflow.JetaFlow(1, t_cb, p_hpfp).calc_h()
    while abs(t_cb-t_r1) > 1e-6:
        jetaflow = fuelflow.JetaFlow(qm_cb+qm_t, t0, p0)
        i+=1
        P_lpfp, _ = jetaflow.pump_hydraulic(p_lpfp, eta_lpfp)
        jetaflow.mix_flows(fuelflow.JetaFlow(qm_hpfp-qm_t-qm_cb, t_cb, jetaflow.p))
        jetaflow.heat_fixed_power(Q_fohe)
        P_hpfp, _ = jetaflow.pump_hydraulic(p_hpfp, eta_hpfp)
        jetaflow.sat_p()
        jetaflow.split_flows(qm_cb)
        jetaflow.heat_fixed_power(Q_idg)
        jetaflow.split_flows(qm_t)
        t_r1 = jetaflow.t
        qm_t = qm_t + (t_r1 - t_cb)/300
    qm_r = qm_hpfp - qm_t - qm_cb
    return qm_r, qm_t, P_lpfp, P_hpfp, i

def get_dh(qm_cb, t_cb, p_hpfp, t0, p0):
    cb = fuelflow.H2Flow(qm_cb, t_cb, p_hpfp, 1)
    f0 = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    dH = cb.qm * cb.calc_h() - f0.qm * f0.calc_h()
    return dH

def h2pump(t_cb, t_hx, eta_hpfp, p_hpfp, qm_cb, t0, p0):
    qm_r1 = 0.5
    h_mix = fuelflow.H2Flow(1, t_hx, p_hpfp, 1).calc_h()
    dH = get_dh(qm_cb, t_cb, p_hpfp, t0, p0)
    t_r1 = 300
    i = 0
    t_mix = 0
    while abs(t_cb - t_r1) + abs(t_mix-t_hx) > 1e-6:
        i+=1
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
        P_hpfp, _ = h2flow.pump_hydraulic(p_hpfp, eta_hpfp)
        H_0 = h2flow.calc_h()*h2flow.qm
        t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r1, t_cb, p_hpfp, 1))
        H_mix = h_mix*h2flow.qm
        h2flow.heat_fixed_power(dH-P_hpfp)
        h2flow.split_flows(qm_cb)
        t_r1 = h2flow.t
        h_r = h2flow.calc_h()
        qm_r1 = (H_mix-H_0)/h_r
    return qm_r1, P_hpfp, dH-P_hpfp, i

def h2after(t_cb, t_hx, eta_hpfp, p_hpfp, qm_cb, t0, p0):
    qm_r1 = 0.45
    h_mix = fuelflow.H2Flow(1, t_hx, p_hpfp, 1).calc_h()
    dH = get_dh(qm_cb, t_cb, p_hpfp, t0, p0)
    t_r1 = 300
    i = 0
    t_mix = 0
    while abs(t_cb - t_r1) + abs(t_mix-t_hx) > 1e-6:
        i+=1
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
        Q_vap = h2flow.heat_to_saturation()
        P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        H_0 = h2flow.calc_h()*h2flow.qm
        t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r1, t_cb, p_hpfp, 1))
        H_mix = h_mix*h2flow.qm
        h2flow.heat_fixed_power(dH - P_hpfp - Q_vap) 
        h2flow.split_flows(qm_cb)
        t_r1 = h2flow.t
        h_r = h2flow.calc_h()
        qm_r1 = (H_mix-H_0)/h_r
    return qm_r1, P_hpfp, dH-P_hpfp, i

def h2dual(t_cb, t_hx, eta_hpfp, p_hpfp, qm_cb, t0, p0):
    qm_r1 = 0.45
    h_mix = fuelflow.H2Flow(1, t_hx, p_hpfp, 1).calc_h()
    dH = get_dh(qm_cb, t_cb, p_hpfp, t0, p0)
    t_r1 = 300
    i = 0
    
    vap = fuelflow.H2Flow(1, 300, p_hpfp, 1)
    vap.t = vap.sat_t() + 5
    h_vap = vap.calc_h()
    H_0 = fuelflow.H2Flow(qm_cb, t0, p0, 0).calc_h()*qm_cb
    h_r = fuelflow.H2Flow(1, t_cb, p_hpfp, 1).calc_h()
    qm_v = (H_0 - qm_cb*h_vap)/(h_vap-h_r)
    t_mix2 = 0
    while abs(t_cb - t_r1) + abs(t_mix2-t_hx) > 1e-6:
        i+=1
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
        h2flow.mix_flows(fuelflow.H2Flow(qm_v, t_cb, p_hpfp, 1))
        P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        H_0 = h2flow.qm*h2flow.calc_h()
        t_mix2 = h2flow.mix_flows(fuelflow.H2Flow(qm_r1, t_cb, p_hpfp, 1))
        H_mix = h_mix*h2flow.qm
        h2flow.heat_fixed_power(dH - P_hpfp) 
        h2flow.split_flows(qm_cb)
        t_r1 = h2flow.t
        h_r = h2flow.calc_h()
        qm_r1 = (H_mix-H_0)/h_r
    return qm_r1, qm_v, P_hpfp, dH-P_hpfp, i

def h2pre(t_cb, t_hx, eta_hpfp, p_hpfp, qm_cb, t0, p0):
    qm_r1 = 0.02
    dH = get_dh(qm_cb, t_cb, p_hpfp, t0, p0)
    t_r1 = 300
    i = 0
    t_hpfp = 0
    while abs(t_cb - t_r1) + abs(t_hpfp-t_hx)> 1e-6:
        i+=1
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
        h2flow.mix_flows(fuelflow.H2Flow(qm_r1, t_cb, p_hpfp, 1))
        P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        h2flow.heat_fixed_power(dH-P_hpfp) 
        h2flow.split_flows(qm_cb)
        qm_r1 = qm_r1 + (t_hx - t_hpfp)/10000
        t_r1 = h2flow.t
    return qm_r1, P_hpfp, dH-P_hpfp, i


t_bk = 300
t_wu = 250
eta_p = 0.88
p_h = 3e6
qm_cb = 0.1
t0 = 22
p0 = 2e5


qm_r , qm_t, P_lpfp, P_hpfp, i = reference(430, 1, 200000, 5500, 3e6, 0.88, 5e5, 0.83, 0.3, 250, 0.4e5)
print(qm_r , qm_t, P_lpfp, P_hpfp, i)

qm_r1, qm_v, P_hpfp, Q, i = h2dual(t_bk, t_wu, eta_p, p_h, qm_cb, t0, p0)
print(qm_r1, qm_v, P_hpfp, Q, i)

qm_r1, P_hpfp, Q, i = h2pump(t_bk, t_wu, eta_p, p_h, qm_cb, t0, p0)
print(qm_r1, P_hpfp, Q, i)

qm_r1, P_hpfp, Q, i = h2after(t_bk, t_wu, eta_p, p_h, qm_cb, t0, p0)
print(qm_r1, P_hpfp, Q, i)

qm_r1, P_hpfp, Q, i = h2pre(t_bk, t_wu, eta_p, p_h, qm_cb, t0, p0)
print(qm_r1, P_hpfp, Q, i)

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











