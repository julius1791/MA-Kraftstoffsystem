import matplotlib.pyplot as plt
import fuelflow





def dualmix(qm_cb, qm_r, qm_r2, t0, p0, eta_hpfp, Q_gsmt, p_hpfp):
    t_r1 = 320
    i = 0
    t_r0 = 0
    while abs(t_r0 - t_r1) > 1e-6:
        i+=1
        t_r0 = t_r1
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
        t_mix = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_r0, p_hpfp, 1))
        # print("p_mix ", h2flow.p)
        # print("T_mix ", t_mix)
        # print("T_sat ", sat_t(h2flow.p, True))
        P_c, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        t_mix2 = h2flow.mix_flows(fuelflow.H2Flow(qm_r2, t_r0, p_hpfp, 1))
        # print("T_HPFP ", t_hpfp)
        t_fohe = h2flow.heat_fixed_power(Q_gsmt) 
        cb_ff = h2flow.split_flows(qm_cb)
        t_r1 = h2flow.t
    return P_c, t_fohe

qm_cb = 0.1
qm_r = 0.017
qm_r2 = 0.45
t0 = 20
p0 = 2e5
Q_gsmt = 331000
p_hpfp = 3e6
eta_hpfp = 0.88

Q_list = range(280, 450, 10)

P_list = list()
T_list = list()
for i, Q in enumerate(Q_list):
    print(i,Q)
    P, T = dualmix(qm_cb, qm_r, qm_r2, t0, p0, eta_hpfp, Q*1e3, p_hpfp)
    P_list.append(P/1e3)
    T_list.append(T)


fig, axs = plt.subplots(2)
axs[0].plot(Q_list, T_list)
axs[1].plot(Q_list, P_list)

axs[0].set(xlabel="Wärme [kW]", ylabel="Kraftstofftemperatur [K]")
axs[1].set(xlabel="Wärme [kW]", ylabel="Verdichterleistung [kW]")
plt.savefig("Dual_mix.png", dpi=600)