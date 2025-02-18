# -*- coding: utf-8 -*-

import h2flow
import jetaflow
from parallel_comustion import calc_parallel_combustion

tolerance = 1e-3
max_iter = 500
rel_fac = 1/80
v0 = 20

def reference(t_cbt, qm_hpfp, Q_fohe, Q_idg, pt_cbt, eta_hpfp, p_cav, eta_lpfp, qm_cb, t0, p0, v = v0, tolerance = tolerance):
    ht_r = jetaflow.calc_ht(t_cbt, pt_cbt, v)
    p_lpfp = 2e5
    i = 0
    qm_t = 0.1
    qm_t0 = 0.1
    p_hpfp = pt_cbt
    #h0 = fuelflow.JetaFlow(1, t0, p0).calc_h()
    #hcb = fuelflow.JetaFlow(1, t_cb, p_hpfp).calc_h()
    condition_bool = True
    while condition_bool:
        i+=1
        ht_rold = ht_r
        ff_main = jetaflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        P_lpfp, _ = ff_main.pump_hydraulic(p_lpfp, eta_lpfp)
        ff_main.mix_flows(jetaflow.JetaFlow(qm_hpfp-qm_t-qm_cb, jetaflow.calc_t(ht_r, ff_main.p, v), ff_main.p, v))
        ff_main.heat_exchanger(Q_fohe)
        p_sat_pi = jetaflow.sat_p(ff_main.t)
        p_pi = ff_main.p
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        pt_cba = jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
        t_cba = ff_main.t
        ff_main.split_flows(qm_cb)
        ff_main.heat_exchanger(Q_idg)
        ff_main.split_flows(qm_t)
        p_hpfp += (pt_cbt - pt_cba)
        ht_r = ff_main.ht
        qm_t += qm_t0 * (t_cba - t_cbt) * rel_fac
        p_lpfp += (p_sat_pi - p_pi + p_cav)
        #print(p_lpfp, p_hpfp)
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(pt_cbt - pt_cba) < tolerance 
            and abs(ht_r - ht_rold) < tolerance
        )
        if i > max_iter:
            return 0, 0, 0, 0, i
    qm_r = qm_hpfp - qm_t - qm_cb
    #print("saturation margin [Pa]: " + str(p_pi - p_sat_pi))
    return qm_r, qm_t, P_lpfp, P_hpfp, i

def get_dh(qm_cb, t_cb, p_cb, t0, p0, v):
    cb = h2flow.H2Flow(qm_cb, t_cb, p_cb, v, True)
    f0 = h2flow.H2Flow(qm_cb, t0, p0, v, False)
    dH = cb.qm * h2flow.calc_ht(cb.t, cb.p, cb.v) - f0.qm * h2flow.calc_ht(f0.t, f0.p, f0.v)
    return dH

def h2pump(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v) 
    dH = dH0
    h_r = dH/qm_cb0
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    while condition_bool:
        i+=1
        if pcc:
            qm_cb = qm_cb0 + calc_parallel_combustion(dH-P_r-P_hpfp-Q_hx, t_cbt, p_cbt, 344)
            dH = dH0 * qm_cb / qm_cb0
        else:
            qm_cb = qm_cb0
        h_rold = h_r
        ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_hpfp, v, True), p_hpfp, v, True)
        t_hxa, _ = ff_main.mix_flows(ff_r)
        ff_main.heat_exchanger(dH-P_hpfp-P_r)
        ff_cb = ff_main.split_flows(qm_cb)
        P_r, _ = ff_main.pump_hydraulic(p_hpfp, eta_r)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += p_cbt - p_cba
        h_r = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
        qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(t_hxa-t_hxt) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - h_rold) < tolerance
        )
        if i > max_iter:
            return 0, 0, 0, 0, i
    return qm_r, P_hpfp, P_r, dH-P_hpfp, i

def h2after(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    while condition_bool:
        i+=1
        if pcc:
            qm_cb = qm_cb0 + calc_parallel_combustion(dH-P_r-P_hpfp-Q_hx, t_cbt, p_cbt, 344)
            dH = dH0 * qm_cb / qm_cb0
        else:
            qm_cb = qm_cb0
        h_rold = h_r
        ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
        Q_vap = ff_main.heat_to_saturation()
        P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_hpfp, v, True), p_hpfp, v, True)
        t_hxa, _ = ff_main.mix_flows(ff_r)
        ff_main.heat_exchanger(dH - P_hpfp - Q_vap  -P_r) 
        ff_cb = ff_main.split_flows(qm_cb)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += (p_cbt - p_cba)
        P_r, _ = ff_main.pump_hydraulic(p_hpfp, eta_r)
        h_r = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
        qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
        #print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(t_hxa-t_hxt) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - h_rold) < tolerance
        )
        if i > max_iter:
            return 0, 0, 0, 0, i
    return qm_r, P_hpfp, P_r, dH-P_hpfp, i

def h2dual(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_v = qm_cb0 * (h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v))
    qm_v0 = qm_v
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v)) - qm_v
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    h_v = h_r
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    
    while condition_bool:
        i+=1
        if pcc:
            qm_cb = qm_cb0 + calc_parallel_combustion(dH-P_r-P_hpfp-Q_hx, t_cbt, p_cbt, 344)
            dH = dH0 * qm_cb / qm_cb0
        else:
            qm_cb = qm_cb0
        h_rold = h_r
        h_vold = h_v
        ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
        ff_v = h2flow.H2Flow(qm_v, h2flow.calc_t(h_v, p_hpfp, v, True), p_hpfp, v, True)
        t_va, _ = ff_main.mix_flows(ff_v)
        t_vt = h2flow.sat_t(ff_main.p) + 10
        P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_hpfp, v, True), p_hpfp, v, True)
        t_hxa, _ = ff_main.mix_flows(ff_r)
        ff_main.heat_exchanger(dH - P_hpfp - P_r) 
        ff_cb = ff_main.split_flows(qm_cb)
        h_v = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
        ff_main.split_flows(qm_v)
        P_r, _ = ff_main.pump_hydraulic(p_hpfp, eta_r)
        h_r = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += (p_cbt - p_cba)*0.6
        qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
        qm_v += qm_v0*(t_vt - t_va) * rel_fac
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(t_va-t_vt) < tolerance
            and abs(t_hxa-t_hxt) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - h_rold) < tolerance
            and abs(h_v - h_vold) < tolerance)
        if i > max_iter:
            return 0, 0, 0, 0, 0, i
    return qm_r, qm_v, P_hpfp, P_r, dH-P_hpfp, i


def h2pre(t_cbt, t_hxt, eta_hpfp, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = qm_cb0 * (h2flow.calc_ht(h2flow.sat_t(p0)+5, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(h2flow.sat_t(p0)+5, p0, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    while condition_bool:
        i+=1
        if pcc:
            qm_cb = qm_cb0 + calc_parallel_combustion(dH-P_r-P_hpfp-Q_hx, t_cbt, p_cbt, 344)
            dH = dH0 * qm_cb / qm_cb0
        else:
            qm_cb = qm_cb0
        h_rold = h_r
        ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
        ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_hpfp, v, True), p_hpfp, v, True)
        t_v, _ = ff_main.mix_flows(ff_r)
        P_hpfp, t_hxa = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        ff_main.heat_exchanger(dH-P_hpfp-P_r) 
        ff_cb = ff_main.split_flows(qm_cb)
        h_r = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += p_cbt - p_cba
        qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
        #qm_r += + qm_r0 * (h2flow.sat_t(p0)+5 - t_v)/50
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            # and abs(h2flow.sat_t(p0)+5 - t_v) < tolerance 
            and abs(t_hxa-t_hxt) < tolerance
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - h_rold) < tolerance
        ) 
        if i > max_iter:
            return 0, 0, 0, i
    return qm_r, P_hpfp, dH-P_hpfp, i





if __name__ == "__main__":

    t_bk = 280
    t_wu = 250
    eta_p = 0.92
    eta_r = 0.9
    p_cb = 1.5e6+168.9e3
    qm_cb = 0.1
    t0 = 22
    p0 = 4.2e5
    
    
    print("reference")
    print("qmr qmt Plp Php i")
    qm_r , qm_t, P_lpfp, P_hpfp, i = reference(430, 1, 200000, 5500, p_cb, 0.88, 2e4, 0.83, 0.3, 250, 0.4e5)
    print(round(qm_r, 4) , round(qm_t, 4), round(P_lpfp/1000, 3), round(P_hpfp/1000, 3), i)
    
    # print("\nh2dual")
    # print("qmr qmv Php P_r Q i")
    # qm_r1, qm_v, P_hpfp, P_r, Q, i = h2dual(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)
    # print(round(qm_r1, 4), round(qm_v, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2pump")
    # print("qmr Php P_r Q i")
    # qm_r1, P_hpfp, P_r, Q, i = h2pump(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2pump_pcc")
    # print("qmr Php P_r Q i")
    # qm_r1, P_hpfp, P_r, Q, i = h2pump(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, 200e3, True)
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2after")
    # print("qmr Php P_r Q i")
    # qm_r1, P_hpfp, P_r, Q, i = h2after(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2pre")
    # print("qmr Php Q i")
    # qm_r1, P_hpfp, Q, i = h2pre(t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0)
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(Q/1000, 3), i)







