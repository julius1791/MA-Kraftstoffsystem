# -*- coding: utf-8 -*-

import fuelflow
from parallel_comustion import calc_parallel_combustion

tolerance = 1e-6
max_iter = 200
v0 = 40

def reference(t_cb, qm_hpfp, Q_fohe, Q_idg, p_hpfp, eta_hpfp, p_lpfp, eta_lpfp, qm_cb, t0, p0, v = v0, tolerance = tolerance):
    t_r1 = 0
    i = 0
    qm_t = 0.1
    #h0 = fuelflow.JetaFlow(1, t0, p0).calc_h()
    #hcb = fuelflow.JetaFlow(1, t_cb, p_hpfp).calc_h()
    while abs(t_cb-t_r1) > tolerance:
        jetaflow = fuelflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        i+=1
        P_lpfp, _ = jetaflow.pump_hydraulic(p_lpfp, eta_lpfp)
        jetaflow.mix_flows(fuelflow.JetaFlow(qm_hpfp-qm_t-qm_cb, t_cb, jetaflow.p, v))
        jetaflow.heat_fixed_power(Q_fohe)
        P_hpfp, _ = jetaflow.pump_hydraulic(p_hpfp, eta_hpfp)
        jetaflow.sat_p()
        jetaflow.split_flows(qm_cb)
        jetaflow.heat_fixed_power(Q_idg)
        jetaflow.split_flows(qm_t)
        t_r1 = jetaflow.t
        qm_t = qm_t + (t_r1 - t_cb)/300
        if i > max_iter:
            return 0, 0, 0, 0, 200
    qm_r = qm_hpfp - qm_t - qm_cb
    return qm_r, qm_t, P_lpfp, P_hpfp, i

def get_dh(qm_cb, t_cb, p_cb, t0, p0, v):
    cb = fuelflow.H2Flow(qm_cb, t_cb, p_cb, v, 1)
    f0 = fuelflow.H2Flow(qm_cb, t0, p0, v, 0)
    dH = cb.qm * cb.calc_h() - f0.qm * f0.calc_h()
    return dH

def h2pump(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = 0.5
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
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, v, 0)
        P_hpfp, _ = h2flow.pump_hydraulic(p_hpfp, eta_hpfp)
        ff_r = fuelflow.H2Flow(qm_r, 200, p_hpfp, v, 1)
        ff_r.raise_to_h(h_r)
        t_hxa, _ = h2flow.mix_flows(ff_r)
        h2flow.heat_fixed_power(dH-P_hpfp-P_r)
        ff_cb = h2flow.split_flows(qm_cb)
        P_r, _ = h2flow.compress_thermic(p_hpfp, eta_r)
        ff_cb.reduce_pressure(ff_cb.p-168.9e3)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += p_cbt - p_cba
        h_r = h2flow.calc_h()
        qm_r += (t_hxt-t_hxa)/200
        print(i, abs(t_cbt - t_cba), abs(t_hxa-t_hxt), abs(p_cbt - p_cba), abs(h_r - h_rold))
        condition_bool = (
            abs(t_cbt - t_cba) > tolerance 
            and abs(t_hxa-t_hxt) > tolerance 
            and abs(p_cbt - p_cba) > tolerance 
            and abs(h_r - h_rold) > tolerance) 
        if i > max_iter:
            return 0, 0, 0, 0, 200
    print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
    return qm_r, P_hpfp, P_r, dH-P_hpfp, i

def h2after(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = 0.5
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
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, v, 0)
        Q_vap = h2flow.heat_to_saturation()
        P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        ff_r = fuelflow.H2Flow(qm_r, 200, p_hpfp, v, 1)
        ff_r.raise_to_h(h_r)
        t_hxa, _ = h2flow.mix_flows(ff_r)
        h2flow.heat_fixed_power(dH - P_hpfp - Q_vap-P_r) 
        ff_cb = h2flow.split_flows(qm_cb)
        ff_cb.reduce_pressure(ff_cb.p-168.9e3)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += (p_cbt - p_cba)
        P_r, _ = h2flow.compress_thermic(p_hpfp, eta_r)
        h_r = h2flow.calc_h()
        qm_r += (t_hxt-t_hxa)/400
        condition_bool = (
            abs(t_cbt - t_cba) > tolerance 
            and abs(t_hxa-t_hxt) > tolerance 
            and abs(p_cbt - p_cba) > tolerance 
            and abs(h_r - h_rold) > tolerance) 
        if i > max_iter:
            return 0, 0, 0, 0, 200
    print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
    return qm_r, P_hpfp, P_r, dH-P_hpfp, i

def h2dual(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = 0.5
    qm_v = 0.02
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
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, v, 0)
        ff_v = fuelflow.H2Flow(qm_v, t0, p0, v, 0)
        ff_v.raise_to_h(h_v)
        t_va, _ = h2flow.mix_flows(ff_v)
        t_vt = h2flow.sat_t() + 10
        qm_v += -(t_va - t_vt)/2000
        P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        ff_r = fuelflow.H2Flow(qm_r, 200, p_hpfp, v, 1)
        ff_r.raise_to_h(h_r)
        t_hxa, _ = h2flow.mix_flows(ff_r)
        h2flow.heat_fixed_power(dH - P_hpfp - P_r) 
        ff_cb = h2flow.split_flows(qm_cb)
        ff_cb.reduce_pressure(ff_cb.p-168.9e3)
        h_v = h2flow.calc_h()
        h2flow.split_flows(qm_v)
        P_r, _ = h2flow.compress_thermic(p_hpfp, eta_r)
        h_r = h2flow.calc_h()
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += (p_cbt - p_cba)*0.6
        qm_r += (t_hxt-t_hxa)/800
        condition_bool = (
            abs(t_cbt - t_cba) > tolerance 
            and abs(t_va-t_vt) > tolerance
            and abs(t_hxa-t_hxt) > tolerance 
            and abs(p_cbt - p_cba) > tolerance 
            and abs(h_r - h_rold) > tolerance
            and abs(h_v - h_vold) > tolerance)
        if i > max_iter:
            return 0, 0, 0, 0, 200
        if i > max_iter:
            return 0, 0, 0, 0, 0, 200
    print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
    return qm_r, qm_v, P_hpfp, P_r, dH-P_hpfp, i

def h2pre(t_cbt, t_hxt, eta_hpfp, p_cbt, qm_cb0, t0, p0, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance):
    qm_r = 0.07
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
        h2flow = fuelflow.H2Flow(qm_cb, t0, p0, v, 0)
        ff_r = fuelflow.H2Flow(qm_r, 200, p0, v, 1)
        ff_r.raise_to_h(h_r)
        h2flow.mix_flows(ff_r)
        P_hpfp, t_hxa = h2flow.compress_thermic(p_hpfp, eta_hpfp)
        h2flow.heat_fixed_power(dH-P_hpfp-P_r) 
        ff_cb = h2flow.split_flows(qm_cb)
        h_r = h2flow.calc_h()
        ff_cb.reduce_pressure(ff_cb.p-168.9e3)
        t_cba = ff_cb.t
        p_cba = ff_cb.p
        p_hpfp += p_cbt - p_cba
        qm_r += (t_hxt-t_hxa)/1600
        condition_bool = (
            abs(t_cbt - t_cba) > tolerance 
            and abs(t_hxa-t_hxt) > tolerance 
            and abs(p_cbt - p_cba) > tolerance 
            and abs(h_r - h_rold) > tolerance) 
        if i > max_iter:
            return 0, 0, 0, 0, 200
    print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
    return qm_r, P_hpfp, dH-P_hpfp, i





if __name__ == "__main__":

    t_bk = 400
    t_wu = 250
    eta_p = 0.88
    eta_r = 0.8
    p_cb = 1.5e6
    qm_cb = 0.1
    t0 = 22
    p0 = 4.2e5
    
    
    # print("reference")
    # print("qmr qmt Plp Php i")
    # qm_r , qm_t, P_lpfp, P_hpfp, i = reference(430, 1, 200000, 5500, p_cb, 0.88, 5e5, 0.83, 0.3, 250, 0.4e5)
    # print(round(qm_r, 4) , round(qm_t, 4), round(P_lpfp/1000, 3), round(P_hpfp/1000, 3), i)
    
    # print("\nh2dual")
    # print("qmr qmv Php P_r Q i")
    # qm_r1, qm_v, P_hpfp, P_r, Q, i = h2dual(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)
    # print(round(qm_r1, 4), round(qm_v, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    print("\nh2pump")
    print("qmr Php P_r Q i")
    qm_r1, P_hpfp, P_r, Q, i = h2pump(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)
    print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
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







