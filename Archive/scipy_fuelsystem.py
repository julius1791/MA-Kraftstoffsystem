# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:21:13 2025

@author: Julius
"""

import scipy.optimize
import fuelflow


def get_dh(qm_cb, t_cb, p_cb, t0, p0):
    cb = fuelflow.H2Flow(qm_cb, t_cb, p_cb, 1)
    f0 = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    dH = cb.qm * cb.calc_h() - f0.qm * f0.calc_h()
    return dH




def h2dual(x, *k):
    (t_rt, p_rt, qm_v, qm_r, p_hpfp) = x
    (t_cbt, t_hxt, eta_hpfp, p_cbt, qm_cb, t0, p0) = k
    dH = get_dh(qm_cb, t_cbt, p_cbt, t0, p0)

    h2flow = fuelflow.H2Flow(qm_cb, t0, p0, 0)
    t_va, _ = h2flow.mix_flows(fuelflow.H2Flow(qm_v, t_cbt, p_cbt, 1))
    t_vt = h2flow.sat_t() + 10
    P_hpfp, t_hpfp = h2flow.compress_thermic(p_hpfp, eta_hpfp)
    t_hxa, _ = h2flow.mix_flows(fuelflow.H2Flow(qm_r, t_rt, p_rt, 1))
    h2flow.heat_fixed_power(dH - P_hpfp) 
    ff_cb = h2flow.split_flows(qm_cb)
    ff_cb.reduce_pressure(ff_cb.p-168.9e3)
    #t_cba = ff_cb.t
    p_cba = ff_cb.p
    t_ra = h2flow.t
    p_ra = h2flow.p
    
    f = [t_hxt - t_hxa]
    #f.append(t_cbt - t_cba)
    f.append(p_cbt - p_cba)
    f.append(t_rt - t_ra)
    f.append(p_rt - p_ra)
    f.append(t_vt - t_va)
    
    return f

t_bk = 270
t_wu = 200
eta_p = 0.88
p_cb = 1331391
qm_cb = 0.1
t0 = 22
p0 = 4.2e5

k = ((t_bk, t_wu, eta_p, p_cb, qm_cb, t0, p0))
print("\nh2dual")
print("qmr qmv Php Q i")

bnds = ((50, 500), (1e5, 3e6), (0, 1), (0, 1), (1e5, 3e6))

out = scipy.optimize.newton(h2dual, [260, 1.4e6, 0.03, 0.4, 2e6], args=k, tol=1e-04)
