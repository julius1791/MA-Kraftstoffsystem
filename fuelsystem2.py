# -*- coding: utf-8 -*-

import h2flow
import jetaflow
from parallel_comustion import calc_parallel_combustion
import csv
import os
import time

# results folder
results_dir = "results2"

tolerance = 1e-3
max_iter = 500
rel_fac = 1/750
rel_fac_2 = 1/2
v0 = 0

# hydrogen lower heating value (250 K) / (200 K)
lhv_h2_250 = 117.24e6  # J/kg
lhv_h2_200 = (
    lhv_h2_250 - h2flow.h2.calc_H2_enthalpy(250, 1e6) 
    + h2flow.h2.calc_H2_enthalpy(200, 1e6)
)

# jeta lower heating value (288 K)
lhv_jeta_288 = 43.10e6  # J/kg


def save_results(
        filename, arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0,
        tpr_hx, 
        Q_hx, pcc, v, P_hpfp, P_r, Q, Q_phc, qm_cb, qm_pch, qm_r, qm_v, p_hpfp,
        i, duration
    ):
    path = os.path.join(os.getcwd(), results_dir, filename)
    with open(path, "w", newline='') as f:
        filewriter = csv.writer(f)
        filewriter.writerow([
            "architecture", "t_cbt", "t_hxt", "eta_hpfp", "eta_r", "p_cbt",
            "qm_cb0", "t0", "p0", "tpr_hx", "Q_hx", "pcc", "v"
            ])
        filewriter.writerow([
            arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt,
            qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v
            ])
        filewriter.writerow(["P_hpfp", "P_rv", "p_hpfp"])
        filewriter.writerow([P_hpfp, P_r, p_hpfp])
        filewriter.writerow(["Q", "Q_phc"])
        filewriter.writerow([Q, Q_phc])
        filewriter.writerow(["qm_cb", "qm_phc/qm_t", "qm_r", "qm_v"])
        filewriter.writerow([qm_cb0, qm_pch, qm_r, qm_v])
        filewriter.writerow([
            "number of iterations", i, "Execution time: ", duration
            ])
    return

def save_failed(
        filename, arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v, exception
    ):
    if len(filename) > 1:
        failed = filename[:-4] + "FAILED" + ".csv"
        path = os.path.join(os.getcwd(), results_dir, failed)
        with open(path, "w", newline='') as f:
            filewriter = csv.writer(f)
            filewriter.writerow(["architecture", "t_cbt", "t_hxt", "eta_hpfp", "eta_r", "p_cbt", "qm_cb0", "t0", "p0", "tpr_hx", "Q_hx", "pcc", "v"])
            filewriter.writerow([arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v])
            filewriter.writerow(["FAILED TO CONVERGE"])
            filewriter.writerow([exception])
    return

def reference(
        t_cbt, qm_hpfp, Q_fohe, Q_idg, tpr_fohe, tpr_idg, pt_cbt, eta_hpfp, p_cav,
        eta_lpfp, qm_cb, t0, p0, v = v0, tolerance = tolerance, p_lpfp0=5e5, 
        filename = ""
):
    start = time.time()
    # correct combustion chamber fuel mass flow for temperature
    lhv_tcbt = (
        lhv_jeta_288 
        + jetaflow.jeta.calc_jeta_enthalpy(t_cbt, 1e6) 
        - jetaflow.jeta.calc_jeta_enthalpy(288, 1e6)
        )
    qm_cb = qm_cb * (lhv_jeta_288 / lhv_tcbt)
    # calculate intitial guess for recirculation flow enthalpy
    ht_r = jetaflow.calc_ht(t_cbt, pt_cbt, v)
    
    # guess initial value for fuel mass flow return to tank
    qm_t = 0.1
    qm_t0 = qm_t
    
    # set initial value for pump pressures
    p_hpfp = pt_cbt
    p_lpfp = p_lpfp0
    
    # init loop variables
    i = 0
    condition_bool = True
    while condition_bool:
        i+=1
        ht_rold = ht_r
        # fuel from boost pump
        ff_main = jetaflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        # calculation of lp fuel pump
        P_lpfp, _ = ff_main.pump_hydraulic(p_lpfp, eta_lpfp)
        if qm_cb + qm_t > qm_hpfp:
            raise ValueError("Mass flow exceeds HPFP limit. QM_T: " + str(qm_t))
        # mix of main flow with recirculation flow
        ff_r = jetaflow.JetaFlow(qm_hpfp-qm_t-qm_cb,
            jetaflow.calc_t(ht_r, ff_main.p, v), ff_main.p, v
            )
        ff_main.mix_flows(ff_r)
        # calculation of primary heat exchanger
        ff_main.heat_exchanger(Q_fohe, tpr_fohe)
        # saturation pressure and static pressure at hp pump inlet
        p_sat_pi = jetaflow.sat_p(ff_main.t)
        p_pi = ff_main.p
        # calculation of hp fuel pump
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        # total pressure and temperature at hp pump outlet
        pt_cba = jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
        t_cba = ff_main.t
        p_sat_cb = jetaflow.sat_p(t_cba)
        ff_main.split_flows(qm_cb)
        ff_main.heat_exchanger(Q_idg, tpr_idg)
        ht_r = ff_main.ht
        
        p_hpfp_old = p_hpfp
        p_hpfp += (pt_cbt - pt_cba) * rel_fac_2
        
        qm_t += qm_t0 * (t_cba - t_cbt) * rel_fac * 4
        qm_t = max(0, qm_t)
        qm_t = min(qm_t, 0.7)
        
        if p_hpfp > 1.1 * p_hpfp_old :
            p_hpfp = 1.1*p_hpfp_old
        elif p_hpfp < 0.9 * p_hpfp_old:
            p_hpfp = 0.9*p_hpfp_old 
        
        print(i, t_cbt - t_cba, pt_cbt - pt_cba, ht_r - ht_rold)
        
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(pt_cbt - pt_cba) < tolerance 
            and abs(ht_r - ht_rold) < tolerance
        )
        if i > max_iter:
            print("failed to converge")
            return
    qm_r = qm_hpfp - qm_t - qm_cb
    print("saturation margin hpp inlet [bar]: " + str((p_pi - p_sat_pi)/1e5))
    print("saturation margin injector [bar]: " + str((pt_cba - 168.9e3 - p_sat_cb)/1e5))
    stop = time.time()
    if qm_r < 0 or qm_t < 0:
        raise ValueError("Solution includes negative mass flow")
    save_results(
        filename, "reference", t_cbt, float("nan"), eta_hpfp, eta_r, pt_cbt, 
        qm_cb, t0, p0, tpr_fohe, float("nan"), float("nan"), v, 
        P_hpfp, 0, Q_idg + Q_fohe, 0, qm_cb, qm_t,
        qm_r, p_hpfp, 0, i,
        stop-start
    )
    return

def get_dh(qm_cb, t_cb, p_cb, t0, p0, v):
    cb = h2flow.H2Flow(qm_cb, t_cb, p_cb, v, True)
    f0 = h2flow.H2Flow(qm_cb, t0, p0, v, False)
    dH = cb.qm * h2flow.calc_ht(cb.t, cb.p, cb.v) - f0.qm * h2flow.calc_ht(f0.t, f0.p, f0.v)
    return dH

def h2pump(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance, filename=""):
    start = time.time()
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v) 
    dH = dH0
    h_r = dH/qm_cb0
    ht_cb = h2flow.calc_ht(t_cbt, p_cbt, v)
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    try:
        while condition_bool:
            i+=1
            if pcc:
                qm_cb = qm_cb0 + calc_parallel_combustion(max(0, dH-P_r-P_hpfp-Q_hx), t_cbt, p_cbt, 344)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
            h_rold = h_r
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(ht_cb, p_cbt, v, True), p_cbt, v, True)
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            t_hxa, _ = ff_main.mix_flows(ff_r)
            ff_main.heat_exchanger(dH-P_hpfp-P_r, tpr_hx)
            ff_cb = ff_main.split_flows(qm_cb)
            ht_cb = h2flow.calc_ht(t_cbt, ff_main.p, v)
            p_cba = h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
            
            t_cba = ff_cb.t
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            h_r = h2flow.calc_ht(ff_main.t, ff_main.p, ff_main.v)
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
            
            
            # print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
            # print(qm_r)
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
            )
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_hx, Q_hx, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_hx, qm_cb, qm_cb-qm_cb0, qm_r, 0, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v, e)
        return
    return

def h2after(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, tpr_vap=0.96, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance, filename=""):
    start = time.time()
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    ht_cb = h2flow.calc_ht(t_cbt, p_cbt, v)
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    try:
        while condition_bool:
            i+=1
            if pcc:
                qm_cb = qm_cb0 + calc_parallel_combustion(max(0, dH-P_r-P_hpfp-Q_hx), t_cbt, p_cbt, 344)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
            h_rold = h_r
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            Q_vap = ff_main.heat_to_saturation(tpr_vap)
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(ht_cb, p_cbt, v, True), p_cbt, v, True)
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            t_hxa, _ = ff_main.mix_flows(ff_r)
            ff_main.heat_exchanger(dH - P_hpfp  -P_r, tpr_hx) 
            ff_main.heat_exchanger(- Q_vap, tpr_vap)
            ff_cb = ff_main.split_flows(qm_cb)
            ht_cb = h2flow.calc_ht(t_cbt, ff_main.p, v)
            t_cba = ff_cb.t
            p_cba = h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0.001, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
            )
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_hx, Q_hx, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_hx, qm_cb, qm_cb-qm_cb0, qm_r, 0, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v, e)
        return
    return

def h2dual(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance, filename=""):
    start = time.time()
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    qm_v = qm_cb0 * (h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v))
    qm_v0 = qm_v
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v)) - qm_v
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    h_v = h_r
    ht_cb = h2flow.calc_ht(t_cbt, p_cbt, v)
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    try:
        while condition_bool:
            i+=1
            if pcc:
                qm_cb = qm_cb0 + calc_parallel_combustion(max(0, dH-P_r-P_hpfp-Q_hx), t_cbt, p_cbt, 344)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
            h_rold = h_r
            h_vold = h_v
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            ff_v = h2flow.H2Flow(qm_v, h2flow.calc_t(ht_cb, p_cbt, v, True), p_cbt, v, True)
            t_va, _ = ff_main.mix_flows(ff_v)
            t_vt = h2flow.sat_t(ff_main.p) + 10
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(ht_cb, p_cbt, v, True), p_cbt, v, True)
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            t_hxa, _ = ff_main.mix_flows(ff_r)
            ff_main.heat_exchanger(dH - P_hpfp - P_r, tpr_hx) 
            p_cba = h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
            t_cba = ff_main.t
            ht_cb = h2flow.calc_ht(t_cbt, ff_main.p, v)
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_v += qm_v0*(t_vt - t_va) * rel_fac
            qm_r = max(0, qm_r)
            qm_v = max(0, qm_v)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_va-t_vt) < tolerance
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
                and abs(h_v - h_vold) < tolerance)
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0 or qm_v < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_hx, Q_hx, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_hx, qm_cb, qm_cb-qm_cb0, qm_r, qm_v, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v, e)
        return
    
    return 


def h2pre(t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx = 0, pcc=False, v = v0, tolerance = tolerance, filename=""):
    start = time.time()
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    t_mix = t_hxt/(1+((p_cbt/p0)**0.286-1)/eta_hpfp)
    qm_r = qm_cb0 * (h2flow.calc_ht(t_mix, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_mix, p0, v))
    qm_r0 = qm_r
    p_hpfp = p_cbt
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    ht_cb = h2flow.calc_ht(t_cbt, p_cbt, v)
    i = 0
    P_r = 0
    P_hpfp = 0
    condition_bool = True
    try:
        while condition_bool:
            i+=1
            if pcc:
                qm_cb = qm_cb0 + calc_parallel_combustion(max(0, dH-P_r-P_hpfp-Q_hx), t_cbt, p_cbt, 344)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
            h_rold = h_r
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(ht_cb, p_cbt, v, True), p_cbt, v, True)
            t_v, _ = ff_main.mix_flows(ff_r)
            P_hpfp, t_hxa = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            ff_main.heat_exchanger(dH-P_hpfp-P_r, tpr_hx) 
            t_cba = ff_main.t
            p_cba = h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)
            ht_cb = h2flow.calc_ht(t_cbt, ff_main.p, v)
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance
                and abs(t_hxa-t_hxt) < tolerance
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
            ) 
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "pre", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_hx, Q_hx, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_hx, qm_cb, qm_cb-qm_cb0, 0, qm_r, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "pre", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, tpr_hx, p0, Q_hx, pcc, v, e)
        return
    return qm_r, P_hpfp, dH-P_hpfp, i





if __name__ == "__main__":

    t_bk = 600
    t_wu = 100
    eta_p = 0.92
    eta_r = 0.9
    p_cb = 1.5e6+168.9e3
    qm_cb = 0.1
    t0 = 22
    p0 = 4.2e5
    
    
    # print("reference")
    # print("qmr qmt Plp Php i")
    # qm_r , qm_t, P_lpfp, P_hpfp, i = reference(470, 1, 200000, 5500, 0.95, 0.98, p_cb, 0.88, 2e4, 0.83, 0.3, 250, 0.4e5)
    # print(round(qm_r, 4) , round(qm_t, 4), round(P_lpfp/1000, 3), round(P_hpfp/1000, 3), i)  
    
    # print("\nh2dual")
    # print("qmr qmv Php P_r Q i")
    # h2dual(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, 0.95, filename="test.csv")
    # print(round(qm_r1, 4), round(qm_v, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    print("\nh2pump")
    print("qmr Php P_r Q i")
    h2pump(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, 0.95, Q_hx= 200e3, pcc=True, filename="test.csv")
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2after")
    # print("qmr Php P_r Q i")
    # h2after(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, 0.95, 0.98)
    #print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(P_r/1000, 3), round(Q/1000, 3), i)
    
    # print("\nh2pre")
    # print("qmr Php Q i")
    # qm_r1, P_hpfp, Q, i = h2pre(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, 0.95)
    # print(round(qm_r1, 4), round(P_hpfp/1000, 3), round(Q/1000, 3), i)







