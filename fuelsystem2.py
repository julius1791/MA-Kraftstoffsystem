# -*- coding: utf-8 -*-

import h2flow
import jetaflow
from parallel_comustion import parallel_combustion
import csv
import time
from operator import itemgetter

# modelling parameters
tolerance = 1e-1
max_iter = 100
rel_fac = 1/250
rel_fac_2 = 3/4

# global velocity setting
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
    with open(filename, "w", newline='') as f:
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
        filewriter.writerow(["qm_cb", "qm_phc/qm_t", "qm_r", "qm_v", "qm"])
        filewriter.writerow([qm_cb0, qm_pch, qm_r, qm_v, qm_cb0+qm_pch])
        filewriter.writerow([
            "number of iterations", i, "Execution time: ", duration
            ])
    return

def save_failed(
        filename, arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v, exception
    ):
    if len(filename) > 1:
        failed = filename[:-4] + "FAILED" + ".csv"
        with open(failed, "w", newline='') as f:
            filewriter = csv.writer(f)
            filewriter.writerow(["architecture", "t_cbt", "t_hxt", "eta_hpfp", "eta_r", "p_cbt", "qm_cb0", "t0", "p0", "tpr_hx", "Q_hx", "pcc", "v"])
            filewriter.writerow([arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_hx, Q_hx, pcc, v])
            filewriter.writerow(["FAILED TO CONVERGE"])
            filewriter.writerow([exception])
    return

def reference(params, t_cbt, p_cbt, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, Q_idg = itemgetter("Q_fohe","tpr_fohe", "Q_idg")(params)
    p_lpfp, eta_lpfp, eta_hpfp = itemgetter("p_lpfp", "eta_lpfp", "eta_hpfp")(params)
    qm_cb0, qm_hpfp = itemgetter("qm_cb0", "qm_hpfp")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    # correct combustion chamber fuel mass flow for temperature
    lhv_tcbt = (
        lhv_jeta_288 
        + jetaflow.jeta.calc_jeta_enthalpy(t_cbt, p_cbt) 
        - jetaflow.jeta.calc_jeta_enthalpy(288, p_cbt)
    )
    qm_cb = qm_cb0 * (lhv_jeta_288 / lhv_tcbt)
    
    # calculate intitial values for independent variables
    h_r = jetaflow.calc_ht(t_cbt, p_cbt, v)
    
    qm_t = 0.1
    qm_t0 = qm_t
    
    p_hpfp = p_cbt/tpr_fohe + dp_l + dp_inj
    
    # init loop variables
    i = 0
    condition_bool = True
    
    while condition_bool:
        # fuel from boost pump
        ff_main = jetaflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        
        # calculation of lp fuel pump
        P_lpfp, _ = ff_main.pump_hydraulic(p_lpfp, eta_lpfp)
        
        # detect negative mass flow and stop calculation
        if qm_cb + qm_t > qm_hpfp:
            raise ValueError("Mass flow exceeds HPFP limit. QM_T: " + str(qm_t))
            
        # initialise recirculation fuel flow
        ff_r = jetaflow.JetaFlow(qm_hpfp-qm_t-qm_cb,
            jetaflow.calc_t(h_r, ff_main.p, v), ff_main.p, v
        )
        
        # mix recirculation flow into main flow
        ff_main.mix_flows(ff_r)
        
        # primary heat exchanger
        ff_main.heat_exchanger(Q_fohe, tpr_fohe)
        
        # saturation pressure and static pressure at hp pump inlet
        # p_sat_pi = jetaflow.sat_p(ff_main.t)
        # p_pi = ff_main.p
        
        # calculation of hp fuel pump
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        
        # apply pressure loss
        ff_main.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
        
        # split off combustion chamber flow
        ff_cb = ff_main.split_flows(qm_cb)
        
        # apply injector pressure loss
        ff_cb.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_inj)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
        
        # vapour pressure at injector
        # p_sat_cb = jetaflow.sat_p(ff_cb.t)
        
        # stagnation pressure at injector
        t_cba = ff_cb.t+ff_cb.v**2/(2*jetaflow.jeta.calc_jeta_cp(ff_cb.t, ff_cb.p))
        p_cba = jetaflow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
        
        # idg heat exchanger
        ff_main.heat_exchanger(Q_idg, 1)
        
        # calculate recirculated mass flow
        qm_r = qm_hpfp - qm_t - qm_cb
        
        # calculate independent variables for next iteration
        h_rold = h_r
        h_r = (1+qm_r/qm_cb)*ff_main.ht-qm_r/qm_cb*h_r
        
        p_hpfp_old = p_hpfp
        p_hpfp += (p_cbt - p_cba) * rel_fac_2
        
        qm_t += qm_t0 * (t_cba - t_cbt) * rel_fac * 4
        qm_t = max(0, qm_t)
        qm_t = min(qm_t, 0.7)
        
        if p_hpfp > 1.1 * p_hpfp_old :
            p_hpfp = 1.1*p_hpfp_old
        elif p_hpfp < 0.9 * p_hpfp_old:
            p_hpfp = 0.9*p_hpfp_old 
            
            
        
        # print(i, t_cbt - t_cba, p_cbt - p_cba, h_r - h_rold)
        
        
        # check for convergence
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - h_rold) < tolerance
        )
        
        # advance counter and check for iteration limit
        i+=1
        if i > max_iter:
            print("failed to converge")
            return
    
    
    
    # print("saturation margin hpp inlet [bar]: " + str((p_pi - p_sat_pi)/1e5))
    # print("saturation margin injector [bar]: " + str((ff_cb.p - p_sat_cb - dp_inj)/1e5))
   
    stop = time.time()
    
    # filter negative mass flows
    if qm_r < 0 or qm_t < 0:
        raise ValueError("Solution includes negative mass flow")
        
    save_results(
        filename, "reference", t_cbt, float("nan"), eta_hpfp, eta_lpfp, p_cbt, 
        qm_cb, t0, p0, tpr_fohe, float("nan"), float("nan"), v, 
        P_hpfp, P_lpfp, Q_idg + Q_fohe, 0, qm_cb, qm_t, qm_r,
        0, p_hpfp, i,
        stop-start
    )
    return

# calculate difference across fuel system
def get_dh(qm_cb, t_cb, p_cb, t0, p0, v):
    dH = qm_cb * (h2flow.calc_ht(t_cb, p_cb, 0) - h2flow.calc_ht(t0, p0, v))
    return dH

def h2pump(params, t_cbt, t_hxt, p_cbt, pcc=True, corr=True, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb0 = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    if corr:
        # correct combustion chamber fuel mass flow for temperature
        qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    
    # guess starting values for independent variables
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc + dp_inj + dp_l
    
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v) 
    dH = dH0
    h_r = dH/qm_cb0
    p_r = p_cbt
    t_cba = t_cbt
    
    t_phc = 400
    
    P_r = 0
    P_hpfp = 0
       
    i = 0
    condition_bool = True
    try:
        while condition_bool:
            
            # calculate h2 requirements of parallel combustion
            if pcc:
                qm_cb = qm_cb0 + parallel_combustion(max(0, dH-P_r-P_hpfp-Q_fohe), t_cbt, t_hx=t_phc+dT)
                # if i > 10:
                #     dH0 += (t_cbt - t_cba)/1000
                
                dH = dH0 * qm_cb / qm_cb0
                
            else:
                qm_cb = qm_cb0
                dH = dH0
            
            # initialise h2 at engine inlet
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            
            # hpfp calculation
            P_hpfp, t_mfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # initialise recirculation h2 flow
            ff_r =h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # calculte recirculation compression
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix recirculation flow into main flow
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # calculate primary heat exchanger
            t_phc = ff_main.heat_exchanger(Q_fohe, tpr_fohe)
            
            # calculate phc heat exchanger
            ff_main.heat_exchanger(dH-P_hpfp-P_r-Q_fohe, tpr_phc)
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))

            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
        
            
            # advance independent variables
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            p_r = ff_main.p
            if not corr:
                p_r = p_hpfp
            h_rold = h_r
            h_r = (1+qm_r/qm_cb)*ff_main.ht-qm_r/qm_cb*h_r
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
            
           
            # print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
            
            # check for convergence
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
            )
            
            # advance counter and check for iteration limit
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0:
                print(qm_r, qm_cb, qm_cb0)
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_fohe, qm_cb, qm_cb-qm_cb0, qm_r, 0, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    return

def h2after(params, t_cbt, t_hxt, p_cbt, pcc=True, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb0 = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    tpr_vhp, tpr_vlp = itemgetter("tpr_vhp", "tpr_vlp")(params)
    
    start = time.time()
    
    # correct combustion chamber fuel mass flow for temperature
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    
    # initial guess for independent variables
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r0 = qm_r
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc +dp_inj+dp_l
    
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    
    p_r = p_cbt
    
    t_phc = 400
    P_r = 0
    P_hpfp = 0
    
    i = 0
    condition_bool = True
    try:
        while condition_bool:
            # calculate h2 requirements of parallel combustion
            if pcc:
                qm_cb = qm_cb0 + parallel_combustion(max(0, dH-P_r-P_hpfp-Q_fohe), t_cbt, t_hx=t_phc+dT)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
                dH=dH0
            
            # initialise h2 at engine inlet
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            
            # vaporiser LP side
            Q_vap = ff_main.heat_to_saturation(tpr_vlp)
            
            # high pressure compressor
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # initialise recirculation flow
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # calculate recirculation compression
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix recirculation flow into main h2 flow
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # primary heat exchanger
            ff_main.heat_exchanger(Q_fohe, tpr_fohe)
            
            # vaporiser hp side
            t_phc=ff_main.heat_exchanger(- Q_vap, tpr_vhp)
            
            # phc heat exchanger
            ff_main.heat_exchanger(dH-P_hpfp-P_r-Q_fohe, tpr_phc) 
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
            
            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
            
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            h_rold = h_r
            h_r = (1+qm_r/qm_cb)*ff_main.ht-qm_r/qm_cb*h_r
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0.001, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
                
            # print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance
            )
            
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_fohe, qm_cb, qm_cb-qm_cb0, qm_r, 0, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    return

def h2dual(params, t_cbt, t_hxt, p_cbt, pcc=True, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb0 = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    # correct combustion chamber fuel flow for temperature
    qm_cb0 = qm_cb0 * lhv_h2_200 / (lhv_h2_200 - h2flow.h2.calc_H2_enthalpy(200, 1e6) + h2flow.h2.calc_H2_enthalpy(t_cbt, 1e6))
    
    # guess inital values for independent variables
    qm_v = qm_cb0 * (h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v))
    qm_v0 = qm_v
    
    qm_r = qm_cb0 * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v)) - qm_v
    qm_r0 = qm_r
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc+dp_inj+dp_l
    
    dH0 = get_dh(qm_cb0, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb0
    t_phc = 400
    p_r = p_hpfp
    
    P_r = 0
    P_hpfp = 0
    
    i = 0
    condition_bool = True
    try:
        while condition_bool:
            
            if pcc:
                qm_cb = qm_cb0 + parallel_combustion(max(0, dH-P_r-P_hpfp-Q_fohe), t_cbt, t_hx=t_phc+dT)
                dH = dH0 * qm_cb / qm_cb0
            else:
                qm_cb = qm_cb0
                dH = dH0
            # intialise h2 flow
            ff_main = h2flow.H2Flow(qm_cb, t0, p0, v, False)
            
            # initialise recirculation flow (vaporisation)
            ff_v = h2flow.H2Flow(qm_v, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # mix vaporisation
            t_va, _ = ff_main.mix_flows(ff_v)
            
            # identify target temperature
            t_vt = h2flow.sat_t(ff_main.p) + 3
            
            # primary compressor
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # initialise second recirculation flow
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # calculate recirculation compressor
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix second recirculation
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # primary heat exchanger
            t_phc = ff_main.heat_exchanger(Q_fohe, tpr_fohe)
            
            # phc heat exchanger
            ff_main.heat_exchanger(dH-P_hpfp-P_r-Q_fohe, tpr_phc) 
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
            
            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
            
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            p_r = ff_main.p
            h_r = (1+(qm_r+qm_v)/qm_cb)*ff_main.ht-(qm_r+qm_v)/qm_cb*h_r
            h_rold = h_r
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_v += qm_v0*(t_vt - t_va) * rel_fac
            qm_r = max(0, qm_r)
            qm_v = max(0, qm_v)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
                
            # print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - h_rold)
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_va-t_vt) < tolerance
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - h_rold) < tolerance)
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_cb - qm_cb0 < 0 or qm_v < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, dH-P_hpfp-P_r,
                dH-P_hpfp-P_r-Q_fohe, qm_cb, qm_cb-qm_cb0, qm_r, qm_v, p_hpfp, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    
    return 




if __name__ == "__main__":

    t_bk = 600  
    t_wu = 100
    
    p_bk = 1.33e6
    
    qm_cb = 0.10998
    eta_p = 0.154
    eta_v = 0.71
    t0 = 25.2
    p0 = 3.45e5
    Q_fohe = 159e3
    tpr_fohe = 0.95
    tpr_phc = 0.98
    dT = 20
    dp_l = 260e3
    dp_inj = 168.9e3
    ref_params = {
        "p0": 180e3,"t0": 270, "Q_fohe": 122e3,"tpr_fohe": 0.95, "Q_idg": 5e3,
        "eta_lpfp": 0.6, "eta_hpfp": 0.73, "p_lpfp": 930e3, "qm_hpfp": 1.113,
        "qm_cb0": 0.31305, "dp_l": 68e3, "dp_inj":300e3
    }
    pump_params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_p, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
    }
    brewer_params = {
        "p0": 344.7e3,"t0": t0, "Q_fohe": 0,"tpr_fohe": 1, 
        "tpr_phc": 1, "dT": 0,"eta_r": 0.71, "eta_hpfp": 0.154, 
        "qm_cb0": 0.166, "dp_l": 30e3, "dp_inj":168.9e3+45.5e3
    }
    after_params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_v, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj, "tpr_vhp":1, 
        "tpr_vlp": 1
    }
    dual_params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_v, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
    }
    
    
    # print("reference")
    # reference(ref_params, 399.15, p_bk, filename="ref.csv")
    
    # print("\nh2dual")
    # h2dual(dual_params, t_bk, t_wu, p_bk, pcc=True, filename="dual.csv")
    
    # print("\nh2pump")
    # h2pump(pump_params, t_bk, t_wu, p_bk, pcc=True, filename="pump.csv")
    
    # print("\nh2after")
    # h2after(after_params, t_bk, t_wu, p_bk, pcc=True, filename="after.csv")
    
    print("\nh2pump")
    h2pump(brewer_params, 264, 200, 1516.2e3, pcc=False, corr=False, filename="brewer.csv")
    







