# -*- coding: utf-8 -*-

import h2flow
import jetaflow
from parallel_comustion import parallel_combustion
import csv
import time
from operator import itemgetter
import os

Verbose = False

# modelling parameters
tolerance = 1e-1
max_iter = 100
rel_fac = 1/200
rel_fac_2 = 0.95

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

eta_pot = 0.74

def save_results(
        filename, arch, t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb0, t0, p0,
        tpr_hx, 
        Q_hx, pcc, v, P_hpfp, P_r, Q, Q_phc, qm_cb, qm_pch, qm_r, qm_v, p_hpfp,
        m_z, i, duration
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
        filewriter.writerow(["qm_cb", "qm_phc/qm_t", "qm_r", "qm_v", "qm", "m_z"])
        filewriter.writerow([qm_cb0, qm_pch, qm_r, qm_v, qm_cb0+qm_pch, m_z])
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
    qm_cb, qm_hpfp = itemgetter("qm_cb0", "qm_hpfp")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    
    # calculate intitial values for independent variables
    h_r = jetaflow.calc_ht(t_cbt, p_cbt, v)
    
    qm_t0 = 0.1
    qm_t = qm_t0
    
    p_hpfp = p_cbt/tpr_fohe + dp_l + dp_inj
    
    # init loop variables
    i = 0
    condition_bool = True
    P_hpfp = 0
    P_lpfp = 0
    
    while condition_bool:
        
        # detect negative mass flow and stop calculation
        if qm_cb + qm_t > qm_hpfp:
            raise ValueError("Mass flow exceeds HPFP limit. QM_T: " + str(qm_t))
            
        # calculate recirculated mass flow
        qm_r = qm_hpfp - qm_t - qm_cb
        
        # init fuel from boost pump
        ff_main = jetaflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        
        # initialise recirculation fuel flow
        ff_r = jetaflow.JetaFlow(qm_r,
            jetaflow.calc_t(h_r, p_lpfp, v), p_lpfp, v
        )
        
        # calculation of lp fuel pump
        P_lpfp, _ = ff_main.pump_hydraulic(p_lpfp, eta_lpfp)
        
        # mix recirculation flow into main flow
        ff_main.mix_flows(ff_r)
        
        # primary heat exchanger
        ff_main.heat_exchanger(Q_fohe, tpr_fohe)
        
        # apply pipe pressure loss
        ff_main.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
        
        # calculation of hp fuel pump
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        
        # split off combustion chamber flow
        ff_cb = ff_main.split_flows(qm_cb)
        
        # apply injector pressure loss
        ff_cb.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_inj)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
        
        # idg heat exchanger
        ff_main.heat_exchanger(Q_idg, 1)

        
        # calculate combustor state
        t_cba = ff_cb.t+ff_cb.v**2/(2*jetaflow.jeta.calc_jeta_cp(ff_cb.t, ff_cb.p))
        p_cba = jetaflow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
        
        
            
        if Verbose:  
            print(qm_t)
            print(i, t_cbt - t_cba, p_cbt - p_cba, h_r - ff_main.ht)
        
        
        # check for convergence
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - ff_main.ht) < tolerance
        )
        
        # advance counter and check for iteration limit
        i+=1
        if i > max_iter:
            print("failed to converge")
            return
        
        # calculate independent variables for next iteration
        h_r = (1+qm_r/qm_cb)*ff_main.ht-qm_r/qm_cb*h_r
        
        p_hpfp_old = p_hpfp
        p_hpfp += (p_cbt - p_cba) * rel_fac_2
        
        qm_t += qm_t0 * (t_cba - t_cbt) * rel_fac * 1.5
        qm_t = max(0, qm_t)
        qm_t = min(qm_t, 0.7)
        
        if p_hpfp > 1.1 * p_hpfp_old :
            p_hpfp = 1.1*p_hpfp_old
        elif p_hpfp < 0.9 * p_hpfp_old:
            p_hpfp = 0.9*p_hpfp_old 
    
    
    
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
        0, p_hpfp, 0, i,
        stop-start
    )
    return


# second jet-a fuel system to determine actual heat demand
def reference2(params, t_cbt, p_cbt, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    tpr_fohe, Q_idg = itemgetter("tpr_fohe", "Q_idg")(params)
    p_lpfp, eta_lpfp, eta_hpfp = itemgetter("p_lpfp", "eta_lpfp", "eta_hpfp")(params)
    qm_cb, qm_hpfp = itemgetter("qm_cb0", "qm_hpfp")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    # no fuel return to tank
    qm_t = 0
    qm_r = qm_hpfp - qm_cb
    
    # calculate intitial values for independent variables
    h_r = jetaflow.calc_ht(t_cbt, p_cbt, v)
    
    p_hpfp = p_cbt/tpr_fohe + dp_l + dp_inj
    
    # init loop variables
    i = 0
    condition_bool = True
    P_hpfp = 0
    P_lpfp = 0
    
    DH = qm_cb*(jetaflow.jeta.calc_jeta_enthalpy(t_cbt, p_cbt)-jetaflow.jeta.calc_jeta_enthalpy(t0, p0))
    
    
    while condition_bool:
        
        # calculate heat demand
        Q_fohe = DH-P_hpfp-P_lpfp-Q_idg
        
        # init fuel from boost pump
        ff_main = jetaflow.JetaFlow(qm_cb+qm_t, t0, p0, v)
        
        # initialise recirculation fuel flow
        ff_r = jetaflow.JetaFlow(qm_r,
            jetaflow.calc_t(h_r, p_lpfp, v), p_lpfp, v
        )
        
        # calculation of lp fuel pump
        P_lpfp, _ = ff_main.pump_hydraulic(p_lpfp, eta_lpfp)
        
        # mix recirculation flow into main flow
        ff_main.mix_flows(ff_r)
        
        # primary heat exchanger
        ff_main.heat_exchanger(Q_fohe, tpr_fohe)
        
        # apply pressure loss
        ff_main.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
        
        # calculation of hp fuel pump
        P_hpfp, _ = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
        
        # split off combustion chamber flow
        ff_cb = ff_main.split_flows(qm_cb)
        
        # apply injector pressure loss
        ff_cb.heat_exchanger(0, (jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_inj)/jetaflow.calc_pt(ff_main.t, ff_main.p, ff_main.v))

        # idg heat exchanger
        ff_main.heat_exchanger(Q_idg, 1)
        
        # calculate combustor state
        t_cba = ff_cb.t+ff_cb.v**2/(2*jetaflow.jeta.calc_jeta_cp(ff_cb.t, ff_cb.p))
        p_cba = jetaflow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
        
        

            
        if Verbose:
            print(i, t_cbt - t_cba, p_cbt - p_cba, h_r - ff_main.ht)
        
        
        # check for convergence
        condition_bool = not (
            abs(t_cbt - t_cba) < tolerance 
            and abs(p_cbt - p_cba) < tolerance 
            and abs(h_r - ff_main.ht) < tolerance
        )
        
        # advance counter and check for iteration limit
        i+=1
        if i > max_iter:
            print("failed to converge")
            return
    
        # calculate independent variables for next iteration
        h_r = (1+qm_r/qm_cb)*ff_main.ht-qm_r/qm_cb*h_r
        
        p_hpfp_old = p_hpfp
        p_hpfp += (p_cbt - p_cba) * rel_fac_2
        
        if p_hpfp > 1.1 * p_hpfp_old :
            p_hpfp = 1.1*p_hpfp_old
        elif p_hpfp < 0.9 * p_hpfp_old:
            p_hpfp = 0.9*p_hpfp_old 
    
    
    # print("saturation margin hpp inlet [bar]: " + str((p_pi - p_sat_pi)/1e5))
    # print("saturation margin injector [bar]: " + str((ff_cb.p - p_sat_cb - dp_inj)/1e5))
   
    stop = time.time()
    
    # filter negative mass flows
    if qm_r < 0 or qm_t < 0:
        raise ValueError("Solution includes negative mass flow")
        
    save_results(
        filename, "reference2", t_cbt, float("nan"), eta_hpfp, eta_lpfp, p_cbt, 
        qm_cb, t0, p0, tpr_fohe, float("nan"), float("nan"), v, 
        P_hpfp, P_lpfp, Q_idg + Q_fohe, 0, qm_cb, qm_t, qm_r,
        0, p_hpfp, 0, i,
        stop-start
    )
    return

# calculate difference across fuel system
def get_dh(qm_cb, t_cb, p_cb, t0, p0, v):
    dH = qm_cb * (h2flow.calc_ht(t_cb, p_cb, 0) - h2flow.calc_ht(t0, p0, v))
    return dH

def h2pump(params, t_cbt, t_hxt, p_cbt, pcc=True, Brewer = False, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    # guess starting values for independent variables
    qm_r0 = qm_cb * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r = qm_r0
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc + dp_inj + dp_l
    
    dH0 = get_dh(qm_cb, t_cbt, p_cbt, t0, p0, v) 
    dH = dH0
    t_r = t_cbt
    p_r = p_cbt
    
    qm_phc = 0
    P_r = 0
    P_hpfp = 0
       
    i = 0
    condition_bool = True
    try:
        while condition_bool:
            
            
            # initialise h2 at engine inlet
            ff_main = h2flow.H2Flow(qm_cb + qm_phc, t0, p0, v, False)
            
            # initialise recirculation h2 flow
            if Brewer:
                ff_r =h2flow.H2Flow(qm_r, t_r, p_hpfp, v, True)
            else:
                ff_r =h2flow.H2Flow(qm_r, t_r, p_r, v, True)
            
            # hpfp calculation
            P_hpfp, t_mfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # calculte recirculation compression
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix recirculation flow into main flow
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # calculate additional heat demand
            Q_phc = dH-P_hpfp-P_r-Q_fohe
            
            # calculate primary heat exchanger
            t_phc = ff_main.heat_exchanger(Q_fohe, tpr_fohe) + dT
            
            # calculate phc heat exchanger
            ff_main.heat_exchanger(Q_phc, tpr_phc)
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
            
            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
        
            # calculate h2 and bleed air requirements of parallel combustion
            if pcc:
                qm_phc, m_z = parallel_combustion(max(0, Q_phc), t_hx=t_phc)
                dH = dH0 * (qm_cb+qm_phc) / qm_cb
            else:
                dH = dH0
                m_z = 0
                qm_phc = 0
                
            # calculate combustor state
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
        
            if Verbose:
                print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, t_r - ff_main.t)
                print(qm_cb, qm_phc)
            
            # check for convergence
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(t_r - ff_main.t) < tolerance
            )
            
            # advance counter and check for iteration limit
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
                
            # advance dependent variables
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            p_r = ff_main.p

            t_r = ff_main.t
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            
            # prevent negative mass flow
            qm_r = max(0, qm_r)
            
            # limit rapid changes of hp pump pressure
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
                
                
                
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_phc < 0:
                print(qm_r, qm_cb, qm_phc)
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, Q_phc + Q_fohe,
                Q_phc, qm_cb, qm_phc, qm_r, 0, p_hpfp, m_z, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "pump", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    return

def h2after(params, t_cbt, t_hxt, p_cbt, pcc=True, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    tpr_vhp, tpr_vlp = itemgetter("tpr_vhp", "tpr_vlp")(params)
    
    start = time.time()
    
    # initial guess for independent variables
    qm_r0 = qm_cb * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v))
    qm_r = qm_r0
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc +dp_inj+dp_l
    
    dH0 = get_dh(qm_cb, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    t_r = t_cbt
    
    p_r = p_cbt
    
    t_phc = 400
    P_r = 0
    P_hpfp = 0
    qm_phc = 0
    
    i = 0
    condition_bool = True
    try:
        while condition_bool:
            
            
            # initialise h2 at engine inlet
            ff_main = h2flow.H2Flow(qm_cb+qm_phc, t0, p0, v, False)
            
            # initialise recirculation flow
            ff_r = h2flow.H2Flow(qm_r, t_r, p_r, v, True)
            
            # vaporiser LP side
            Q_vap = ff_main.heat_to_saturation(tpr_vlp)
            
            # high pressure compressor
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # calculate recirculation compression
            P_r, _ = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix recirculation flow into main h2 flow
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # calculate additional heat demand
            Q_phc = dH-P_hpfp-P_r-Q_fohe
            
            # primary heat exchanger
            ff_main.heat_exchanger(Q_fohe, tpr_fohe)
            
            # vaporiser hp side
            t_phc=ff_main.heat_exchanger(- Q_vap, tpr_vhp) + dT
            
            # phc heat exchanger
            ff_main.heat_exchanger(Q_phc, tpr_phc) 
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
            
            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
            
            # calculate h2 and bleed air requirements of parallel combustion
            if pcc:
                qm_phc, m_z = parallel_combustion(max(0, Q_phc), t_hx=t_phc)
                dH = dH0 * (qm_cb+qm_phc) / qm_cb
            else:
                dH = dH0
                m_z = 0
                qm_phc = 0
            
                
            # calculate combustor state
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
            
            if Verbose:
                print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, t_r - ff_main.t)
                print(qm_cb, qm_phc, qm_r)
                print(dH-Q_fohe-Q_phc-P_hpfp-P_r)
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(t_r - ff_main.t) < tolerance
            )
            
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
                
             # advance dependent variables  
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            t_r = ff_main.t
            p_r = ff_main.p
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_r = max(0.001, qm_r)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
                

        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_phc < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, Q_phc + Q_fohe,
                Q_phc, qm_cb, qm_phc, qm_r, 0, p_hpfp, m_z, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "after", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    return

def h2dual(params, t_cbt, t_hxt, p_cbt, pcc=True, filename="", v=v0, tolerance=tolerance):
    # unpack params dict
    p0, t0 = itemgetter("p0","t0")(params)
    Q_fohe, tpr_fohe, tpr_phc, dT = itemgetter("Q_fohe","tpr_fohe", "tpr_phc", "dT")(params)
    eta_r, eta_hpfp = itemgetter("eta_r", "eta_hpfp")(params)
    qm_cb = itemgetter("qm_cb0")(params)
    dp_l, dp_inj = itemgetter("dp_l", "dp_inj")(params)
    
    start = time.time()
    
    # guess inital values for independent variables
    qm_v0 = qm_cb * (h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(h2flow.sat_t(p0)+10, p0, v))
    qm_v = qm_v0
    
    qm_r0 = qm_cb * (h2flow.calc_ht(t_hxt, p_cbt, v)-h2flow.calc_ht(t0, p0, v))/(h2flow.calc_ht(t_cbt, p_cbt, v)-h2flow.calc_ht(t_hxt, p_cbt, v)) - qm_v
    qm_r = qm_r0
    
    p_hpfp = p_cbt/tpr_fohe/tpr_phc+dp_inj+dp_l
    
    dH0 = get_dh(qm_cb, t_cbt, p_cbt, t0, p0, v)
    dH = dH0
    h_r = dH/qm_cb
    p_r = p_hpfp
    
    qm_phc=0
    P_r = 0
    P_hpfp = 0
    
    i = 0
    condition_bool = True
    try:
        while condition_bool:
        
            # intialise h2 flow
            ff_main = h2flow.H2Flow(qm_cb+qm_phc, t0, p0, v, False)
            
            # initialise recirculation flow (vaporisation)
            ff_v = h2flow.H2Flow(qm_v, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # initialise second recirculation flow
            ff_r = h2flow.H2Flow(qm_r, h2flow.calc_t(h_r, p_r, v, True), p_r, v, True)
            
            # mix vaporisation
            t_va, _ = ff_main.mix_flows(ff_v)
            
            # identify target temperature
            t_vt = h2flow.sat_t(ff_main.p) + 3
            
            # primary compressor
            P_hpfp, t_hpfp = ff_main.pump_hydraulic(p_hpfp, eta_hpfp)
            
            # calculate recirculation compressor
            P_r, t_r = ff_r.pump_hydraulic(p_hpfp, eta_r)
            
            # mix second recirculation
            t_hxa, _ = ff_main.mix_flows(ff_r)
            
            # calculate additional heat demand
            Q_phc = dH-P_hpfp-P_r-Q_fohe
            
            # primary heat exchanger
            t_phc = ff_main.heat_exchanger(Q_fohe, tpr_fohe) + dT
            
            # phc heat exchanger
            ff_main.heat_exchanger(Q_phc, tpr_phc) 
            
            # apply pressure loss
            ff_main.heat_exchanger(0, (h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v)-dp_l)/h2flow.calc_pt(ff_main.t, ff_main.p, ff_main.v))
            
            # split off combustion chamber h2 flow
            ff_cb = ff_main.split_flows(qm_cb)
            
            # apply injector pressure loss
            ff_cb.heat_exchanger(0, (h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)-dp_inj)/h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v))
            
            # calculate h2 and bleed air requirements of parallel combustion
            if pcc:
                qm_phc, m_z = parallel_combustion(max(0, Q_phc), t_hx=t_phc)
                dH = dH0 * (qm_cb+qm_phc) / qm_cb
            else:
                dH = dH0
                m_z = 0
                qm_phc = 0
                
            # calculate combustor state
            p_cba = h2flow.calc_pt(ff_cb.t, ff_cb.p, ff_cb.v)
            t_cba = ff_cb.t + ff_cb.v**2/(2*h2flow.h2.calc_H2_cp(ff_cb.t, ff_cb.p))
             
            if Verbose:
                print(t_cbt - t_cba, t_hxa-t_hxt, p_cbt - p_cba, h_r - ff_main.ht)
                print(dH-Q_phc-Q_fohe-P_r-P_hpfp)
                print(qm_r, qm_v)
            
            condition_bool = not (
                abs(t_cbt - t_cba) < tolerance 
                and abs(t_va-t_vt) < tolerance
                and abs(t_hxa-t_hxt) < tolerance 
                and abs(p_cbt - p_cba) < tolerance 
                and abs(h_r - ff_main.ht) < tolerance)
            
            # advance counter and check for iteration limit
            i+=1
            if i > max_iter:
                raise Exception("Exceeded max iterations")
                
            # advance dependent variables
            p_hpfp_old = p_hpfp
            p_hpfp += (p_cbt - p_cba) * rel_fac_2
            
            p_r = ff_main.p
            h_r = (1+(qm_r+qm_v)/qm_cb)*ff_main.ht-(qm_r+qm_v)/qm_cb*h_r
            
            qm_r += qm_r0*(t_hxt-t_hxa) * rel_fac
            qm_v += qm_v0*(t_vt - t_va) * rel_fac
            qm_r = max(0, qm_r)
            qm_v = max(0, qm_v)
            
            if p_hpfp > 1.1 * p_hpfp_old :
                p_hpfp = 1.1*p_hpfp_old
            elif p_hpfp < 0.9 * p_hpfp_old:
                p_hpfp = 0.9*p_hpfp_old
                
                
                
        stop = time.time()
        if len(filename) > 1:
            if qm_r < 0 or qm_phc < 0 or qm_v < 0:
                raise ValueError("Solution includes negative mass flow")
            save_results(
                filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb,
                t0, p0, tpr_fohe, Q_fohe, pcc, v, P_hpfp, P_r, Q_phc + Q_fohe,
                Q_phc, qm_cb, qm_phc, qm_r, qm_v, p_hpfp, m_z, i,
                stop-start
            )
    except Exception as e:
        print("Failed to converge: " + filename[:-4] + "FAILED" + ".csv")
        print("Number of iterations: " + str(i))
        save_failed(filename, "dual", t_cbt, t_hxt, eta_hpfp, eta_r, p_cbt, qm_cb, t0, p0, tpr_fohe, Q_fohe, pcc, v, e)
        return
    
    return 




if __name__ == "__main__":

    t_bk = 300  
    t_wu = 160
    
    p_bk = 1.33e6
    
    qm_cb = 0.11
    eta_p = 0.154
    eta_v = 0.71
    t0 = 25.2
    p0 = 3.45e5
    Q_fohe = 149e3
    tpr_fohe = 0.95
    tpr_phc = 0.98
    dT = 30
    dp_l = 260e3
    dp_inj = 168.9e3
    ref_params = {
        "p0": 180e3,"t0": 270, "Q_fohe": 112e3,"tpr_fohe": 0.95, "Q_idg": 5e3,
        "eta_lpfp": 0.6, "eta_hpfp": 0.73, "p_lpfp": 930e3, "qm_hpfp": 1.11,
        "qm_cb0": 0.313, "dp_l": 68e3, "dp_inj":300e3
    }
    pump_params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_p, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
    }
    brewer_params = {
        "p0": 344.7e3,"t0": t0, "Q_fohe": 0,"tpr_fohe": 1, 
        "tpr_phc": 1, "dT": 0,"eta_r": 0.71, "eta_hpfp": 0.154, 
        "qm_cb0": 0.166, "dp_l": 31.1e3, "dp_inj":168.9e3+45.5e3
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
    
    folder = "single_results"
    print("reference")
    reference(ref_params, 399, p_bk, filename= os.path.join(folder, "ref.csv"))
    print("reference2")
    reference2(ref_params, 399, p_bk, filename=os.path.join(folder, "ref2.csv"))
    
    # print("\nh2dual")
    # h2dual(dual_params, t_bk, t_wu, p_bk, pcc=True, filename=os.path.join(folder, "dual.csv"))
    
    # print("\nh2pump")
    # h2pump(pump_params, t_bk, t_wu, p_bk, pcc=True, filename=os.path.join(folder, "pump.csv"))
    
    # print("\nh2after")
    # h2after(after_params, t_bk, t_wu, p_bk, pcc=True, filename=os.path.join(folder, "after.csv"))
    
    # print("\nh2pump")
    # h2pump(brewer_params, 264, 200, 1516.2e3, pcc=False, Brewer = True, filename=os.path.join(folder, "brewer.csv"))
    







