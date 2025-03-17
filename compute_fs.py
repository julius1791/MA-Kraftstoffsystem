# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
import time
from gather_data import get_Data
from pathlib import Path


func = "dual"

# results folder
results_dir = "results"

# subfolder
folders = ["pump", "after", "dual"]


t_bk_list = range(150, 410, 10)
t_wu_list = range(100, 290, 10)

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

pump_params = {
    "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
    "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_p, 
    "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
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

def create_param_combs(t_bk_list, t_wu_list, folder, results_dir):

    param_combinations = []
    
    foldername = os.path.join(os.getcwd(), results_dir, folder)
    Path(foldername).mkdir(parents=True, exist_ok=True)
    
    for t_bk in t_bk_list:
        for t_wu in t_wu_list:
            if t_wu + 20 > t_bk:
                continue
            filename = str(t_bk) + str(t_wu) + ".csv"
            path = os.path.join(results_dir, folder, filename)
            param_combinations.append([t_bk, t_wu, path])
    return param_combinations



start = time.time()
if func == "pump":
    folder = folders[0]
    param_combinations = create_param_combs(t_bk_list, t_wu_list, folders[0], results_dir)
    Parallel(n_jobs=20)(
        delayed(fuelsystem2.h2pump)(pump_params, t_bk, t_wu, p_bk, filename=path) for t_bk, t_wu, path in param_combinations
    )
elif func == "after":
    folder = folders[1]
    param_combinations = create_param_combs(t_bk_list, t_wu_list, folders[1], results_dir)
    Parallel(n_jobs=20)(
        delayed(fuelsystem2.h2after)(after_params, t_bk, t_wu, p_bk, filename=path) for t_bk, t_wu, path in param_combinations
    )
elif func == "dual":
    folder = folders[2]
    param_combinations = create_param_combs(t_bk_list, t_wu_list, folders[2], results_dir)
    Parallel(n_jobs=20)(
        delayed(fuelsystem2.h2dual)(dual_params, t_bk, t_wu, p_bk, filename=path) for t_bk, t_wu, path in param_combinations
    )

        
stop = time.time()
get_Data(results_dir, folder)
print(stop-start)
