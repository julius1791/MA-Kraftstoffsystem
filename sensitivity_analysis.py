# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
import csv
import json
import numpy as np
from pathlib import Path

func = "dual"

# results folder (!CHANGE IN fuelsystem2!!!)
results_dir = "sensitivity"

# subfolder
folders = ["pump", "after", "dual"]


t_bk = 270
t_wu = 220
p_bk = 1.33e6

qm_cb = 0.10998
qm_cbs = qm_cb+0.005
eta_p = 0.154
eta_ps = 0.146
eta_v = 0.71
eta_vs = 0.67
t0 = 25.2
p0 = 4.2e5
Q_fohe = 159e3
Q_fohes = Q_fohe*0.95
tpr_fohe = 0.95
tpr_fohes = 1-(1-tpr_fohe)*1.05
tpr_phc = 0.98
tpr_phcs = 1-(1-tpr_phc)*1.05
dT = 20
dp_l = 260e3
dp_ls = dp_l*1.05
dp_inj = 168.9e3
dp_injs = dp_inj*1.05
tpr_vhps = 0.98
tpr_vlps = 0.99




if func == "pump":
    folder = folders[0]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj"]
    par_unit_list = ["kW", "-", "-", "-", "-", "kg/s", "kPa", "kPa"]
    par_value_list = [Q_fohes, tpr_fohes, tpr_phcs, eta_vs, eta_ps, qm_cbs, dp_ls, dp_injs]
    params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_p, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
    }
elif func == "after":
    folder = folders[1]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj", "tpr_vhp", "tpr_vlp"]
    par_unit_list = ["kW", "-", "-", "-", "-", "kg/s", "kPa", "kPa", "-", "-"]
    par_value_list = [Q_fohes, tpr_fohes, tpr_phcs, eta_vs, eta_vs, qm_cbs, dp_ls, dp_injs, tpr_vhps, tpr_vlps]
    params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_v, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj, "tpr_vhp":1, 
        "tpr_vlp": 1
    }
elif func == "dual":
    folder = folders[2]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj"]
    par_unit_list = ["kW", "-", "-", "-", "-", "kg/s", "kPa", "kPa"]
    par_value_list = [Q_fohes, tpr_fohes, tpr_phcs, eta_vs, eta_vs, qm_cbs, dp_ls, dp_injs]
    params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_v, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj
    }



foldername = os.path.join(os.getcwd(), results_dir, folder)
Path(foldername).mkdir(parents=True, exist_ok=True)

param_comb_list = []
path_list = []
step_size = []
defname = os.path.join(foldername, "default.csv")

for i in range(len(par_name_list)):

    filename = str(i)+ ".csv"
    path = os.path.join(foldername, filename)
    param_comb = params.copy()
    param_comb[par_name_list[i]] = par_value_list[i]
    param_comb_list.append(param_comb)
    path_list.append(path)
    step_size.append(par_value_list[i]-params[par_name_list[i]])

if func == "pump":
    Parallel(n_jobs=i)(delayed(fuelsystem2.h2pump)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    fuelsystem2.h2pump(params, t_bk, t_wu, p_bk, filename=defname)
elif func == "after":
    Parallel(n_jobs=i)(delayed(fuelsystem2.h2after)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    fuelsystem2.h2after(params, t_bk, t_wu, p_bk, filename=defname)
elif func == "dual":
    Parallel(n_jobs=i)(delayed(fuelsystem2.h2dual)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    fuelsystem2.h2dual(params, t_bk, t_wu, p_bk, filename=defname)

with open(defname, newline='') as def_file:
    data_def = csv.reader(def_file)
    _ = next(data_def)
    _ = next(data_def)
    _ = next(data_def)
    P_row = next(data_def)
    _ = next(data_def)
    Q_row = next(data_def)
    P_m0 = float(P_row[0])
    Q0 = float(Q_row[0])
    P_r0 = float(P_row[1])
    

jacobian = np.zeros([len(par_name_list), 3])
absolute_difference = np.zeros([len(par_name_list), 3])
relative_difference = np.zeros([len(par_name_list), 3])
for i in range(len(par_name_list)):
    filename = str(i)+ ".csv"
    path = os.path.join(foldername, filename)
    with open(path, newline="") as f:
        data_f = csv.reader(f)
        par_names = next(data_f)
        par_vals = next(data_f)
        _ = next(data_f)
        P_row = next(data_f)
        _ = next(data_f)
        Q_row = next(data_f)

        P_mu = float(P_row[0])
        Qu = float(Q_row[0])
        P_ru = float(P_row[1])
        
    dQ = Qu-Q0
    dP_m = P_mu - P_m0
    dP_r = P_ru - P_r0
    if step_size[i] > 1e3:
        step_size[i] = step_size[i]/1e3
    gQ = dQ/step_size[i]
    gP_m = dP_m/step_size[i]
    gP_r = dP_r/step_size[i]
    jacobian[i, 0] = gP_m
    jacobian[i, 1] = gP_r
    jacobian[i, 2] = gQ
    absolute_difference[i, 0] = dP_m
    absolute_difference[i, 1] = dP_r
    absolute_difference[i, 2] = dQ
    relative_difference[i, 0] = dP_m/P_m0
    try:
        relative_difference[i, 1] = dP_r/P_r0
    except:
        #pre init
        pass
    relative_difference[i, 2] = dQ/Q0
    
params = dict()
for i, name, value in zip(range(len(par_names)), par_names, par_vals):
    params.update({name: value})
content = {"Parameters": params, "Absolute Differences": absolute_difference.tolist(), "Jacobian": jacobian.tolist(), "Parameter step": step_size, "Parameter names": par_name_list}
jsonfn = os.path.join(results_dir, folder + ".json")
csvfn = os.path.join(results_dir, folder + ".csv")
with open(jsonfn, "w") as jsonout:
    jsonout.write(json.dumps(content, indent=4))
with open(csvfn, 'w', newline="") as csvout:
    csvwriter = csv.writer(csvout)
    csvwriter.writerow(["Variable", "Einheit", "Schrittweite", "dP_mfp/P_mfp [%]", "dP_r/P_r [%]", "dQ/Q [%]"])
    for i in range(len(par_name_list)):
        csvwriter.writerow([par_name_list[i], par_unit_list[i], step_size[i], relative_difference[i, 0]*1e2, relative_difference[i, 1]*1e2, relative_difference[i, 2]*1e2])
        
    
    
    


