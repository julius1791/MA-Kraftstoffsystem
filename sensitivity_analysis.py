# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
import csv
import json
import numpy as np
from pathlib import Path


t_bk = 300
t_wu = 200
eta_p = 0.92
eta_r = 0.9
p_cb = 1.5e6+168.9e3
qm_cb = 0.10998
t0 = 22
p0 = 4.2e5
pcc = True
tpr_hx = 0.95
Q_hx = 200e3
v = 0


# results folder (!CHANGE IN fuelsystem2!!!)
results_dir = "sensitivity"


# subfolder
folders = ["pump", "after", "pre", "dual"]
folder = folders[3]

foldername = os.path.join(os.getcwd(), results_dir, folder)
Path(foldername).mkdir(parents=True, exist_ok=True)

pc_list = []
par_name_list = ["v", "Q_hx", "eta_p", "eta_r", "p_cb", "qm_cb", "tpr_hx"]
par_unit_list = ["[m/s]", "[kW]", "[-]", "[-]", "[kPa]", "[g/s]", "[-]"]
step_size = [2, 20, 0.02, 0.02, 200e3, 0.01, 0.01]

for i in range(len(par_name_list)):
    par_name = par_name_list[i]

    filename = str(i)+ ".csv"
    path = os.path.join(foldername, filename)
    pc = [v, Q_hx, eta_p, eta_r, p_cb, qm_cb, tpr_hx]
    pc[i] = pc[i] + step_size[i]
    pc.append(path)
    pc_list.append(pc)

    
Parallel(n_jobs=7)(
                # t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, tpr_hx, Q_hx=200e3, pcc=pcc, filename=path
    delayed(fuelsystem2.h2dual)(t_bk, t_wu, pc[2], pc[3], pc[4], pc[5], t0, p0, pc[6], Q_hx=pc[1], v=pc[0], pcc=pcc, filename=pc[7]) for pc in pc_list
)
defname = os.path.join(foldername, "default.csv")
fuelsystem2.h2dual(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, tpr_hx, Q_hx=Q_hx, pcc=pcc, filename=defname, v=v)
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
    

pc = [v, Q_hx, eta_p, eta_r, p_cb, qm_cb, tpr_hx]

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
    
par_vals[10] = Q_hx
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
        
    
    
    


