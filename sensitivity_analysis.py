# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
from gather_data import get_Data
import csv
import json
import numpy as np


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

step_size = 1.01

# results folder (!CHANGE IN fuelsystem2!!!)
results_dir = "sensitivity"

# subfolder
folders = ["pump", "after", "pre", "dual"]
folder = folders[3]
pc_list = []
par_name_list = ["t_bk", "t_wu", "eta_p", "eta_r", "p_cb", "qm_cb", "tpr_hx", "Q_hx"]

for i in range(0, 8):
    par_name = par_name_list[i]

    filename = str(i)+ "u.csv"
    path = os.path.join(folder, filename)
    pc = [t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, tpr_hx, Q_hx]
    pc[i] = pc[i]*step_size
    pc.append(path)
    pc_list.append(pc)
    
    filename = str(i) + "l.csv"
    path = os.path.join(folder, filename)
    pc = [t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, tpr_hx, Q_hx]
    pc[i] = pc[i]/step_size
    pc.append(path)
    pc_list.append(pc)
    
Parallel(n_jobs=16)(
                # t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, tpr_hx, Q_hx=200e3, pcc=pcc, filename=path
    delayed(fuelsystem2.h2dual)(pc[0], pc[1], pc[2], pc[3], pc[4], pc[5], t0, p0, pc[6], Q_hx=pc[7], pcc=pcc, filename=pc[8]) for pc in pc_list
)

pc = [t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, tpr_hx, Q_hx]

jacobian = np.zeros([len(par_name_list), 3])
absolute_difference = np.zeros([len(par_name_list), 3])
for i in range(0,8):
    filename = str(i)+ "u.csv"
    path = os.path.join(results_dir, folder, filename)
    with open(path, newline="") as upper:
        data_f = csv.reader(upper)
        par_names = next(data_f)
        par_vals = next(data_f)
        _ = next(data_f)
        P_row = next(data_f)
        _ = next(data_f)
        Q_row = next(data_f)

        P_mu = float(P_row[0])
        Qu = float(Q_row[0])
        P_ru = float(P_row[1])
        
    filename = str(i)+ "l.csv"
    path = os.path.join(results_dir, folder, filename)
    with open(path, newline="") as lower:
        data_f = csv.reader(lower)
        par_names = next(data_f)
        par_vals = next(data_f)
        _ = next(data_f)
        P_row = next(data_f)
        _ = next(data_f)
        Q_row = next(data_f)

        P_ml = float(P_row[0])
        Ql = float(Q_row[0])
        P_rl = float(P_row[1])
    dQ = Qu-Ql
    dP_m = P_mu - P_ml
    dP_r = P_ru - P_rl
    par_upper = pc[i]*step_size
    par_lower = pc[i]/step_size
    dpar = par_upper - par_lower
    gQ = dQ/dpar
    gP_m = dP_m/dpar
    gP_r = dP_r/dpar
    jacobian[i, 0] = gP_m
    jacobian[i, 1] = gP_r
    jacobian[i, 2] = gQ
    absolute_difference[i, 0] = dP_m
    absolute_difference[i, 1] = dP_r
    absolute_difference[i, 2] = dQ
    
par_vals[10] = Q_hx
params = dict()
for i, name, value in zip(range(len(par_names)), par_names, par_vals):
    if i == 6:
        continue
    params.update({name: value})
content = {"Parameters": params, "Absolute Differences": absolute_difference.tolist(), "Jacobian": jacobian.tolist()}
jsonfn = os.path.join(results_dir, folder + ".json")
with open(jsonfn, "w") as jsonout:
    jsonout.write(json.dumps(content, indent=4))
    
    


