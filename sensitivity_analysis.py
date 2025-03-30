# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import fuelsystem
import os
from joblib import Parallel, delayed
import csv
import json
import numpy as np
from pathlib import Path
import math

func = "reference"

# results folder (!CHANGE IN fuelsystem2!!!)
results_dir = "sensitivity"

# subfolder
folders = ["pump", "after", "dual", "reference"]
fac = 1.1

t_bk = 300
t_wu = 160
p_bk = 1.33e6
t_bk_jeta = 399.15

qm_cb = 0.11
qm_cbs = qm_cb
eta_p = 0.154
eta_ps = eta_p-0.008
eta_v = 0.71
eta_vs = eta_v-0.04
t0 = 25.2
p0 = 3.45e5
Q_fohe = 149e3
Q_fohes = Q_fohe*0.95
tpr_fohe = 0.95
tpr_fohes = tpr_fohe-0.01
tpr_phc = 0.98
tpr_phcs = tpr_phc-0.004
dT = 20
dp_l = 260e3
dp_ls = dp_l*1.2
dp_inj = 169e3
dp_injs = dp_inj*1.05
tpr_vhps = 0.995
tpr_vlps = 0.99


def sigfig(num, sigfigs):
    """round number to specified significant figures"""
    # filter math range error if num == 0
    if num == float("inf"):
        return num
    if num != 0:
        power = math.ceil(math.log(abs(num), 10))
        rounded = round(num/10**power, sigfigs)*10**power
    else:
        rounded = 0
    return rounded

def fixfloat(string):
    if string == "+inf":
        return string
    sep_index = string.find(".")
    string = string[:sep_index] + string[sep_index+1:]
    if "99999" in string:
        index = string.find("99999")
        string = string[:index-1] + str(int(string[index-1])+1)
    if "00000" in string:
        index = string.find("00000")
        string = string[:index]
    if len(string) < sep_index+1:
        while len(string) < sep_index:
            string += "0"
    else:
        string = string[:sep_index] + "." + string[sep_index:]
    while len(string) < 5:
        string += "0"
    if string[1] == "0" and len(string) == 5 and string[2] == ".":
        string += "0"
    return string



if func == "reference":
    folder = folders[3]
    par_hr_list = ["$\dot{Q}_\mathrm{FOHE}$", "$1-\pi_\mathrm{FOHE}$", "$p_\mathrm{LPFP}$", "$\eta_\mathrm{LPFP}$", "$\eta_\mathrm{HPFP}$", "$\Delta p_\mathrm{L}$", "$\Delta p_\mathrm{inj}$"]
    par_name_list = ["Q_fohe", "tpr_fohe", "p_lpfp", "eta_lpfp", "eta_hpfp", "dp_l", "dp_inj"]
    par_unit_list = ["\si{\kilo\W}", "-", "-", "-", "-", "\si{\kilo\Pa}", "\si{\kilo\Pa}"]
    par_value_list = [112e3-117e3*0.05, tpr_fohes, 930e3*0.9, 0.6-0.03, 0.73-0.04, 68e3*1.2, 300e3*1.1]
    params = {
        "p0": 180e3,"t0": 270, "Q_fohe": 122e3,"tpr_fohe": 0.95, "Q_idg": 5e3,
        "eta_lpfp": 0.6, "eta_hpfp": 0.73, "p_lpfp": 930e3, "qm_hpfp": 1.11,
        "qm_cb0": 0.313, "dp_l": 68e3, "dp_inj":300e3
    }
if func == "pump":
    folder = folders[0]
    par_hr_list = ["$\dot{Q}_\mathrm{FOHE}$", "$1-\pi_\mathrm{FOHE}$", "$1-\pi_\mathrm{PHC}$", "$\eta_\mathrm{RV}$", "$\eta_\mathrm{HPFP}$", "$\dot{m}_\mathrm{BK,0}$", "$\Delta p_\mathrm{L}$", "$\Delta p_\mathrm{inj}$"]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj"]
    par_unit_list = ["\si{\kilo\W}", "-", "-", "-", "-", "\si{\kg\per\s}", "\si{\kilo\Pa}", "\si{\kilo\Pa}"]
    par_value_list = [Q_fohes, tpr_fohes, tpr_phcs, eta_vs, eta_ps, qm_cbs, dp_ls, dp_injs]
    params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_p, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj": dp_inj
    }
elif func == "after":
    folder = folders[1]
    par_hr_list = ["$\dot{Q}_\mathrm{FOHE}$", "$1-\pi_\mathrm{FOHE}$", "$1-\pi_\mathrm{PHC}$", "$\eta_\mathrm{RV}$", "$\eta_\mathrm{HPFC}$", "$\dot{m}_\mathrm{BK,0}$", "$\Delta p_\mathrm{L}$", "$\Delta p_\mathrm{inj}$", "$1-\pi_\mathrm{V, HP}$", "$1-\pi_\mathrm{V, LP}$"]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj", "tpr_vhp", "tpr_vlp"]
    par_unit_list = ["\si{\kilo\W}", "-", "-", "-", "-", "\si{\kg\per\s}", "\si{\kilo\Pa}", "\si{\kilo\Pa}", "-", "-"]
    par_value_list = [Q_fohes, tpr_fohes, tpr_phcs, eta_vs, eta_vs, qm_cbs, dp_ls, dp_injs, tpr_vhps, tpr_vlps]
    params = {
        "p0": p0,"t0": t0, "Q_fohe": Q_fohe,"tpr_fohe": tpr_fohe, 
        "tpr_phc": tpr_phc, "dT": dT,"eta_r": eta_v, "eta_hpfp": eta_v, 
        "qm_cb0": qm_cb, "dp_l": dp_l, "dp_inj":dp_inj, "tpr_vhp":1, 
        "tpr_vlp": 1
    }
elif func == "dual":
    folder = folders[2]
    par_hr_list = ["$\dot{Q}_\mathrm{FOHE}$", "$1-\pi_\mathrm{FOHE}$", "$1-\pi_\mathrm{PHC}$", "$\eta_\mathrm{RV}$", "$\eta_\mathrm{HPFC}$", "$\dot{m}_\mathrm{BK,0}$", "$\Delta p_\mathrm{L}$", "$\Delta p_\mathrm{inj}$"]
    par_name_list = ["Q_fohe", "tpr_fohe", "tpr_phc", "eta_r", "eta_hpfp", "qm_cb0", "dp_l", "dp_inj"]
    par_unit_list = ["\si{\kilo\W}", "-", "-", "-", "-", "\si{\kg\per\s}", "\si{\kilo\Pa}", "\si{\kilo\Pa}"]
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
rel_step_size = []
defname = os.path.join(foldername, "default.csv")

for i in range(len(par_name_list)):

    filename = str(i)+ ".csv"
    path = os.path.join(foldername, filename)
    param_comb = params.copy()
    param_comb[par_name_list[i]] = par_value_list[i]
    param_comb_list.append(param_comb)
    path_list.append(path)
    if not "pi" in par_hr_list[i]:
        step_size.append(par_value_list[i]-params[par_name_list[i]])
        rel_step_size.append(step_size[i]/params[par_name_list[i]])
        print(par_name_list[i])
    elif not "V," in par_hr_list[i]:
        step_size.append(-par_value_list[i]+params[par_name_list[i]])
        rel_step_size.append(step_size[i]/(1-params[par_name_list[i]]))
    else:
        step_size.append(-par_value_list[i]+params[par_name_list[i]])
        rel_step_size.append(float("inf"))
param_comb_list.append(params)
path_list.append(defname)


if False:
    if func == "reference":
        Parallel(n_jobs=i+1)(delayed(fuelsystem.reference)(pc, t_bk_jeta, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    if func == "pump":
        Parallel(n_jobs=i+1)(delayed(fuelsystem.h2pump)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    elif func == "after":
        Parallel(n_jobs=i+1)(delayed(fuelsystem.h2after)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))
    elif func == "dual":
        Parallel(n_jobs=i+1)(delayed(fuelsystem.h2dual)(pc, t_bk, t_wu, p_bk, filename=path) for pc, path in zip(param_comb_list, path_list))

with open(defname, newline='') as def_file:
    data_def = csv.reader(def_file)
    _ = next(data_def)
    _ = next(data_def)
    _ = next(data_def)
    P_row = next(data_def)
    _ = next(data_def)
    Q_row = next(data_def)
P_m0 = float(P_row[0])
Q0 = float(Q_row[1])
if Q0 == 0:
    Q0 = 1e3
P_r0 = float(P_row[1])
    

jacobian = np.zeros([len(par_name_list), 3])
absolute_difference = np.zeros([len(par_name_list), 3])
relative_difference = np.zeros([len(par_name_list), 3])
deltaHP = []
deltaRV = []
deltaQ = []
deltaP = []
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
    Qu = float(Q_row[1])
    P_ru = float(P_row[1])
        
    dQ = Qu-Q0
    dP_m = P_mu - P_m0
    dP_r = P_ru - P_r0
    
    if abs(step_size[i]) > 1e3:
        step_size[i] = step_size[i]/1e3
    # gQ = dQ/step_size[i]
    # gP_m = dP_m/step_size[i]
    # gP_r = dP_r/step_size[i]
    # jacobian[i, 0] = gP_m
    # jacobian[i, 1] = gP_r
    # jacobian[i, 2] = gQ
    absolute_difference[i, 0] = dP_m
    absolute_difference[i, 1] = dP_r
    absolute_difference[i, 2] = dQ
    relative_difference[i, 0] = dP_m/P_m0
    # print(par_name_list[i], dP_m, dP_r, P_m0, P_r0)
    # print(round(dP_m/P_m0*100, 3), round(dP_r/P_r0*100, 3))
    deltaHP.append(sigfig(dP_m/P_m0*100, 3))
    deltaRV.append(sigfig(dP_r/P_r0*100, 3))
    deltaQ.append(sigfig(dQ/Q0*100, 3))
    deltaP.append(sigfig((dP_m+dP_r)/(P_m0+P_r0)*100, 3))
    
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
# jsonfn = os.path.join(results_dir, folder + ".json")
# csvfn = os.path.join(results_dir, folder + ".csv")
# with open(jsonfn, "w") as jsonout:
#     jsonout.write(json.dumps(content, indent=4))
# with open(csvfn, 'w', newline="") as csvout:
#     csvwriter = csv.writer(csvout)
#     csvwriter.writerow(["Variable", "Einheit", "Schrittweite", "dP_mfp/P_mfp [%]", "dP_r/P_r [%]", "dQ/Q [%]"])
#     for i in range(len(par_name_list)):
#         csvwriter.writerow([par_name_list[i], par_unit_list[i], step_size[i], relative_difference[i, 0]*1e2, relative_difference[i, 1]*1e2, relative_difference[i, 2]*1e2])


txtfn = os.path.join(results_dir, folder + ".txt")
with open(txtfn, "w") as f:
    for i in range(len(par_name_list)):
        if par_name_list[i] == "qm_cb0":
            continue
        sz_string = str(sigfig(step_size[i], 3)) if step_size[i] < 0 else "+" + str(sigfig(step_size[i], 3))
        rsz_string = str(sigfig(rel_step_size[i], 3)*100) if rel_step_size[i] < 0 else "+" + str(sigfig(rel_step_size[i], 3)*100)
        dP_string = str(sigfig(deltaP[i], 3)) if deltaP[i] < 0 else "+" + str(sigfig(deltaP[i], 3))
        dQ_string = str(sigfig(deltaQ[i], 3)) if deltaQ[i] < 0 else "+" + str(sigfig(deltaQ[i], 3))
        rsz_string = fixfloat(rsz_string)
        sz_string = fixfloat(sz_string)
        dP_string = fixfloat(dP_string)
        dQ_string = fixfloat(dQ_string)
        
        f.write((par_hr_list[i] + " & " + par_unit_list[i] + " & " + sz_string + " & " + rsz_string + " & " + dP_string + " \\\ \hline \n").replace(".", ",").replace("+inf", "-"))
        
    
    
    


