# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
import time
from gather_data import get_Data

t_bk_list = range(200, 220, 20)
t_wu_list = range(100, 300, 20)
eta_p = 0.92
eta_r = 0.9
p_cb = 1.5e6+168.9e3
qm_cb = 0.10998
t0 = 22
p0 = 4.2e5
pcc = True
tpr_hx = 0.95


# results folder
results_dir = "results2"

# subfolder
folders = ["pump2", "after2", "pre2", "dual2"]
functions = [fuelsystem2.h2pump, fuelsystem2.h2after, fuelsystem2.h2pre, fuelsystem2.h2dual]
param_combinations = []

for folder, func in zip(folders, functions):
    for t_bk in t_bk_list:
        for t_wu in t_wu_list:
            if t_wu + 20 > t_bk:
                continue
            filename = str(t_bk) + str(t_wu) + ".csv"
            path = os.path.join(folder, filename)
            param_combinations.append([t_bk, t_wu, path])

    start = time.time()
    Parallel(n_jobs=20)(
        delayed(func)(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, tpr_hx, Q_hx=200e3, pcc=pcc, filename=path) for t_bk, t_wu, path in param_combinations
    )
    stop = time.time()
    get_Data(results_dir, folder)
    print(stop-start)