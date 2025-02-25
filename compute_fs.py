# -*- coding: utf-8 -*-
import fuelsystem2
import os
from joblib import Parallel, delayed
import time

t_bk_list = range(200, 620, 20)
t_wu_list = range(100, 300, 20)
eta_p = 0.92
eta_r = 0.9
p_cb = 1.5e6+168.9e3
qm_cb = 0.1
t0 = 22
p0 = 4.2e5
pcc = True
tpr_hx = 0.95

folder = "after"
param_combinations = []

for t_bk in t_bk_list:
    for t_wu in t_wu_list:
        if t_wu + 20 > t_bk:
            continue
        filename = str(t_bk) + str(t_wu) + ".csv"
        path = os.path.join(folder, filename)
        param_combinations.append([t_bk, t_wu, path])

start = time.time()
Parallel(n_jobs=20)(
    delayed(fuelsystem2.h2after)(t_bk, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0, tpr_hx, Q_hx=200e3, pcc=pcc, filename=path) for t_bk, t_wu, path in param_combinations
)
stop = time.time()
print(stop-start)