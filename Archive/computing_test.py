# -*- coding: utf-8 -*-

import ray
import matplotlib.pyplot as plt
import fuelsystem as fs

t_bk = range(250, 400, 50)

t_wu = 250
eta_p = 0.88
eta_r = 0.9
p_cb = 1.5e6
qm_cb = 0.1
t0 = 22
p0 = 4.2e5

ray.init()

@ray.remote
def pump(t, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0):
    return fs.h2pump(t, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0)


parameters = list()
results = list()
for t in t_bk:
    results.append(pump.remote(t, t_wu, eta_p, eta_r, p_cb, qm_cb, t0, p0))
                   
for result in results:
    print(ray.get(result))