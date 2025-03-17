# -*- coding: utf-8 -*-
import os
import pathlib
import csv
import json

folder = os.path.join(os.getcwd(), "results2")

subfolders = ["pump", "after", "pre", "dual"]

def get_Data(folder, subfolder, write=True):  
    sf_path = os.path.join(folder, subfolder)
    # find files in target directory
    files = [f for f in os.listdir(sf_path) if os.path.isfile(
        os.path.join(pathlib.Path().resolve(), sf_path, f))]
    P_mfp = list()
    P_r = list()
    t_bk = list()
    t_wu = list()
    qm_r = list()
    qm_v = list()
    qm_cb = list()
    qm_phc = list()
    Q = list()
    params = dict()
    for file in files:
        # only import .csv files
        if file.split(".")[1] == "csv":
            filename = os.path.join(sf_path, file)
            if "FAILED" == filename[-10:-4]:
                continue
            with open(filename, newline="") as f:
                data_f = csv.reader(f)
                param_names = next(data_f)
                param_row = next(data_f)
                _ = next(data_f)
                P_row = next(data_f)
                _ = next(data_f)
                Q_row = next(data_f)
                _ = next(data_f)
                qm_row = next(data_f)
                t_bk.append(float(param_row[1]))
                t_wu.append(float(param_row[2]))
                P_mfp.append(float(P_row[0]))
                Q.append(float(Q_row[0]))
                qm_r.append(float(qm_row[2]))
                qm_cb.append(float(qm_row[0]))
                qm_phc.append(float(qm_row[1]))
                if param_row[0] == "dual":
                    qm_v.append(float(qm_row[3]))
                P_r.append(float(P_row[1]))
                params = dict()
                for i, name, value in zip(range(len(param_names)), param_names, param_row):
                    if i in [1,2,6]:
                        continue
                    params.update({name: value})
    content = {"Parameters": params, "t_bk": t_bk, "t_wu": t_wu, "P_mfp": P_mfp, "P_r": P_r, "Q": Q, "qm_cb": qm_cb, "qm_r": qm_r, "qm_v": qm_v, "qm_phc": qm_phc}
    json_object = json.dumps(content, indent=4)
    json_fn = os.path.join(folder, subfolder+".json")
    if write:
        with open(json_fn, "w") as jsonout:
            jsonout.write(json_object)
    return content


if __name__ == "__main__":
    data = dict()
    for subfolder in subfolders:
        content = get_Data(folder, subfolder)
        name = subfolder[:-5]
        data.update({name: content})

    
    

