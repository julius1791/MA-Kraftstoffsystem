# -*- coding: utf-8 -*-
import os
import pathlib
import csv
import json

folder = os.path.join(os.getcwd(), "results")

subfolders = ["pump_test", "after_test", "pre_test", "dual_test"]

def get_Data(folder, subfolder):  
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
    for file in files:
        # only import .csv files
        if file.split(".")[1] == "csv":
            filename = os.path.join(sf_path, file)
            with open(filename, newline="") as f:
                data_f = csv.reader(f)
                _ = next(data_f)
                param_row = next(data_f)
                _ = next(data_f)
                P_row = next(data_f)
                _ = next(data_f)
                _ = next(data_f)
                _ = next(data_f)
                qm_row = next(data_f)
                t_bk.append(float(param_row[0]))
                t_wu.append(float(param_row[1]))
                P_mfp.append(float(P_row[0]))
                qm_r.append(float(qm_row[1]))
                if subfolder[:3] == "dual":
                    qm_v.append(float(qm_row[2]))
                if subfolder[:3] != "pre":
                    P_r.append(float(P_row[1]))
    content = {"t_bk": t_bk, "t_wu": t_wu, "P_mfp": P_mfp, "P_r": P_r, "qm_r": qm_r, "qm_v": qm_v}
    json_object = json.dumps(content, indent=4)
    json_fn = os.path.join(folder, subfolder+".json")
    with open(json_fn, "w") as jsonout:
        jsonout.write(json_object)
    return content


if __name__ == "__main__":
    data = dict()
    for subfolder in subfolders:
        content = get_Data(folder, subfolder)
        name = subfolder[:-5]
        data.update({name: content})

    
    

