# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pathlib
import os


        
def interpolate(X, Y, x: float):
    """
    Parameters
    ----------
    X : array_like 
        sorted list of Physical properties forming the abscissa of the interpolation
    Y : array_like 
        sorted list of Physical properties forming the ordinate of the interpolation
    x : float 
        Point of interest on the abscissa

    Returns
    -------
    y : float
        Point of interest on the ordinate

    """
    "find the nearest point above point of interest"
    x1 = X[X > x].min()     
    "find the nearest point below point of interest"
    x0 = X[X < x].max()
    "get indeces of the known points"
    id_1 = np.where(X == x1)
    id_0 = np.where(X == x0)
    "get ordinates of the known points"
    y1 = Y[id_1]
    y0 = Y[id_0]
    "apply linear interpolation"
    y = y0 + (y1-y0)/(x1-x0)*(x-x0)
    return y

def import_tpp():
    files = [f for f in os.listdir(os.path.join(pathlib.Path().resolve(), "stoffdaten")) if os.path.isfile(os.path.join(pathlib.Path().resolve(), "stoffdaten", f))]
    tpp = dict()
    for file in files:
        name = file.split(".")[0]
        content = import_chemdata(file)
        print(content)
        tpp.update({name: content})
    return tpp

def import_chemdata(filename: str):
    """
    Parameters
    ----------
    filename : str
        filename of physical property data file (must be placed in MA Kraftstoffsystem\stoffdaten)

    Returns
    -------
    df : pandas DataFrame
        DataFrame containing the physical property data read from the text file

    """
    path = os.path.join(pathlib.Path().resolve(), "stoffdaten", filename)
    df = pd.read_csv(path, delimiter=",")
    return df

tpp = import_tpp()
print(list(tpp["h2sat"]))
print(tpp["h2sat"][["Temperature (K)", "Pressure (MPa)"]])

