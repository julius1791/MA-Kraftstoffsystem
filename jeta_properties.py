import os
import pathlib
import json
import pandas as pd
import math

# global variables
prop_dir = "material properties"
R = 8.3145                                                          # Jmol-1K

                         

# find files in target directory
files = [f for f in os.listdir(os.path.join(pathlib.Path().resolve(), prop_dir)) if os.path.isfile(
    os.path.join(pathlib.Path().resolve(), prop_dir, f))]
# init dictionary for tpp
tpp = dict()
# iterate list of files
for file in files:
    # only import .csv files
    if file.split(".")[1] == "csv":
        # use filename as index (without extension)
        name = file.split(".")[0]
        # feed file contents into pandas dataframe
        content = pd.read_csv(os.path.join(pathlib.Path().resolve(), prop_dir, file), delimiter=",")
        tpp.update({name: content})


                                   
# import data from json
with open(os.path.join(
        pathlib.Path().resolve(), prop_dir, "jeta.json")) as f:
    jeta_props = json.load(f)
M_jeta = float(jeta_props["molecular_weight"])                  # gmol-1
R_jeta = R/M_jeta                                               # kJ/kgK  

# iterate over split string, remove brackets and convert to float
a_list = list()
for a in jeta_props["a_coefficients"].split(" "):
    a = a.replace("]", "")
    a = a.replace("[", "")
    a_list.append(float(a))

###################################################
# ## c_p  (isobare Wärmekapazität)  [J/(kg*K)] ## #
###################################################

def calc_jeta_cp(t: float, p: float):
    """function for calculating specific heat capacity of jet a"""
    # import specific gas constant and polynomial coefficients
    
    # init variables for calculation loop
    i = 0
    cp = 0
    exponent = -2
    
    # iterate over polynomial coefficients from json data
    for a in a_list:     
        # sum up polynomial for specific heat capacity
        cp += a*pow(t, exponent)
        # prepare next iteration
        exponent += 1
        i += 1
        
    # multipy with specific gas constant + unit conversion
    cp = cp * R_jeta*1000                                           # J/kgK
    return cp

#############################################
# ########## Enthalpie h [J/kg] ########### #
#############################################

def calc_jeta_enthalpy(t: float, p: float):

    rho = calc_jeta_density(t, p)
    
    # sum up enthalpy components
    h = -a_list[0]/t
    h += a_list[1]*math.log(t)
    h += a_list[2]*t
    h += a_list[3]*pow(t, 2)/2
    h += a_list[4]*pow(t, 3)/3
    h += a_list[5]*pow(t, 4)/4
    h += a_list[6]*pow(t, 5)/5
    # unit conversion
    h =  h*R_jeta*1000                                              # J/kg                    
    h += p/rho
    return h

#############################################
# ########## Entropie s [J/kgK] ########### #
#############################################
    
def calc_jeta_entropy(t: float, p: float):
    
    # sum up enthalpy components
    s =  -a_list[0]*pow(t, -2)/2
    s += -a_list[1]/t
    s += a_list[2]*math.log(t)
    s += a_list[3]*t
    s += a_list[4]*pow(t, 2)/2
    s += a_list[5]*pow(t, 3)/3
    s += a_list[6]*pow(t, 4)/4
    # unit conversion
    s =  s*R_jeta*1000                                              # J/kgK
    return s

#############################################
# ########## Dichte rho [kg/m3] ########### #
#############################################

def calc_jeta_density(t: float, p: float):
    """
    Calculate the density of jet fuel given temperature and pressure
    
    Parameters
    ----------
    t : float 
        Temperature (K)
    p : float 
        Pressure (Pa)

    Returns
    -------
    rho : float
        density of jet fuel (kg/m3)

    """
    # import thermo-physical properties of jet fuel from .csv 
    jeta_tpp = tpp["jeta"]
    # feed the relevant data into numpy arrays
    temp = jeta_tpp[["Temperature (K)"]].to_numpy()
    press = jeta_tpp[["Pressure (Mpa)"]].to_numpy()
    rho_arr = jeta_tpp[["Density (kg/m3)"]].to_numpy()
    # pressure unit conversion (Pa to MPa)
    p = p/10e6
    
    # use the 4 closest points to point of interest to 
    # bilinearly inter-/extrapolate
    
    # find closest Temperature values
    t0, t1 = find_t(t, temp)
    # select the pressure and density values corresponding to the selected temperature values
    press_t1 = press[temp == t1]
    press_t0 = press[temp == t0]
    rho_t1 = rho_arr[temp == t1]
    rho_t0 = rho_arr[temp == t0]
    # find closest pressure values
    t0p0, t0p1 = find_p(p, press_t0)
    t1p0, t1p1 = find_p(p, press_t1)
    # find the density at the previously defined temperature and pressure
    rho_t0p1 = rho_t0[press_t0 == t0p1]
    rho_t0p0 = rho_t0[press_t0 == t0p0]
    rho_t1p1 = rho_t1[press_t1 == t1p1]
    rho_t1p0 = rho_t1[press_t1 == t1p0]
    # interpolate along the pressure axis for the two temperature values
    rhot1 = rho_t1p1 + (rho_t1p1-rho_t1p0)/(t1p1-t1p0)*(p-t1p1)
    rhot0 = rho_t0p1 + (rho_t0p1-rho_t0p0)/(t0p1-t0p0)*(p-t0p1)
    # interpolate along the temperature axis
    rho = rhot1 + (rhot1-rhot0)/(t1-t0)*(t-t1)
    return rho[0]

def find_t(t: float, temp):
    """
    Return neighbouring temperature values from list or two largest/smallest values when extrapolating

    Parameters
    ----------
    t : float
        Temperature (K)
    temp : numpy array
        Temperatures for which data is available (K)

    Returns
    -------
    t0 : float
        smaller temperature value  (K)
    t1 : float
        larger temperature value (K)
    """
    # when extrapolating an exception is triggered. Then, use two smallest or largest temperature values for extrapolation
    # otherwise when interpolating use the closest temperature values to the point of interest 
    try:
        # when interpolating set the upper temperature to the smallest value larger than the point of interest
        t1 = temp[temp > t].min()
    # when no temperature larger than the point of interest exists, extrapolation occurs and an exception is triggered
    except:
        # set the upper temperature to the largest temperature value
        t1 = temp.max()
        # set the lower temperature to the next largest temperature value
        t0 = temp[temp < t1].max()
        return t0, t1
    try:
        t0 = temp[temp < t].max()
    except:
        t0 = temp.min()
        t1 = temp[temp > t0].min()
        return t0, t1
    return t0, t1


def find_p(p: float, press):
    """
    Return neighbouring pressure values from list or two largest/smallest values when extrapolating

    Parameters
    ----------
    p : float
        pressure (MPa)
    press : numpy array
        pressures for which data is available (MPa)

    Returns
    -------
    p0 : float
        smaller pressure value (MPa)
    p1 : float
        larger pressure value (MPa)
    """
    # when extrapolating an exception is triggered. Then, use two smallest or largest pressure values for extrapolation
    # otherwise when interpolating use the closest pressure values to the point of interest 
    try:
        # when interpolating set the upper pressure to the smallest value larger than the point of interest
        p1 = press[press > p].min()
    # when no pressure larger than the point of interest exists, extrapolation occurs and an exception is triggered
    except:
        # set the upper pressure to the largest temperature value
        p1 = press.max()
        # set the lower pressure to the next largest temperature value
        p0 = press[press < p1].max()
        return p0, p1
    try:
        p0 = press[press < p].max()
    except:
        p0 = press.min()
        p1 = press[press > p0].min()
        return p0, p1
    return p0, p1


