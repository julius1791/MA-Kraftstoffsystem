# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pathlib
import os
from h2_stoffmodell import calculate_H2_cp as H2_cp
from h2_stoffmodell import calculate_H2_enthalpy as H2_enth
from h2_stoffmodell import calculate_H2_rho as H2_rho
import matplotlib.pyplot as plt
import math
import json


def H2_h(t, p, sat):
    if sat == 0:
        h = H2_enth(t, p)
    elif sat == 1:
        h = H2_enth(t, p)
    else:
        t_sat = sat_t(p, True)
        h_v = H2_h(t_sat + 0.1, p)
        h_l = H2_h(t_sat - 0.1, p)
        h = sat*h_v + (1-sat)*h_l
    return h

class Jeta_Flow:
    def __init__(self, mass_flow: float, temperature: float, pressure: float):
        """
        Init method for jeta

        Parameters
        ----------
        mass_flow : float
            mass flow of jet a (kg/s)
        temperature : float
            temperature (K)
        pressure : float
            pressure (Pa)
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure
        
    def heat(self, Q_dot: float):
        """
        Method for heat added to jet fuel in heat exchanger

        Parameters
        ----------
        Q_dot : float
            Absolute thermal power added (W)

        Returns
        -------
        t1 : float
            final temperature (K)

        """
        # specific heat
        q = Q_dot/self.qm
        # initial specific enthalpy
        h0 = jeta_properties(self.t, self.p)[1]
        # final enthalpy
        h1 = h0 + q/1000
        # final temperature
        t1 = jeta_find_t_for_h(h1, self.p)
        self.t = t1
        return t1
    
    def pump(self, p1: float, eta: float):
        """
        Method for raising pressure of jet fuel

        Parameters
        ----------
        p1 : float
            final pressure (Pa)
        eta : float
            pump efficiency (N/A)

        Returns
        -------
        P : float
            pump power (W)
        t1 : float
            final temperature (K)

        """
        # initial density and specific enthalpy
        rho = jeta_density(self.t, self.p)
        h0 = jeta_properties(self.t, self.p)[1]
        # reversible pump work
        a_rev = (p1-self.p)/rho/1000
        # pump work
        a = a_rev/eta
        # final enthalpy and temperature
        h1 = a + h0
        t1 = jeta_find_t_for_h(h1, p1)
        self.p = p1
        self.t = t1
        # pump power
        P = a*self.qm
        return P, t1
    
    def mix(self, secondary_flow):
        """
        Method for mixing two jeta flows into one

        Parameters
        ----------
        secondary_flow : Jeta_Flow
            secondary jeta flow 

        Returns
        -------
        t1 : float
            final temperature (K)

        """
        # find initial enthalpy 
        H0_1 = jeta_properties(self.t, self.p)[1]*self.qm
        H0_2 = jeta_properties(secondary_flow.t, secondary_flow.p)[1]*secondary_flow.qm
        
        # calculate final enthalpy
        H1 = H0_1 + H0_2
        # final mass flow
        qm1 = (self.qm+secondary_flow.qm)
        # final specific enthalpy
        h1 = H1/qm1
        # final pressure
        p1 = min(self.p, secondary_flow.p)
        # final temperature
        t1 = jeta_find_t_for_h(h1, p1)
        self.t = t1
        self.p = p1
        self.qm = qm1
        return t1
    
    def split(self, qm_split: float):
        """
        Method for splitting one flow instance into two separate flows

        Parameters
        ----------
        qm_split : float
            mass flow of secondary flow

        Returns
        -------
        secondary_flow : Jeta_Flow
            secondary flow split off from main jet a flow

        """
        qm0 = self.qm
        self.qm = qm0 - qm_split
        # secondary flow inherits pressure and temperature from primary flow
        secondary_flow = Jeta_Flow(qm_split, self.t, self.p)
        return secondary_flow
    
class H2_Flow:
    def __init__(self, mass_flow: float, temperature: float, pressure: float, saturation: float):
        """
        Init method for hydrogen

        Parameters
        ----------
        mass_flow : float
            mass flow of h2 (kg/s)
        temperature : float
            temperature (K)
        pressure : float
            pressure (Pa)
        saturation : float
            physical state of hydrogen 0 -> super cooled liquid, 1 -> overheated gas, 0~1 -> saturation state
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure
        self.sat = saturation
        
    def heat(self, Q_dot: float):
        """
        Method for heat added to hydrogen in heat exchanger

        Parameters
        ----------
        Q_dot : float
            Absolute thermal power added (W)

        Returns
        -------
        t1 : float
            final temperature (K)

        """
        # specific heat
        q = Q_dot/self.qm
        # initial specific enthalpy
        h0 = H2_h(self.t, self.p)
        # final enthalpy
        h1 = h0 + q
        # final temperature
        t1, saturation = H2_find_t_for_h(h1, self.p)
        self.t = t1
        self.sat = saturation
        return t1
    
    def vaporiser(self):
        """
        Calculate the power needed to vaporise the hydrogen

        Returns
        -------
        Q_dot : float
            Thermal power required for vaporisation

        """
        h0 = H2_h(self.t, self.p)
        t1 = sat_t(self.p, True) + 0.1
        h1 = H2_h(t1, self.p)
        Q_dot = h1 - h0
        self.t = t1
        self.sat = 1
        return Q_dot
    
    def compressor(self, p1: float, eta: float):
        """
        

        Parameters
        ----------
        p1 : float
            pressure (Pa)
        eta : float
            isentropic compressor efficiency

        Returns
        -------
        P : float
            compressor power (W)
        t1 : float
            final temperature (K)

        """
        
        P = 0
        t1 = 0
        return P, t1
    
    def pump(self, p1: float, eta: float):
        """
        Method for raising pressure of hydrogen

        Parameters
        ----------
        p1 : float
            final pressure (Pa)
        eta : float
            pump efficiency (N/A)

        Returns
        -------
        P : float
            pump power (W)
        t1 : float
            final temperature (K)

        """
        # initial density and specific enthalpy
        rho = jeta_density(self.t, self.p)
        h0 = jeta_properties(self.t, self.p)[1]
        # reversible pump work
        a_rev = (p1-self.p)/rho/1000
        # pump work
        a = a_rev/eta
        # final enthalpy and temperature
        h1 = a + h0
        t1 = jeta_find_t_for_h(h1, p1)
        self.p = p1
        self.t = t1
        # pump power
        P = a*self.qm
        return P, t1
    
    def mix(self, secondary_flow):
        """
        Method for mixing two hydrogen flows into one

        Parameters
        ----------
        secondary_flow : H2_Flow
            secondary hydrogen flow 

        Returns
        -------
        t1 : float
            final temperature (K)

        """
        # find initial enthalpy 
        H0_1 = H2_h(self.t, self.p)*self.qm
        H0_2 = H2_h(secondary_flow.t, secondary_flow.p)*secondary_flow.qm
        
        # calculate final enthalpy
        H1 = H0_1 + H0_2
        # final mass flow
        qm1 = (self.qm+secondary_flow.qm)
        # final specific enthalpy
        h1 = H1/qm1
        # final pressure
        p1 = min(self.p, secondary_flow.p)
        # final temperature
        t1 = H2_find_t_for_h(h1, p1)
        self.t = t1
        self.p = p1
        self.qm = qm1
        return t1
    
    def split(self, ratio: float):
        """
        Method for splitting one flow instance into two separate flows

        Parameters
        ----------
        ratio : float
            ratio between mass flow of initial flow and secondary outflow

        Returns
        -------
        secondary_flow : H2_Flow
            secondary flow split off from main hydrogen flow

        """
        qm0 = self.qm
        # mass flow of the split flows
        qm1 = qm0*(1-ratio)
        qm2 = qm0*ratio
        self.qm = qm1
        # secondary flow inherits pressure and temperature from primary flow
        secondary_flow = H2_Flow(qm2, self.t, self.p)
        return secondary_flow

def jeta_antoine():
    """
    Generate antoine equation constants for jet-a

    Returns
    -------
    A, B, C : float
        Antoine equation parameters
    """
    A = 8.81923182000836
    B = 1374.12563
    C = 502.76012
    return A, B, C

def H2_antoine():
    """
    Generate antoine equation constants for hydrogen

    Returns
    -------
    A, B, C : float
        Antoine equation parameters
    """
    A = 8.54314
    B = 99.395
    C = 7.726
    return A, B, C

def sat_p(t: float, substance: bool):
    """
    Find saturation pressure for a given Temperature

    Parameters
    ----------
    t : float
        Temperature (K)
    substance : bool 
        True -> H2
        False -> jet a

    Returns
    -------
    p : float
        Saturation pressure (Pa)

    """
    # define antoine equation constants
    if substance:
        A, B, C = H2_antoine()
    else:
        A, B, C = jeta_antoine()
    # apply antoine equation
    p = 10**(A-(B/(t + C)))
    return p


def sat_t(p: float, substance: bool):
    """
    Calculate saturation temperature for given pressure

    Parameters
    ----------
    p : float
        pressure (Pa)
    substance : bool 
        True -> H2
        False -> jet a

    Returns
    -------
    t : float
        temperature (K)

    """
    # define antoine equation constants
    if substance:
        A, B, C = H2_antoine()
    else:
        A, B, C = jeta_antoine()
    # apply antoine equation
    t = B/(A-math.log(p, 10))-C
    return t


def interpolate(X, Y, x: float):
    """
    Linear interpolation along axis X
    
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
    # find the nearest point above point of interest
    x1 = X[X > x].min()
    # find the nearest point below point of interest
    x0 = X[X < x].max()
    # get indeces of the known points
    id_1 = np.where(X == x1)
    id_0 = np.where(X == x0)
    # get ordinates of the known points
    y1 = Y[id_1]
    y0 = Y[id_0]
    # apply linear interpolation
    y = y0 + (y1-y0)/(x1-x0)*(x-x0)
    return y


def import_tpp():
    """
    Import all thermophysical property data stored on .csv files in "stoffdaten" directory 

    Returns
    -------
    tpp : dict
        dictionary indexed by filename with thermo-physical properties

    """
    # find files in target directory
    files = [f for f in os.listdir(os.path.join(pathlib.Path().resolve(), "stoffdaten")) if os.path.isfile(
        os.path.join(pathlib.Path().resolve(), "stoffdaten", f))]
    # init dictionary for tpp
    tpp = dict()
    # iterate list of files
    for file in files:
        # only import .csv files
        if file.split(".")[1] == "csv":
            # use filename as index (without extension)
            name = file.split(".")[0]
            # feed file contents into pandas dataframe
            content = pd.read_csv(os.path.join(pathlib.Path().resolve(), "stoffdaten", file), delimiter=",")
            tpp.update({name: content})
    return tpp



def H2_find_t_for_h(h_t: float, p: float):
    """
    Find the temperature corresponding to an enthalpy and pressure of hydrogen

    Parameters
    ----------
    h_t : float
        target enthalpy (kJ/kg)
    p : float
        pressure (Pa)

    Returns
    -------
    t : float
        Temperature
    saturation : float
        physical state of hydrogen 0 -> super cooled liquid, 1 -> overheated gas, 0~1 -> saturation state

    """
    # find saturation temperature
    t_sat = sat_t(p, True)
    p_crit = 1.2964e6
    if p > p_crit:
        t = 300
        tmax = 2000
        tmin = 20
    else:
        # determine bounds and initial temperature depending on the physical state
        # determine if the hydrogen is an overheated gas
        if h_t > H2_h(t_sat + 0.1, p):
            tmin = t_sat + 0.1
            t = tmin + 10
            tmax = 2000
            saturation = 1
    
        # determine if the hydrogen is in the saturation regime
        elif h_t > H2_h(t_sat - 0.1, p):
            h_v = H2_h(t_sat + 0.1, p)
            h_l = H2_h(t_sat - 0.1, p)
            saturation = (h_t - h_l)/(h_v-h_l)
            t = t_sat
            
            return t, saturation
    
        # else it is a supercooled liquid
        else:
            saturation = 0
            tmax = t_sat - 0.1
            t = tmax - 5
            tmin = 1
    
    # loop until target enthalpy is achieved
    while abs(H2_h(t, p)-h_t) > 1e-9:
        # calculate temperature step asssuming constant heat capacity
        dt = (h_t-H2_h(t, p))/H2_cp(t, p)

        t = t + dt

        # ensure that the new calculated temperature is within the bounds
        # for the current state
        if t > tmax:
            t = tmax
        if t < tmin:
            t = tmin
    return t, saturation


def jeta_find_t_for_h(h_t: float, p: float):
    """
    Find the temperature corresponding to an enthalpy and pressure of jet fuel

    Parameters
    ----------
    h_t : float
        target enthalpy (kJ/kg)
    p : float
        pressure (Pa)

    Returns
    -------
    t : float
        Temperature

    """
    t = 300
    cp, h = jeta_properties(t, p)
    # loop until target enthalpy is achieved
    while abs(h-h_t) > 1e-9:
        cp, h = jeta_properties(t, p)
        # calculate temperature step asssuming constant heat capacity
        dt = (h_t-h)/cp

        t = t + dt
    return t


def jeta_properties(t: float, p: float):
    """
    Calculate Heat capacity and enthalpy of jet fuel for a given temperature

    Parameters
    ----------
    t : float
        Temperature (K)
    p : float
        Pressure (Pa)

    Returns
    -------
    cp : float
        specific heat capacity of jet fuel (kJ/kgK)
    h : float
        enthalpy of jet fuel (kJ/kg)

    """
    # define universal gas constant
    R = 8.3145      # Jmol-1K
    # import data from json
    with open(os.path.join(pathlib.Path().resolve(), "stoffdaten", "jeta.json")) as f:
        jeta = json.load(f)
    M_jeta = float(jeta["molecular_weight"])        # gmol-1
    R_jeta = R/M_jeta                               # kJ/kgK
    
    # init variables for loop
    a_list = list()
    i = 0
    cp = 0
    exponent = -2
    
    # iterate over polynomial coefficients from json data
    for a in jeta["a_coefficients"].split(" "):
        # remove bracket curls, convert to float and append to list
        a = a.replace("]", "")
        a_list.append(float(a.replace("[", "")))
        
        # sum up polynomial for specific heat capacity
        cp = cp + a_list[i]*pow(t, exponent)
        
        # prepare next iteration
        exponent += 1
        i += 1
        
    # multipy with specific gas constant
    cp = cp * R_jeta
    
    
    rho = jeta_density(t, p)
    
    # sum up enthalpy components
    h = -a_list[0]/t
    h += a_list[1]*math.log(t)
    h += a_list[2]*t
    h += a_list[3]*pow(t, 2)/2
    h += a_list[4]*pow(t, 3)/3
    h += a_list[5]*pow(t, 4)/4
    h += a_list[6]*pow(t, 5)/5
    h = h*R_jeta
    h += p/rho/1000

    return cp, h

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


def jeta_density(t: float, p: float):
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
    jeta = import_tpp()["jeta"]
    # feed the relevant data into numpy arrays
    temp = jeta[["Temperature (K)"]].to_numpy()
    press = jeta[["Pressure (Mpa)"]].to_numpy()
    rho_arr = jeta[["Density (kg/m3)"]].to_numpy()
    # pressure unit conversion (Pa to MPa)
    p = p/10e6
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



# jeta_props = jeta_properties(300, 1e5)
# dh = jeta_properties(400, 1e5)[1]-jeta_properties(300, 1e5)[1]
# print(jeta_props)
# print(dh)
# jeta_props = jeta_properties(400, 1e5)
# print(jeta_props)
# tpp = import_tpp()

# h2_2_cp = tpp['h2_2bar']["Cp (J/g*K)"].to_numpy()
# h2_2_t = tpp['h2_2bar']["Temperature (K)"].to_numpy()
# h2_30_cp = tpp['h2_30bar']["Cp (J/g*K)"].to_numpy()
# h2_30_t = tpp['h2_30bar']["Temperature (K)"].to_numpy()

# diff_2 = list()
# diff_30 = list()
# cp_2 = list()
# cp_30 = list()
# for i in range(np.prod(h2_2_cp.shape)):
#     cp_2.append(H2_cp(h2_2_t[i], 2*10**5)/1000)
#     diff_2.append(abs(1 - 1/(h2_2_cp[i]/cp_2[i])))
# for i in range(np.prod(h2_30_cp.shape)):
#     cp_30.append(H2_cp(h2_30_t[i], 30*10**5)/1000)
#     diff_30.append(abs(1 - 1/(h2_30_cp[i]/cp_30[i])))
# fig, axs = plt.subplots(2)


# axs[0].plot(h2_2_t, h2_2_cp)
# axs[0].plot(h2_2_t, cp_2)
# axs[1].plot(h2_30_t, h2_30_cp)
# axs[1].plot(h2_30_t, cp_30)

# axs[0].set(ylim=(10, 14))
# axs[1].set(ylim=(10, 36))
# axs[0].title.set_text("2 bar")
# axs[1].title.set_text("30 bar")
# for i in range(2):
#     axs[i].set(xlim=(0, 200), xlabel="Temperature [K]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
#     axs[i].legend(["nist.gov", "stoffmodell"])
# fig.tight_layout()
# plt.savefig("h2stoffmodell_abweichungen.png", dpi=600)

# min_p = 10**5
# max_p = 2*10**6
# p_increment = 10**4
# cp_t30 = list()
# cp_t40 = list()
# cp_t300 = list()
# cp_t600 = list()
# p_list = list()

# for i in range(min_p, max_p, p_increment):
#     p_list.append(i)
#     cp_t30.append(H2_cp(30, i)/1000)
#     cp_t40.append(H2_cp(40, i)/1000)
#     cp_t300.append(H2_cp(300, i)/1000)
#     cp_t600.append(H2_cp(600, i)/1000)

# fig, ax = plt.subplots()
# ax.plot(p_list, cp_t30)
# ax.plot(p_list, cp_t40)
# ax.plot(p_list, cp_t300)
# ax.plot(p_list, cp_t600)

# ax.legend(["T = 30K", "T = 40K", "T = 300K", "T = 600K"])
# ax.set(xlabel="Pressure [Pa]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
# ax.set_xlim(xmin=min_p, xmax=max_p)
# plt.savefig("isothermen_cp", dpi=600)

# min_t = 10
# max_t = 200
# t_increment = 1
# cp_p1 = list()
# cp_p10 = list()
# cp_p15 = list()
# cp_p20 = list()
# cp_p30 = list()
# h_p1 = list()
# h_p10 = list()
# h_p15 = list()
# h_p20 = list()
# h_p30 = list()
# t_list = list()

# for i in range(min_t, max_t, t_increment):
#     i = i/2
#     t_list.append(i)
#     cp_p1.append(H2_cp(i, 10**5)/1000)
#     cp_p10.append(H2_cp(i, 10**6)/1000)
#     cp_p15.append(H2_cp(i, 1.5*10**6)/1000)
#     cp_p20.append(H2_cp(i, 2*10**6)/1000)
#     cp_p30.append(H2_cp(i, 3*10**6)/1000)
#     h_p1.append(H2_h(i, 10**5)/1000)
#     h_p10.append(H2_h(i, 10**6)/1000)
#     h_p15.append(H2_h(i, 1.5*10**6)/1000)
#     h_p20.append(H2_h(i, 2*10**6)/1000)
#     h_p30.append(H2_h(i, 3*10**6)/1000)

# fig, ax = plt.subplots()
# ax.plot(t_list, cp_p1)
# ax.plot(t_list, cp_p10)
# ax.plot(t_list, cp_p15)
# ax.plot(t_list, cp_p20)
# ax.plot(t_list, cp_p30)

# ax.legend(["p = 1 bar", "p = 10 bar", "p = 15 bar",  "p = 20 bar", "p = 30 bar"])
# ax.set(xlabel="Temperature [K]", ylabel=r"$\ c_p [\frac{kJ}{kgK}]$ ")
# ax.set_xlim(xmin=min_t, xmax=max_t/2)
# plt.savefig("isobaren_cp", dpi=600)

# fig, ax = plt.subplots()
# ax.plot(t_list, h_p1)
# ax.plot(t_list, h_p10)
# ax.plot(t_list, h_p15)
# ax.plot(t_list, h_p20)
# ax.plot(t_list, h_p30)

# ax.legend(["p = 1 bar", "p = 10 bar", "p = 15 bar",  "p = 20 bar", "p = 30 bar"])
# ax.set(xlabel="Temperature [K]", ylabel=r"$\ h [\frac{kJ}{kg}]$ ")
# ax.set_xlim(xmin=min_t, xmax=max_t/2)
# plt.savefig("isobaren_h", dpi=600)

# t_hot = H2_find_t_for_h(H2_h(80, 1e5), 1e5)
# print(t_hot)



# ff1 = Jeta_Flow(1, 300, 1e5)
# ff1.heat(50000)
# print(ff1.t)
# ff1.heat(-50000)
# print(ff1.t)

# print(jeta_properties(300, 1e5)[0])
# print(jeta_properties(300, 3*1e6)[1]-jeta_properties(300, 1e6)[1])
# work = ff1.pump(3e6, 0.8)
# print(ff1.p, ff1.t)
# print(work)




t_r1 = 420
qm_cb = 0.3
qm_r = 0.3
qm_t = 0.1
t0 = 250
p0 = 0.4e5
p_lpfp = 3e5
eta_lpfp = 0.83
Q_fohe = 200000
p_hpfp = 3e6
eta_hpfp = 0.88
Q_idg = 5500

t_r0 = 1000
i = 0
while abs(t_r0 - t_r1) > 1e-6:
    jetaflow = Jeta_Flow(qm_cb+qm_t, t0, p0)
    i+=1
    t_r0 = t_r1
    t_lpfp = jetaflow.pump(p_lpfp, eta_lpfp)
    t_mix = jetaflow.mix(Jeta_Flow(qm_r, t_r0, jetaflow.p))
    t_fohe = jetaflow.heat(Q_fohe)
    t_hpfp = jetaflow.pump(p_hpfp, eta_hpfp)
    cb_ff = jetaflow.split(qm_cb)
    t_idg = jetaflow.heat(Q_idg)
    to_tank_ff = jetaflow.split(qm_t)
    t_r1 = jetaflow.t