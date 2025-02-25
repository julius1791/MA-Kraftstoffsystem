# -*- coding: utf-8 -*-

import h2_properties as h2
import math


# modelling parameters
tolerance = 1e-6    # -
max_iter = 100       # -

# hydrogen critical pressure and temerature
p_crit = 1.2964e6   # Pa
t_crit = 33         # K


def antoine():
    """antoine equation parameters for hydrogen"""
    A = 8.54314
    B = 99.395
    C = 7.726
    return A, B, C


###############################################################################
######################## Class for hydrogen fuel flow #########################
###############################################################################

class H2Flow:
    def __init__(
            self, mass_flow: float, temperature: float,
            pressure: float, velocity: float, vapour:bool
        ):
        """
        FuelFlow constructor

        Parameters
        ----------
        mass_flow : float
            mass flow of h2 (kg/s)
        temperature : float
            static temperature (K)
        pressure : float
            static pressure (Pa)
        velocity : float
            velocity of fluid flow (m/s)
        vapour : bool
            physical state of hydrogen True -> vapour, False -> liquid
            
        enthalpy : float
            specific total enthalpy
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure   
        self.v = velocity
        self.vap = vapour
        self.ht = calc_ht(temperature, pressure, velocity)
        
    def heat_exchanger(self, Q_dot:float, tpr: float):
        """
        Method for heat added in heat exchanger

        Parameters
        ----------
        Q_dot : float
            absolute thermal power added (W)
        tpr : float
            total pressure ratio across the heat exchanger (-)

        Returns
        -------
        t1 : float
            final total temperature (K)
        """
        if not self.vap:
            raise ValueError("Expected vapour, not liquid h2")
        # initial specific total enthalpy
        ht0 = self.ht                                                   # J/kg  
        pt0 = calc_pt(self.t, self.p, self.v)               # Pa
        # specific heat
        q = Q_dot/self.qm                                               # J/kg
        
        # final specific total enthalpy
        ht1 = ht0 + q                                                   # J/kg
        
        self.ht = ht1
        
        t1 = calc_t(ht1, self.p, self.v, True)                          # K
        self.t = t1
        
        pt1 = pt0 * tpr
        t1, p1 = apply_total_pressure(pt1, ht1, t1, self.v)
        self.t = t1
        self.p = p1
        
        return t1
    
    def heat_to_saturation(self, tpr):
        """
        Method for adding heat to fuel in heat exchanger, creating a 
        super heated vapour

        Returns
        -------
        Q_dot : float
            thermal power required for vaporisation (W)
        tpr : float
            total pressure ratio across the heat exchanger (-)

        """
        if self.vap:
            raise ValueError("Expected liquid h2")
        
        pt0 = calc_pt(self.t, self.p, self.v)                           # Pa
        ht0 = calc_ht(self.t, self.p, self.v)                           # J/kg
        t1 = sat_t(self.p) + 5
        self.t = t1
        ht1 = calc_ht(t1, self.p, self.v)                               # J/kg
        pt1 = pt0 * tpr
        Q_dot = (ht1 - ht0)*self.qm                                     # W
        self.t, self.p = apply_total_pressure(pt1, ht1, t1, self.v)    
        self.vap = True
        self.ht = ht1
        
        return Q_dot
    
    def pump_hydraulic(self, p1: float, eta: float):
        """
        Method for raising pressure of fuel flow in a hydraulic pump

        Parameters
        ----------
        p1 : float
            final static pressure (Pa)
        eta : float
            pump efficiency (-)

        Returns
        -------
        P : float
            pump power (W)
        t1 : float
            final static temperature (K)
        """     
        # function is meant for pumping into the supercritical state
        if p1 < p_crit:
            raise ValueError(
                "Expected supercritical final pressure. Got p1 = " + str(p1)
            )
        
        # calculate initial specific entropy and specific total enthalpy
        s = h2.calc_H2_entropy(self.t, self.p)
        ht_0 = self.ht
        # guess isentropic temperature
        ts_guess = self.t * (p1/self.p)**0.286
        # find static temperature that yields 
        # the same entropy at higher pressure
        ts1 = find_t_for_s(s, p1, ts_guess)
        # find the specific total enthalpy of the isentropic point 
        ht_s1 = calc_ht(ts1, p1, self.v)
        # find reversible pump work
        a_rev = ht_s1 - ht_0 
        # calculate actual compression work and final total enthalpy
        a = a_rev/eta
        ht_1 = ht_0 + a
        
        # set attributes
        self.p = p1
        self.ht = ht_1
        self.vap = True
        
        # find final static temperature
        t1 = calc_t(ht_1, p1, self.v, True)
        self.t = t1
        
        # calculate compressor power
        P = a*self.qm
        
        return P, t1
    
    def mix_flows(self, secondary_flow):
        """
        Method for mixing two hydrogen fuel flow instances into one
        
        Requirements:
            - Resulting flow must have sufficent enthalpy to be superheated

        Parameters
        ----------
        secondary_flow : H2Flow
            secondary fuel flow 

        Returns
        -------
        float
            final static temperature (K)
        """
        if type(self) != type(secondary_flow):
            raise ValueError("Can't mix flows of different fuels")
        
        # find initial total enthalpy 
        Ht_1 = self.ht*self.qm                                           # J
        Ht_2 = secondary_flow.ht*secondary_flow.qm                       # J
        pt_1 = calc_pt(self.t, self.p, self.v)
        pt_2 = calc_pt(secondary_flow.t, secondary_flow.p, secondary_flow.v)
        
        # final total pressure
        pt1 = min(pt_1, pt_2)
        
        # final total enthalpy
        Ht = Ht_1 + Ht_2                                                # J
        # final mass flow
        qm = (self.qm+secondary_flow.qm)                                # kg/s
        # final specific total enthalpy
        ht = Ht/qm                                                      # J/kg
        
        self.vap = True
        self.qm = qm
        self.ht = ht
        
        t1, p1 = apply_total_pressure(pt1, ht, self.t, self.v)
        self.p = p1
        self.t = t1
        return t1, p1
    
    def split_flows(self, qm: float):
        """
        Method for splitting one fuel flow instance into two separate flows

        Parameters
        ----------
        qm : float
            mass flow of split flow

        Returns
        -------
        FuelFlow
            secondary flow split off from main flow
        """
        # splitting flows is possible regardless physical state
        if qm > self.qm:
            message = "Can't split " + str(qm) + " kg/s because "
            message += "original flow is only " + str(self.qm) + " kg/s"
            raise ValueError(message)
        # secondary flow inherits pressure and enthalpy from primary flow
        self.qm = self.qm - qm
        
        return H2Flow(qm, self.t, self.p, self.v, self.vap)
    
    def reduce_pressure(self, p1):      
        """
        Adiabatic pressure reduction (for instance in a throttle)

        Parameters
        ----------
        p1 : float
            final static pressure (Pa)

        Returns
        -------
        t1 : float
            final static temperature (K)
        """
        if p1 > self.p:
            message = "The selected pressure " + str(p1)
            message += " is greater than initial pressure."
            raise ValueError(message)                                 
        self.p = p1
        t1 = calc_t(self.ht, self.p, self.v, self.vap)
        self.t = t1
        return t1
    
    def copy(self):
        return H2Flow(self.qm, self.t, self.p, self.v, self.vap)
        






def find_t_for_s(s, p, t):
    """
    Find the total temperature corresponding to a given specific entropy and
    static pressure

    Parameters
    ----------
    s : float
        specific entropy (J/kgK)
    p : float
        final static pressure
    t : float
        starting temperature (K)

    Returns
    -------
    t1 : float
        final static temperature (K)
    """
    t1 = t
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1
        t0 = t1
        s_a = h2.calc_H2_entropy(t0, p)
        t1 = t0*math.exp((s-s_a)/h2.calc_H2_cp(t0, p))
        t1 = 0.6*t1 + 0.4*t0
        s_a = h2.calc_H2_entropy(t1, p)
        condition_bool = not abs(s_a-s)/(s_a+s) < tolerance
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1

def apply_total_pressure(pt, ht, t1, v):
    """
    Calculate static temperature and pressure given total pressure and 
    specific total enthalpy, while keeping flow velocity constant

    Parameters
    ----------
    pt : float
        total pressure (Pa)
    ht : float
        specific total enthalpy (J/kg)
    t1 : float
        starting temperature (K)

    Returns
    -------
    t1 : float
        static temperature (K)
    p1 : float
        static pressure (Pa)

    """
    # set starting values for static temperature and pressure
    p1 = pt
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1
        # roll over
        t0, p0 = t1, p1
        # calculate specific total enthalpy 
        ht_a = calc_ht(t0, p0, v)
        # calculate total pressure
        pt_a = p0 + v**2 * h2.calc_H2_density(t0, p0) / 2
        
        # calculate static pressure and temperature for next iteration
        p1 += pt - pt_a
        t1 += (ht-ht_a)/min(h2.calc_H2_cp(t0, p0), 8e3) / 8
        t1 = max(16, t1)
        condition_bool = not (
            abs(ht-ht_a)/(ht+ht_a) < tolerance
            and abs(pt-pt_a)/(pt+pt_a) < tolerance
            )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1, p1

def calc_pt(t, p, v):
    pt = p + h2.calc_H2_density(t, p) * v**2/2
    return pt

def calc_t(ht, p, v, vap):
    
    t1, tmax, tmin = define_calc_t_bounds(ht, p, v, vap)
    
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1 
        t0 = t1
        
        ht_a = calc_ht(t0, p, v)
        t1 += (ht - ht_a)/(h2.calc_H2_cp(t0, p))
        
        t1 = loop_check_bounds(t1, tmax, tmin)
        
        condition_bool = not (
            abs(ht-ht_a)/(ht+ht_a) < tolerance
        )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1
        
def calc_tp(ht, pt, v, vap):
    
    t1, tmax, tmin = define_calc_t_bounds(ht, pt, v, vap)
    p1 = pt
    
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1 
        t0 = t1
        p0 = p1
        
        ht_a = calc_ht(t0, p0, v)
        pt_a = p1 + h2.calc_H2_density(t0, p0) * v**2/2
        t1 += (ht - ht_a)/(h2.calc_H2_cp(t0, p0))
        p1 += pt-pt_a
        
        t1 = loop_check_bounds(t1, tmax, tmin)
        
        condition_bool = not (
            abs(ht-ht_a)/(ht+ht_a) < tolerance
            and abs(pt-pt_a)/(pt+pt_a) < tolerance
        )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1

def define_calc_t_bounds(ht, p, v, vap):
    if p < p_crit:
        t_sat = sat_t(p)
        if vap:
            ht_vap = calc_ht(t_sat+1, p, v)
            if ht < ht_vap:
                raise ValueError("Enthalpy too low for vapour")
            tmin = t_sat + 1
            tmax = 2000
            t1 = t_sat + 10
        else:
            ht_l = calc_ht(t_sat-1, p, v)
            if ht > ht_l:
                raise ValueError("Enthalpy too high for liquid")
            tmax = t_sat - 1
            tmin = 16
            t1 = 20
                
    else:
        tmin = 16
        tmax = 2000
        t1 = 40
    return t1, tmax, tmin

def loop_check_bounds(t ,tmax, tmin):
    """Ensure that total temperature remains within specified limits"""
    if t > tmax:
        print("warning tmax")
        t = tmax - 1
    if t < tmin:
        print("warning tmin")
        t = tmin + 1
    return t

def sat_t(p):
    """calculate saturation temperature (K) given static pressure (Pa)"""
    if p > p_crit:
        return t_crit
    # antoine equation
    A, B, C = antoine()
    t = B/(A-math.log(p, 10))-C                                # K
    return t

def sat_p(t):
    """calculate saturation pressure (Pa) given static temperature (K)"""
    if t > t_crit:
        return p_crit
    # antoine equation
    A, B, C = antoine()
    p = 10**(A-(B/(t + C)))                                    # Pa
    return p

def calc_ht(t, p, v):
    ht = h2.calc_H2_enthalpy(t, p) + v**2/2
    return ht
        
        
        
        
        