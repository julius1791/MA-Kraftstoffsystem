# -*- coding: utf-8 -*-
import fluid_properties.jeta_properties as jeta
import math

# modelling parameters
tolerance = 1e-1    # -
max_iter = 100      # -

def antoine():
    """antoine equation parameters for jet a"""
    A = 8.81923182000836
    B = 1374.12563
    C = -43.54
    return A, B, C

###############################################################################
######################### Class for Jet-A fuel flow ###########################
###############################################################################

class JetaFlow:
    def __init__(
            self, mass_flow: float,
            temperature: float, pressure: float, velocity: float
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
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure   
        self.v = velocity
        self.ht = calc_ht(temperature, pressure, velocity)
    
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
        ht0 =  calc_ht(self.t, self.p, self.v) 
        ht1 = ht0    
        self.p = p1
        t1 = calc_t(ht1, p1, self.v)          
        self.t = t1
        return t1              
        
    def heat_exchanger(self, Q_dot: float, tpr: float):
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
            final static temperature (K)
        """
        # this method may be used regardless of saturation
        # specific heat
        pt0 = calc_pt(self.t, self.p, self.v)   
        q = Q_dot/self.qm                                               # J/kg
        # initial specific total enthalpy
        ht0 = calc_ht(self.t, self.p, self.v)                           # J/kg
        # final specific total enthalpy
        ht1 = ht0 + q                                                   # J/kg
        self.ht = ht1
        t1 = calc_t(ht1, self.p, self.v)
        self.t = t1
        
        pt1 = pt0 * tpr
        t1, p1 = apply_total_pressure(pt1, ht1, t1, self.v)
        self.t = t1
        self.p = p1
        
        return self.t
        
    def pump_hydraulic(self, p1: float, eta: float):
        """
        Method for raising pressure of fuel flow in a hydraulic pump

        Parameters
        ----------
        p1 : float
            final total pressure (Pa)
        eta : float
            pump efficiency (-)

        Returns
        -------
        P : float
            pump power (W)
        t1 : float
            final total temperature (K)
        """        
        # calculate initial specific entropy and specific total enthalpy
        s = jeta.calc_jeta_entropy(self.t, self.p)
        ht_0 = self.ht
        # find static temperature that yields
        # the same entropy at higher pressure
        ts1 = find_t_for_s(s, p1, self.t)
        # find the specific total enthalpy of the isentropic point 
        ht_s1 = calc_ht(ts1, p1, self.v)
        # find reversible pump work
        a_rev = ht_s1 - ht_0 
        # calculate actual compression work and final total enthalpy
        a = a_rev/eta
        ht_1 = ht_0 + a
        
        self.p = p1
        self.ht = ht_1
        
        # find final static temperature
        t1 = calc_t(ht_1, p1, self.v)
        self.t = t1
        
        # calculate compressor power
        P = a*self.qm
        
        return P, t1
    
    def mix_flows(self, secondary_flow):
        """
        Method for mixing two fuel flow instances into one

        Parameters
        ----------
        secondary_flow : FuelFlow
            secondary fuel flow 

        Returns
        -------
        float
            final total temperature (K)
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
        # splitting flows is possible regardless of saturation state, 
        # however in the saturation regime it is assumed that both  split flows
        # receive the same ratio of liquid and gas
        if qm > self.qm:
            message = "Can't split " + str(qm) + " kg/s because "
            message += "original flow is only " + str(self.qm) + " kg/s"
            raise ValueError(message)
        # secondary flow inherits pressure and temperature from primary flow
        self.qm = self.qm - qm
        
        return JetaFlow(qm, self.t, self.p, self.v)
    
    
    def calc_static(self, tolerance=tolerance, max_iter=max_iter):
        """"calculate static temperature and pressure of fuel flow"""
        ts1 = self.t
        ps1 = self.p
        condition_bool = True
        i = 0
        while condition_bool:
            i+=1
            ts0 = ts1
            ts1 = [
                ts1*0.1 
                + 0.9*(self.t - self.v**2/(2*jeta.calc_jeta_cp(ts0, ps1)))
            ]
            ps1 = self.p - jeta.calc_jeta_density(ts0, ps1)*self.v**2/2
            condition_bool = abs(ts0-ts1) > tolerance
            if i > max_iter:
                raise ValueError("calc_static did not converge")
        return ts1, ps1

    







    
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
        s_a = jeta.calc_jeta_entropy(t0, p)
        t1 = t0*math.exp((s-s_a)/jeta.calc_jeta_cp(t0, p))
        s_a = jeta.calc_jeta_entropy(t1, p)
        condition_bool = abs(s_a-s) > tolerance
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
        pt_a = p0 + v**2 * jeta.calc_jeta_density(t0, p0) / 1.6
        
        # calculate static pressure and temperature for next iteration
        p1 += pt - pt_a
        t1 += (ht-ht_a)/min(jeta.calc_jeta_cp(t0, p0), 8e3) / 4
        condition_bool = not (
            abs(ht-ht_a) < tolerance
            and abs(pt-pt_a) < tolerance
            )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1, p1

def calc_pt(t, p, v):
    pt = p + jeta.calc_jeta_density(t, p) * v**2/2
    return pt

        
def calc_tp(ht, pt, v):
    t1 = 300
    p1 = pt
    
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1 
        t0 = t1
        p0 = p1
        
        ht_a = calc_ht(t0, p0, v)
        pt_a = p1 + jeta.calc_jeta_density(t0, p0) * v**2/2
        t1 += (ht - ht_a)/(jeta.calc_jeta_cp(t0, p0))
        p1 += pt-pt_a
        
        condition_bool = not (
            abs(ht-ht_a) < tolerance
            and abs(pt-pt_a) < tolerance
        )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    p_sat = sat_p(t1)
    if p_sat > p1:
        raise Exception("Vapour pressure exceeds static pressure")
    return t1, p1

def calc_t(ht, p, v):
    t1 = 300
    condition_bool = True
    i = 0
    while condition_bool:
        i += 1 
        t0 = t1
        
        ht_a = calc_ht(t0, p, v)
        t1 += (ht - ht_a)/(jeta.calc_jeta_cp(t0, p))
        
        condition_bool = not (
            abs(ht-ht_a) < tolerance
        )
        if i > max_iter:
            raise Exception("Exceeded number of iterations")
    return t1

def sat_t(p):
    """calculate saturation temperature (K) given static pressure (Pa)"""
    # antoine equation
    A, B, C = antoine()
    t = B/(A-math.log(p, 10))-C                                # K
    return t

def sat_p(t):
    """calculate saturation pressure (Pa) given static temperature (K)"""
    # antoine equation
    A, B, C = antoine()
    p = 10**(A-(B/(t + C)))                                    # Pa
    return p

def calc_ht(t, p, v):
    ht = jeta.calc_jeta_enthalpy(t, p) + v**2/2
    return ht
        
        
        
        
        