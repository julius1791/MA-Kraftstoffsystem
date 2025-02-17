# -*- coding: utf-8 -*-
import jeta_properties as jeta
import math

# modelling parameters
tolerance = 1e-3    # -
max_iter = 50       # -

def saturation_exception(sat, sat_exp):
    """function for generating exception warning of incompatible saturation"""
    message = "Unexpected Saturation. Saturation is " 
    message += str(sat) + ", expected " + str(sat_exp)
    return message

###############################################################################
######## Parent class for both hydrogen and jeta fuel flow instances ##########
###############################################################################

class JetaFlow:
    def __init__(self, mass_flow: float, temperature: float, pressure: float, velocity: float):
        """
        FuelFlow constructor

        Parameters
        ----------
        mass_flow : float
            mass flow of h2 (kg/s)
        temperature : float
            total temperature (K)
        pressure : float
            total pressure (Pa)
        velocity : float
            velocity of fluid flow (m/s)
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure   
        self.v = velocity
    
    def reduce_pressure(self, p1):      
        """
        Adiabatic pressure reduction (for instance in a throttle)

        Parameters
        ----------
        pressure_ratio : float
            total pressure ratio across the component being modelled (Pa)

        Returns
        -------
        t1 : float
            final total temperature (K)
        """
        if p1 > self.p:
            message = "The selected pressure " + str(p1)
            message += " is greater than initial pressure."
            raise ValueError(message)
        h0 =  self.calc_h()     
        h1 = h0                                   
        p1 = p1
        self.p = p1
        t1 = self.raise_to_h(h1)          
        self.t = t1
        return t1              
        
    def heat_fixed_power(self, Q_dot: float):
        """
        Method for heat added in heat exchanger

        Parameters
        ----------
        Q_dot : float
            absolute thermal power added (W)

        Returns
        -------
        t1 : float
            final total temperature (K)
        """
        # this method may be used regardless of saturation
        # specific heat
        q = Q_dot/self.qm                                               # J/kg
        # initial specific total enthalpy
        h0 = self.calc_h()                                              # J/kg
        # final specific total enthalpy
        h1 = h0 + q                                                     # J/kg
        
        self.raise_to_h(h1)
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
        # calculate initial specific entropy and total enthalpy
        s0 = self.calc_s()
        h0 = self.calc_h()
        # initial copy of flow instance to calculate ideal process
        sim = self.copy()
        # find static temperature that yields the same entropy at higher total pressure
        sim.p = p1
        ts1, _ = sim.find_t_for_s(s0)
        sim.t = ts1
        # find the specific total enthalpy of the isentropic point 
        hs1 = sim.calc_h()
        # find reversible pump work
        a_rev = hs1 - h0 
        # calculate actual compression work and final total enthalpy
        a = a_rev/eta
        h1 = h0 + a
        # find final total pressure and temperature
        self.p = p1
        t1 = self.raise_to_h(h1)
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
        # this method may be applied regardless of saturation states
        # find initial total enthalpy 
        H0_1 = self.calc_h()*self.qm                                    # J
        H0_2 = secondary_flow.calc_h()*secondary_flow.qm                # J
        # final total enthalpy
        H1 = H0_1 + H0_2                                                # J
        # final mass flow
        qm1 = (self.qm+secondary_flow.qm)                               # kg/s
        # final specific total enthalpy
        h1 = H1/qm1                                                     # J/kg
        # final total pressure
        p1 = min(self.p, secondary_flow.p)
        self.p = p1
        self.qm = qm1
        t1 = self.raise_to_h(h1)       
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
    
    def sat_t(self):
        """calculate saturation temperature (K) given the instance pressure"""
        # calculate static temperature and pressure
        _, p = self.calc_static()
        # antoine equation
        A, B, C = self.antoine()
        t = B/(A-math.log(p, 10))-C                                # K
        return t
    
    def sat_p(self):
        """calculate saturation pressure (Pa) given the instance temperature"""
        # calculate static temperature and pressure
        t, _ = self.calc_static()
        # antoine equation
        A, B, C = self.antoine()
        p = 10**(A-(B/(t + C)))                                    # Pa
        return p
    
    

    
    def calc_static(self, tolerance=tolerance, max_iter=max_iter):
        """"calculate static temperature and pressure of fuel flow"""
        ts1 = self.t
        ps1 = self.p
        condition_bool = True
        i = 0
        while condition_bool > tolerance:
            i+=1
            ts0 = ts1
            ts1 = ts1*0.1 + 0.9*(self.t - self.v**2/(2*jeta.calc_jeta_cp(ts0, ps1)))
            ps1 = self.p - jeta.calc_jeta_density(ts0, ps1)*self.v**2/2
            condition_bool = abs(ts0-ts1) > tolerance
            if i > max_iter:
                raise ValueError("calc_static did not converge")
        return ts1, ps1
    
    def antoine(self):
        """antoine equation parameters for jet a"""
        A = 8.81923182000836
        B = 1374.12563
        C = -43.54
        return A, B, C
    
    def calc_h(self):
        """calculate specific total enthalpy of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        h = jeta.calc_jeta_enthalpy(t, p) + self.v**2/2                 # J/kg
        return h
    
    def calc_cp(self):
        """calculate specific heat capacity of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        cp = jeta.calc_jeta_cp(t, p)                                    # J/kgK
        return cp

    def calc_s(self):
        """calculate specific entropy of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        s = jeta.calc_jeta_entropy(t, p)
        return s
    
    def calc_rho(self):
        """calculate specific density of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        rho = jeta.calc_jeta_density(t, p)
        return rho
    
    def raise_to_h(self, h, tolerance=tolerance, max_iter=max_iter, verbosity=False):
        """
        Increase enthalpy of instance to given value by calculating the 
        corresponding temperature

        Parameters
        ----------
        h : float
            target specific enthalpy (J/kg)

        Returns
        -------
        float
            final temperature (K)
        """
        # loop until target enthalpy is achieved
        i = 0
        while abs(self.calc_h()-h) > tolerance:
            i += 1
            # calculate total temperature step asssuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()
            self.t = self.t + dt
            if i > max_iter:
                raise ValueError("raise_to_h did not converge")
        if verbosity:
            print("raise_to_h,", i)
        return self.t
    
    def find_t_for_s(self, s, tolerance=tolerance, max_iter=max_iter, verbosity=False):
        """
        Find the temperature corresponding to a given specific entropy and
        instance pressure

        Parameters
        ----------
        s : float
            specific entropy (J/kgK)

        Returns
        -------
        ts : float
            temperature (K)
        """
        sim = self
        # loop until target entropy is achieved
        i = 0
        while abs(sim.calc_s()-s) > tolerance:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            sim.t = sim.t/math.exp((sim.calc_s()-s)/sim.calc_cp())
            if i > max_iter:
                raise ValueError("find_t_for_s did not converge")
        if verbosity:
            print("find_t_for_s,", i)
        return sim.t, 0
    
    def copy(self):
        return JetaFlow(self.qm, self.t, self.p, self.v)