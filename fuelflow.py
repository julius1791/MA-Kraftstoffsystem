import h2_properties as h2
import jeta_properties as jeta
import math

verbosity = False

def saturation_exception(sat, sat_exp):
    """function for generating exception warning of incompatible saturation"""
    message = "Unexpected Saturation. Saturation is " 
    message += str(sat) + ", expected " + str(sat_exp)
    return message

###############################################################################
######## Parent class for both hydrogen and jeta fuel flow instances ##########
###############################################################################

class FuelFlow:
    def __init__(self, mass_flow: float, temperature: float, pressure: float):
        """
        FuelFlow constructor

        Parameters
        ----------
        mass_flow : float
            mass flow of h2 (kg/s)
        temperature : float
            temperature (K)
        pressure : float
            pressure (Pa)
        saturation : float
            physical state of hydrogen 0 -> super cooled liquid, 
            1 -> overheated gas/supercritical, 0~1 -> saturation state
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure   
    
    def reduce_pressure(self, pressure_ratio):      
        """
        Adiabatic pressure reduction (for instance in a HX or throttle)

        Parameters
        ----------
        pressure_ratio : float
            Pressure ratio across the component being modelled

        Returns
        -------
        t1 : float
            Temperature after pressure drop
        """
        if pressure_ratio > 1:
            message = "The selected pressure ratio " + str(pressure_ratio)
            message += " is greater than 1. Expected a value smaller than 1"
            raise ValueError(message)
        h0 =  self.calc_h()     
        h1 = h0                                   
        p0 = self.p
        p1 = p0 * pressure_ratio   
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
            final temperature (K)
        """
        # this method may be used regardless of saturation
        # specific heat
        q = Q_dot/self.qm                                               # J/kg
        # initial specific enthalpy
        h0 = self.calc_h()                                              # J/kg
        # final enthalpy
        h1 = h0 + q                                                     # J/kg
        self.raise_to_h(h1)
        return self.t
    
    def pump_hydraulic(self, p1: float, eta: float):
        """
        Method for raising pressure of fuel flow in a hydraulic pump

        Parameters
        ----------
        p1 : float
            final pressure (Pa)
        eta : float
            pump efficiency 

        Returns
        -------
        P : float
            pump power (W)
        t1 : float
            final temperature (K)
        """
        # a hydraulic pump may only be used on supercooled liquids
        if type(self) == H2Flow:
            if self.sat != 0:
                raise ValueError(saturation_exception(self.sat, 0))
        # initial density and specific enthalpy
        rho = self.calc_rho()
        h0 = self.calc_h()
        # reversible pump work
        a_rev = (p1-self.p)/rho
        # pump work
        a = a_rev/eta
        # final enthalpy and temperature
        h1 = a + h0
        self.p = p1
        t1 = self.raise_to_h(h1)
        # pump power
        P = a*self.qm
        return P, t1
    
    def mix_flows(self, secondary_flow):
        """
        Method for mixing two hydrogen flows into one

        Parameters
        ----------
        secondary_flow : H2_Flow
            secondary hydrogen flow 

        Returns
        -------
        float
            final temperature (K)
        """
        # this method may be applied regardless of saturation states
        # find initial enthalpy 
        H0_1 = self.calc_h()*self.qm                                    # J
        H0_2 = secondary_flow.calc_h()*secondary_flow.qm                # J
        # calculate final enthalpy
        H1 = H0_1 + H0_2                                                # J
        # final mass flow
        qm1 = (self.qm+secondary_flow.qm)                               # kg/s
        # final specific enthalpy
        h1 = H1/qm1                                                     # J/kg
        # final pressure
        p1 = min(self.p, secondary_flow.p)                              # Pa
        # final temperature
        t1 = self.raise_to_h(h1)                                        # K
        self.p = p1
        self.qm = qm1
        return t1
    
    def split_flows(self, ratio: float):
        """
        Method for splitting one flow instance into two separate flows

        Parameters
        ----------
        ratio : float
            ratio between mass flow of initial flow and secondary outflow

        Returns
        -------
        FuelFlow
            secondary flow split off from main flow
        """
        # splitting flows is possible regardless of saturation state, 
        # however in the saturation regime it is assumed that both  split flows
        # receive the same ratio of liquid and gas
        qm0 = self.qm                                                   # kg/s
        # mass flow of the split flows
        qm1 = qm0*(1-ratio)                                             # kg/s
        qm2 = qm0*ratio                                                 # kg/s
        self.qm = qm1
        # secondary flow inherits pressure and temperature from primary flow
        if type(self) == H2Flow:
            return H2Flow(qm2, self.t, self.p, self.sat)
        else:
            return JetaFlow(qm2, self.t, self.p)
    
    def sat_t(self):
        """calculate saturation temperature (K) given the instance pressure"""
        A, B, C = self.antoine()
        t = B/(A-math.log(self.p, 10))-C                                # K
        return t
    
    def sat_p(self):
        """calculate saturation pressure (Pa) given the instance temperature"""
        A, B, C = self.antoine()
        p = 10**(A-(B/(self.t + C)))                                    # Pa
        return p
    
    
        
    ###########################################################################
    ######### Dummy methods that are populated in the child classes ###########
    ###########################################################################
   
    def antoine(self):
        return
    def calc_h(self):
        return
    def calc_cp(self):
        return
    def calc_s(self):
        return
    def calc_rho(self):
        return
    def raise_to_h(self, h):
        return
    def find_t_for_s(self, s):
        return
    
###############################################################################
################## Child class for hydrogen fuel flow #########################
###############################################################################

class H2Flow(FuelFlow):
    def __init__(
            self, mass_flow: float, temperature: float, 
            pressure: float, saturation: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure)
        # hydrogen requires the additional saturation attribute due to it being
        # present in both liquid and gaseous form throughout the fuel system
        self.sat = saturation
    
    def antoine(self):
        """antoine equation parameters for hydrogen"""
        A = 8.54314
        B = 99.395
        C = 7.726
        return A, B, C
    
    def calc_h(self):
        """calculate specific enthalpy of instance"""
        if self.sat == 0 or self.sat == 1:
            h = h2.calc_H2_enthalpy(self.t, self.p)                     # J/kg
        else:
            h_v = h2.calc_H2_enthalpy(self.t + 0.1, self.p)             # J/kg
            h_l = h2.calc_H2_enthalpy(self.t - 0.1, self.p)             # J/kg
            h = self.sat*h_v + (1-self.sat)*h_l                         # J/kg
        return h
    
    def calc_cp(self):
        """calculate specific isobaric heat capacity of instance"""
        return h2.calc_H2_cp(self.t, self.p)                            # J/kgK
    
    def calc_s(self):
        """calculate specific entropy of instance"""
        if self.sat == 0 or self.sat == 1:
            s = h2.calc_H2_entropy(self.t, self.p)                      # J/kgK
        else:
            s_v = h2.calc_H2_entropy(self.t + 0.1, self.p)              # J/kgK
            s_l = h2.calc_H2_entropy(self.t - 0.1, self.p)              # J/kgK
            s = self.sat*s_v + (1-self.sat)*s_l                         # J/kgK
        return s
    
    def calc_rho(self):
        """calculate density of instance"""
        if self.sat == 0 or self.sat == 1:
            rho = h2.calc_H2_density(self.t, self.p)                    # kg/m3
        else:
            rho_v = h2.calc_H2_density(self.t + 0.1, self.p)            # kg/m3
            rho_l = h2.calc_H2_density(self.t - 0.1, self.p)            # kg/m3
            rho = self.sat*rho_v + (1-self.sat)*rho_l                   # kg/m3
        return rho
    
    def raise_to_h(self, h):
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
        p_crit = 1.2964e6                                               # Pa
        
        # find saturation temperature
        t_sat = self.sat_t()                                            # K
        
        # initialise copies of instance in superheated and supercooled form
        vapour = self
        liquid = self
        liquid.t = t_sat - 0.1
        liquid.sat = 0
        vapour.t - t_sat + 0.1
        vapour.sat = 1
        
        # set solver bounds and initial values depending on physical state
        # of instance after raising enthalpy
        
        # supercritical fluid 
        if self.p > p_crit:
            self.t = 300
            tmax = 2000
            tmin = 20
            self.sat = 1
        else:
            # super heated vapour
            if h > vapour.calc_h():
                tmin = t_sat + 0.1
                self.t = tmin + 10
                tmax = 2000
                self.sat = 1
            # saturation regime
            elif h > liquid.calc_h():
                h_l = liquid.calc_h()
                h_v = vapour.calc_h()
                self.sat = (h - h_l)/(h_v-h_l)
                self.t = t_sat
                return self.t
            # supercooled liquid 
            else:
                self.sat = 0
                tmax = t_sat - 0.1
                self.t = tmax - 5
                tmin = 1
        
        # loop until target enthalpy is achieved
        i = 0
        while abs(self.calc_h()-h) > 1e-6:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()
            self.t = self.t + dt

            # ensure that the new calculated temperature is within the bounds
            # for the current state
            if self.t > tmax:
                self.t = tmax
            if self.t < tmin:
                self.t = tmin
        if verbosity:
            print("raise_to_h,", i)
        return self.t
    
    def find_t_for_s(self, s):
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
        sat_s : float
            saturation parameter in final state
        """
        p_crit = 1.2964e6
        sim = self
        
        # find saturation temperature
        t_sat = self.sat_t()
        
        # initialise copies of instance in superheated and supercooled form
        vapour = self
        liquid = self
        liquid.t = t_sat - 0.1
        liquid.sat = 0
        vapour.t - t_sat + 0.1
        vapour.sat = 1
        
        # set solver bounds and initial values depending on physical state
        # of instance after raising enthalpy
        
        # supercritical fluid 
        if sim.p > p_crit:
            sim.t = t_sat + 10
            tmax = 2000
            tmin = 20
            sim.sat = 1
        else:
            # super heated vapour
            if s > vapour.calc_s():
                tmin = t_sat + 0.1
                sim.t = tmin + 10
                tmax = 2000
                sim.sat = 1
            # saturation regime
            elif s > liquid.calc_s():
                s_l = liquid.calc_s()
                s_v = vapour.calc_s()
                sim.sat = (s - s_l)/(s_v-s_l)
                sim.t = t_sat
                return sim.t, sim.sat
            # supercooled liquid 
            else:
                sim.sat = 0
                tmax = t_sat - 0.1
                sim.t = tmax - 5
                tmin = 1
        
        # loop until target entropy is achieved
        i = 0
        while abs(sim.calc_s()-s) > 1e-9:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            sim.t = sim.t/math.exp((sim.calc_s()-s)/sim.calc_cp())
            
            # ensure that the new calculated temperature is within the bounds
            # for the current state
            if sim.t > tmax:
                sim.t = tmax
            if sim.t < tmin:
                sim.t = tmin
        if verbosity:
            print("find_t_for_s,", i)
        return sim.t, sim.sat
    
    def compress_thermic(self, p1: float, eta: float):
        """
        Compress supercritical or gaseous hydrogen and determine
        the final temperature and power required 

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
        if self.sat != 1:
            raise ValueError(saturation_exception(self.sat, 1))
        # calculate initial specific entropy and enthalpy
        s0 = self.calc_s()
        h0 = self.calc_h()
        sim = self
        # find temperature that yields the same entropy at higher pressur
        sim.p = p1
        ts1, sat_s1 = sim.find_t_for_s(s0)
        sim.t = ts1
        sim.sat = sat_s1
        # find the specific enthalpy of the isentropic point 
        hs1 = sim.calc_h()
        # find reversible compressor work
        a_rev = hs1 - h0 
        # calculate actual compression work and final enthalpy
        a = a_rev/eta
        h1 = h0 + a
        # find final temperature
        self.p = p1
        t1 = self.raise_to_h(h1)
        # calculate compressor power
        P = a*self.qm
        return P, t1
    
    def heat_to_saturation(self):
        """
        Method for adding heat to fuel in heat exchanger to full saturation

        Returns
        -------
        Q_dot : float
            thermal power required for vaporisation (W)

        """
        if self.sat == 1:
            raise ValueError(saturation_exception(self.sat, 0))
        
        h0 = self.calc_h()                                              # J/kg
        t1 = self.sat_t() + 0.1                                         # K
        self.t = t1
        h1 = self.calc_h()                                              # J/kg
        Q_dot = (h1 - h0)*self.qm                                       # W
        self.sat  = 1
        return Q_dot

###############################################################################
#################### Child class for jet a fuel flow ##########################
###############################################################################

class JetaFlow(FuelFlow):
    def __init__(self, mass_flow: float, temperature: float, pressure: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure)
    
    def antoine(self):
        """antoine equation parameters for jet a"""
        A = 8.81923182000836
        B = 1374.12563
        C = 502.76012
        return A, B, C
    
    def calc_h(self):
        """calculate specific enthalpy of instance"""
        return jeta.calc_jeta_enthalpy(self.t, self.p)
    
    def calc_cp(self):
        """calculate specific heat capacity of instance"""
        return jeta.calc_jeta_cp(self.t, self.p)

    def calc_s(self):
        """calculate specific entropy of instance"""
        return jeta.calc_jeta_entropy(self.t, self.p)
    
    def calc_rho(self):
        """calculate specific density of instance"""
        return jeta.calc_jeta_density(self.t, self.p)
    
    def raise_to_h(self, h):
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
        while abs(self.calc_h()-h) > 1e-9:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()
            self.t = self.t + dt
        if verbosity:
            print("raise_to_h,", i)
        return self.t
    
    def find_t_for_s(self, s):
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
        while abs(sim.calc_s()-s) > 1e-9:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            sim.t = sim.t/math.exp((sim.calc_s()-s)/sim.calc_cp())
        if verbosity:
            print("find_t_for_s,", i)
        return sim.t