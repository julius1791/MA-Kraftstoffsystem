import h2_properties as h2
import jeta_properties as jeta
import math

t_oil = 350         # K
dT_ref = 105        # K
q_ref = 123105      # J/kg
tpr_ref = 0.01489   # -
tolerance = 1e-6    # -
max_iter = 50       # -
p_crit = 1.2964e6   # Pa


def saturation_exception(sat, sat_exp):
    """function for generating exception warning of incompatible saturation"""
    message = "Unexpected Saturation. Saturation is " 
    message += str(sat) + ", expected " + str(sat_exp)
    return message

###############################################################################
######## Parent class for both hydrogen and jeta fuel flow instances ##########
###############################################################################

class FuelFlow:
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
        saturation : float
            physical state of hydrogen 0 -> super cooled liquid, 
            1 -> overheated gas/supercritical, 0~1 -> saturation state
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
        
        if isinstance(self, H2Flow):
            dT = t_oil - self.t
            tpr = tpr_ref*q/q_ref/dT*dT_ref
            p1 = self.p * (1-tpr)
            self.p = p1
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
        # a hydraulic pump may only be used on supercooled liquids
        if isinstance(self, H2Flow):
            if self.sat != 0:
                raise ValueError(saturation_exception(self.sat, 0))
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
        if isinstance(self, H2Flow):
            return H2Flow(qm, self.t, self.p, self.v, self.sat)
        else:
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
    
    
        
    ###########################################################################
    ######### Dummy methods that are populated in the child classes ###########
    ###########################################################################
   
    def calc_static(self):
        raise NotImplementedError("Use Subclass")
    def antoine(self):
        raise NotImplementedError("Use Subclass")
    def calc_h(self):
        raise NotImplementedError("Use Subclass")
    def calc_s(self):
        raise NotImplementedError("Use Subclass")
    def raise_to_h(self, h):
        raise NotImplementedError("Use Subclass")
    def find_t_for_s(self, s):
        raise NotImplementedError("Use Subclass")
    def copy(self):
        raise NotImplementedError("Use Subclass")
    
    
###############################################################################
################## Child class for hydrogen fuel flow #########################
###############################################################################

class H2Flow(FuelFlow):
    def __init__(
            self, mass_flow: float, temperature: float, 
            pressure: float, velocity: float, saturation: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure, velocity)
        # hydrogen requires the additional saturation attribute due to it being
        # present in both liquid and gaseous form throughout the fuel system
        self.sat = saturation
        
    def calc_static(self, tolerance=tolerance, max_iter=max_iter):
        if self.sat != 1 and self.sat != 0:
            A, B, C = self.antoine()
            condition_bool = True
            i = 0
            ps1 = self.p
            while condition_bool > tolerance:
                ps0 = ps1
                ts1= B/(A-math.log(ps0, 10))-C 
                ps1 = self.p - h2.calc_H2_density(ts1, ps0)*self.v**2/2
                condition_bool = abs(ps0-ps1) > tolerance
                if i > max_iter:
                    raise ValueError("calc_static did not converge")
            return ts1, ps1
        ts1 = self.t
        ps1 = self.p
        condition_bool = True
        i = 0
        while condition_bool > tolerance:
            i+=1
            ts0 = ts1
            ts1 = ts1*0.1 + 0.9*(self.t - self.v**2/(2*h2.calc_H2_cp(ts0, ps1)))
            ps1 = self.p - h2.calc_H2_density(ts0, ps1)*self.v**2/2
            condition_bool = abs(ts0-ts1) > tolerance
            if i > max_iter:
                raise ValueError("calc_static did not converge")
        return ts1, ps1
    
    def antoine(self):
        """antoine equation parameters for hydrogen"""
        A = 8.54314
        B = 99.395
        C = 7.726
        return A, B, C
    
    def calc_h(self):
        """calculate specific total enthalpy of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        if self.sat == 0 or self.sat == 1:
            h = h2.calc_H2_enthalpy(t, p) + self.v**2/2                 # J/kg
        else:
            h_v = h2.calc_H2_enthalpy(t + 0.1, p)                       # J/kg
            h_l = h2.calc_H2_enthalpy(t - 0.1, p)                       # J/kg
            h = self.sat*h_v + (1-self.sat)*h_l + self.v**2/2           # J/kg
        return h
    
    def calc_cp(self):
        """calculate specific isobaric heat capacity of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        return h2.calc_H2_cp(t, p)                                      # J/kgK
    
    def calc_s(self):
        """calculate specific entropy of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        if self.sat == 0 or self.sat == 1:
            s = h2.calc_H2_entropy(t, p)                                # J/kgK
        else:
            s_v = h2.calc_H2_entropy(t + 0.1, p)                        # J/kgK
            s_l = h2.calc_H2_entropy(t - 0.1, p)                        # J/kgK
            s = self.sat*s_v + (1-self.sat)*s_l                         # J/kgK
        return s
    
    def calc_rho(self):
        """calculate density of instance"""
        # calculate static temperature and pressure
        t, p = self.calc_static()
        if self.sat == 0 or self.sat == 1:
            rho = h2.calc_H2_density(t, p)                              # kg/m3
        else:
            rho_v = h2.calc_H2_density(t + 0.1, p)                      # kg/m3
            rho_l = h2.calc_H2_density(t - 0.1, p)                      # kg/m3
            rho = self.sat*rho_v + (1-self.sat)*rho_l                   # kg/m3
        return rho
    
    def raise_to_h(self, h):
        """
        Increase total enthalpy of instance to given value by calculating the 
        corresponding total temperature

        Parameters
        ----------
        h : float
            target specific total enthalpy (J/kg)

        Returns
        -------
        float
            final total temperature (K)
        """
        # set solver bounds and initial values depending on physical state
        # of instance after raising enthalpy
        
        # find saturation temperature and static pressure
        t_sat = self.sat_t()
        _, p = self.calc_static()
        
        # supercritical fluid 
        if self.p > p_crit:
            tmax, tmin = self.init_loop(2000, 15, 300, 1)
            return self.enthalpy_loop(tmax, tmin, h)
        else:
            # initialise copies of instance in superheated and supercooled form
            vapour, liquid = self.init_vapour_liquid()
            # check for super heated vapour
            if h > vapour.calc_h():
                tmin = t_sat + 1 + self.v**2/(2*h2.calc_H2_cp(t_sat, p))
                tmax, tmin = self.init_loop(2000, tmin, tmin + 10, 1)
                return self.enthalpy_loop(tmax, tmin, h)
            # check for saturation 
            elif h > liquid.calc_h():
                h_l = liquid.calc_h()
                h_v = vapour.calc_h()
                self.sat = (h - h_l)/(h_v-h_l)
                self.t = t_sat + self.v**2/(2*h2.calc_H2_cp(t_sat, p))
                return self.t
            # else it must be supercooled liquid 
            else:
                tmax = t_sat - 1 + self.v**2/(2*h2.calc_H2_cp(t_sat, p))
                tmax, tmin = self.init_loop(tmax, 1, tmax - 5, 0)
                return self.enthalpy_loop(tmax, tmin, h)
            
    def init_loop(self, tmax, tmin, t, sat):
        self.t = t
        self.sat = sat
        return tmax, tmin
    
    def enthalpy_loop(self, tmax, tmin, h, tolerance=tolerance, max_iter=max_iter, verbosity=False):
        """Iteratively determine final temperature for given specific total 
        enthalpy and total pressure of supercooled/heated fluid"""
        # loop until target enthalpy is achieved
        i = 0
        while abs(self.calc_h()-h) > tolerance:
            i += 1
            # calculate total temperature step assuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()
            self.t = self.t + dt

            # ensure that the new calculated total temperature is within the 
            # bounds for the current state
            self.loop_check_bounds(tmax, tmin)
            
            # quit loop after max iterations are reached
            if i > max_iter:
                raise ValueError("raise_to_h did not converge")
        if verbosity:
            print("raise_to_h,", i)
        return self.t
    
    def loop_check_bounds(self, tmax, tmin):
        """Ensure that total temperature remains within specified limits"""
        if self.t > tmax:
            self.t = tmax
        if self.t < tmin:
            self.t = tmin
        return
    
    def find_t_for_s(self, s):
        """
        Find the total temperature corresponding to a given specific entropy and
        instance total pressure

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
        # set solver bounds and initial values depending on physical state
        # of instance after raising enthalpy
        
        sim = self.copy()
        # find saturation temperature and static pressure
        t_sat = self.sat_t()
        _, p = self.calc_static()
        
        # supercritical fluid 
        if sim.p > p_crit:
            tmax, tmin = sim.init_loop(2000, 15, t_sat + 10, 1)
            return sim.entropy_loop(tmax, tmin, s)
        else:
            # initialise copies of instance in superheated and supercooled form
            vapour, liquid = self.init_vapour_liquid()
            
            # check for super heated vapour
            if s > vapour.calc_s():
                tmin = t_sat + 1 + self.v**2/(2*h2.calc_H2_cp(t_sat, p))
                tmax, tmin = sim.init_loop(2000, tmin, tmin + 10, 1)
                return sim.entropy_loop(tmax, tmin, s)
            # check for saturation 
            elif s > liquid.calc_s():
                s_l = liquid.calc_s()
                s_v = vapour.calc_s()
                sim.sat = (s - s_l)/(s_v-s_l)
                sim.t = t_sat
                return sim.t, sim.sat
            # else it must be supercooled liquid 
            else:
                tmax = t_sat - 1 + self.v**2/(2*h2.calc_H2_cp(t_sat, p))
                tmax, tmin = sim.init_loop(tmax, 1, tmax - 5, 0)
                return sim.entropy_loop(tmax, tmin, s)
        
        
    def entropy_loop(self, tmax, tmin, s, tolerance=tolerance, max_iter=max_iter, verbosity=False):
        """Iteratively determine final temperature for given specific entropy 
        and total pressure of supercooled/heated fluid within constraints"""
        # loop until target entropy is achieved
        i = 0
        while abs(self.calc_s()-s) > tolerance:
            i += 1
            # calculate temperature step asssuming constant heat capacity
            self.t = self.t/math.exp((self.calc_s()-s)/self.calc_cp())
            
            # ensure that the new calculated temperature is within the bounds
            # for the current state
            self.loop_check_bounds(tmax, tmin)
            
            # quit loop after max iterations are reached
            if i > max_iter:
                raise ValueError("find_t_for_s did not converge")
        if verbosity:
            print("find_t_for_s,", i)
        return self.t, self.sat
    
    def init_vapour_liquid(self):
        """"Initialise H2flow instances of saturated vapour/liquid"""
        # find saturation temperature and static pressure
        t_sat = self.sat_t()
        _, p = self.calc_static()
        
        # initialise copies of instance in superheated and supercooled form
        vapour = self.copy()
        liquid = self.copy()
        liquid.t = t_sat - 1 + self.v**2/(2*h2.calc_H2_cp(t_sat - 1, p))
        liquid.sat = 0
        vapour.t = t_sat + 1 + self.v**2/(2*h2.calc_H2_cp(t_sat + 1, p))
        vapour.sat = 1
        return vapour, liquid
    
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
        _, p = self.calc_static()
        t1 = self.sat_t() + 1 + self.v**2/(2*h2.calc_H2_cp(self.sat_t() + 1, p))  # K
        self.t = t1
        h1 = self.calc_h()                                              # J/kg
        Q_dot = (h1 - h0)*self.qm                                       # W
        self.sat  = 1
        return Q_dot
    
    def copy(self):
        return H2Flow(self.qm, self.t, self.p, self.v, self.sat)


###############################################################################
#################### Child class for jet a fuel flow ##########################
###############################################################################

class JetaFlow(FuelFlow):
    def __init__(self, mass_flow: float, temperature: float, pressure: float, velocity: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure, velocity)
    
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