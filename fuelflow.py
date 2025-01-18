import h2_properties as h2
import jeta_properties as jeta
import math

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
            physical state of hydrogen 0 -> super cooled liquid, 1 -> overheated gas, 0~1 -> saturation state
        """
        self.qm = mass_flow
        self.t = temperature
        self.p = pressure
        
        
    def heat_fixed_power(self, Q_dot: float):
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
        h0 = self.calc_h()
        # final enthalpy
        h1 = h0 + q
        self.raise_to_h(h1)
        return self.t
    
    def heat_to_saturation(self):
        """
        Calculate the power needed to vaporise the fuel

        Returns
        -------
        Q_dot : float
            Thermal power required for vaporisation

        """
        if self.sat == 1:
            raise Exception("Unexpected Saturation. Saturation is ", self.sat, ", expected <1")
        
        h0 = self.calc_h()
        t1 = self.sat_t() + 0.1
        self.t = t1
        h1 = self.calc_h()
        Q_dot = (h1 - h0)*self.qm
        self.sat  = 1
        return Q_dot
    
    def compress_thermic(self, p1: float, eta: float):
        """
        Calculate the power needed to compress supercritical or gaseous hydrogen 
        and the final temperature accounting for non ideal compression

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
            raise Exception("Unexpected Saturation. Saturation is ", self.sat, ", expected 1")
        # calculate initial specific entropy and enthalpy
        s0 = self.calc_s()
        h0 = self.calc_h()
        sim = self
        # find temperature that yields the same entropy at higher pressure (isentropic compression)
        sim.p = p1
        ts1 = sim.find_t_for_s(s0)
        sim.t = ts1
        # find the specific enthalpy of the isentropic point and the reversible compression work
        hs1 = sim.calc_h()
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
    
    def pump_hydraulic(self, p1: float, eta: float):
        """
        Method for raising pressure of fuel flow in a hydraulic pump

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
        if type(self) == H2Flow:
            if self.sat != 0:
                raise Exception("Unexpected Saturation. Saturation is ", self.sat, ", expected 0")
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
        t1 : float
            final temperature (K)

        """
        # find initial enthalpy 
        H0_1 = self.calc_h()*self.qm
        H0_2 = secondary_flow.calc_h()*secondary_flow.qm
        
        # calculate final enthalpy
        H1 = H0_1 + H0_2
        # final mass flow
        qm1 = (self.qm+secondary_flow.qm)
        # final specific enthalpy
        h1 = H1/qm1
        # final pressure
        p1 = min(self.p, secondary_flow.p)
        # final temperature
        t1 = self.raise_to_h(h1)
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
        secondary_flow : FuelFlow
            secondary flow split off from main flow

        """
        qm0 = self.qm
        # mass flow of the split flows
        qm1 = qm0*(1-ratio)
        qm2 = qm0*ratio
        self.qm = qm1
        # secondary flow inherits pressure and temperature from primary flow
        if type(self) == H2Flow:
            return H2Flow(qm2, self.t, self.p, self.sat)
        else:
            return JetaFlow(qm2, self.t, self.p)
    def calc_h(self):
        return 0
    def calc_s(self):
        return
    def calc_cp(self):
        return
    def calc_rho(self):
        return
    def raise_to_h(self, h):
        return self.t
    def sat_t(self):
        A, B, C = self.antoine()
        t = B/(A-math.log(self.p, 10))-C
        return t
    def sat_p(self):
        A, B, C = self.antoine()
        p = 10**(A-(B/(self.t + C)))
        return p
    def antoine(self):
        return 0, 0, 0

class H2Flow(FuelFlow):
    def __init__(self, mass_flow: float, temperature: float, pressure: float, saturation: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure)
        self.sat = saturation
    def antoine(self):
        A = 8.54314
        B = 99.395
        C = 7.726
        return A, B, C
    def calc_h(self):
        if self.sat == 0 or self.sat == 1:
            h = h2.calc_H2_enthalpy(self.t, self.p)
        else:
            h_v = h2.calc_H2_enthalpy(self.t + 0.1, self.p)
            h_l = h2.calc_H2_enthalpy(self.t - 0.1, self.p)
            h = self.sat*h_v + (1-self.sat)*h_l
        return h
    def calc_cp(self):
        return h2.calc_H2_cp(self.t, self.p)
    def calc_s(self):
        if self.sat == 0 or self.sat == 1:
            s = h2.calc_H2_entropy(self.t, self.p)
        else:
            s_v = h2.calc_H2_entropy(self.t + 0.1, self.p)
            s_l = h2.calc_H2_entropy(self.t - 0.1, self.p)
            s = self.sat*s_v + (1-self.sat)*s_l
        return s
    def calc_rho(self):
        if self.sat == 0 or self.sat == 1:
            rho = h2.calc_H2_density(self.t, self.p)
        else:
            rho_v = h2.calc_H2_density(self.t + 0.1, self.p)
            rho_l = h2.calc_H2_density(self.t - 0.1, self.p)
            rho = self.sat*rho_v + (1-self.sat)*rho_l
        return rho
    def raise_to_h(self, h):
        # find saturation temperature
        t_sat = self.sat_t()
        p_crit = 1.2964e6
        vapour = self
        liquid = self
        liquid.t = t_sat - 0.1
        liquid.sat = 0
        vapour.t - t_sat + 0.1
        vapour.sat = 1
        if self.p > p_crit:
            self.t = 300
            tmax = 2000
            tmin = 20
            self.sat = 1
        else:
            # determine bounds and initial temperature depending on the physical state
            # determine if the hydrogen is a superheated gas
            if h > vapour.calc_h():
                tmin = t_sat + 0.1
                self.t = tmin + 10
                tmax = 2000
                self.sat = 1
            # determine if the hydrogen is in the saturation regime
            elif h > liquid.calc_h():
                h_l = liquid.calc_h()
                h_v = vapour.calc_h()
                self.sat = (h - h_l)/(h_v-h_l)
                self.t = t_sat
                return self.t
            # else it is a supercooled liquid 
            else:
                self.sat = 0
                tmax = t_sat - 0.1
                self.t = tmax - 5
                tmin = 1
        
        # loop until target enthalpy is achieved
        while abs(self.calc_h()-h) > 1e-6:
            # calculate temperature step asssuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()

            self.t = self.t + dt

            # ensure that the new calculated temperature is within the bounds
            # for the current state
            if self.t > tmax:
                self.t = tmax
            if self.t < tmin:
                self.t = tmin
        return self.t
    def find_t_for_s(self, s):
        sim = self
        # loop until target entropy is achieved
        while abs(sim.calc_s()-s) > 1e-9:
            # calculate temperature step asssuming constant heat capacity
            sim.t = sim.t/math.exp((sim.calc_s()-s)/sim.calc_cp())
        
        return sim.t
    
class JetaFlow(FuelFlow):
    def __init__(self, mass_flow: float, temperature: float, pressure: float):
        FuelFlow.__init__(self, mass_flow, temperature, pressure)
    def antoine(self):
        A = 8.81923182000836
        B = 1374.12563
        C = 502.76012
        return A, B, C
    def calc_cp(self):
        return jeta.calc_jeta_cp(self.t, self.p)
    def calc_h(self):
        return jeta.calc_jeta_enthalpy(self.t, self.p)
    def calc_s(self):
        return jeta.calc_jeta_entropy(self.t, self.p)
    def calc_rho(self):
        return jeta.calc_jeta_density(self.t, self.p)
    def raise_to_h(self, h):
        # loop until target enthalpy is achieved
        while abs(self.calc_h()-h) > 1e-9:
            # calculate temperature step asssuming constant heat capacity
            dt = (h-self.calc_h())/self.calc_cp()
            self.t = self.t + dt
        return self.t
