import PsychroPy.psychropy as psp
import h2_properties as h2

tolerance = 1e-6
rel_fac = 0.9


# hydrogen combustion mass fractions
w_h2 = 2.02            # g/mol
w_o2 = 16              # g/mol
w_h2o = 18.02          # g/mol

# mass fractions of dry air
w_n2_air = 0.7552
w_o2_air = 0.2314
w_ar_air = 0.0129
w_co2_air = 0.000626

# specific humidity of inlet bleed air
w_h2o_air_0 = 0.000000

# lower heating value (25˚C)
lhv_h2 = 119.96e6       # J/kg
t_ref = 298.15          # K

def calc_parallel_combustion(phi, t_h2_0, t_air_0, p_h2, p_air, Q):
    """
    Calculate the mass flow of hydrogen to a bleed air supplied 
    parallel combustion chamber required to supply the set amount of heat  

    Parameters
    ----------
    phi : float
        fuel - air equivalence ratio (-)
    t_h2_0 : float
        static temperature of hydrogen prior to injection (K)
    t_air_0 : float
        static temperature of inlet bleed air (K)
    p_h2 : float
        static pressure of hydrogen prior to injection (Pa)
    p_air : float
        static pressure of bleed air (Pa)
    Q : float
        Heat demand (W)

    Returns
    -------
    qm_h2 : float
        hydrogen mass flow (kg/s)

    """
    # calculate ratio of dry air to hydrogen given the equivalence ratio before combustion
    w_air_h2_0 = w_o2/w_h2 / (phi*w_o2_air)
    
    ###########################################################################
    
    # only required if combustion gases are cooled to a 
    # temperature that isn't the reference temperature
    
    # # calculate 
    # w_air_h2_1 = w_o2/w_h2 * (1/(phi*w_o2_air)-1)
    
    # # calculate ratio of water produced in combustion to hydrogen
    # w_h2o_h2 = w_h2o/w_h2
    
    # # calculate specific humidity after combustion
    # w_h2o_air_1 = w_h2o_h2/w_air_h2_1 + w_h2o_air*w_air_h2_0/w_air_h2_1 
    
    ###########################################################################
    
    # calculate premix gas temperature
    t0 = calculate_initial_temperature(w_air_h2_0, t_h2_0, t_air_0, p_air, p_h2)
    
    # calculate the mass flow of hydrogen required
    qm_h2 = calculate_enthalpy_change(w_air_h2_0, w_h2o_air_0, t0, p_air, Q)
    return qm_h2

def calculate_initial_temperature(w_air_h2, t_h2_0, t_air_0, p_air, p_h2):
    """
    Calculate the initial temperature of the reactant gases after
    eliquibrium has been established

    Parameters
    ----------
    w_air_h2 : float
        mass ratio of dry air to hydrogen (-)
    t_h2_0 : float
        static temperature of hydrogen prior to injection (K)
    t_air_0 : float
        static temperature of inlet bleed air (K)
    p_h2 : float
        static pressure of hydrogen prior to injection (Pa)
    p_air : float
        static pressure of bleed air (Pa)

    Returns
    -------
    t1 : float
        Initial temperature in combustion chamber (K)

    """
    # arbitrary initial temperature
    t1 = 300    
    
    # init humid air instances
    air0 = psp.MoistAir()
    air1 = psp.MoistAir()
    
    condition_bool = True
    while condition_bool:
        # calculate change in specific enthalpy of hydrogen
        dh_h2 = h2.calc_H2_enthalpy(t_h2_0, p_h2) - h2.calc_H2_enthalpy(t1, p_air)
        
        # calculate change in specific enthalpy of air
        air0.TW = t_air_0, w_h2o_air_0
        air1.TW = t1, w_h2o_air_0
        dh_air = (air0.specific_enthalpy - air1.specific_enthalpy)*1e3
        
        # calculate specific heat capacity of air
        cp_air = dh_air/(t_air_0-t1)
        
        # calculate specific enthalpy of the mixture
        dh = (dh_h2 + dh_air * w_air_h2)/(1 + w_air_h2)
        
        # calculate specific heat capacity of the mixture
        cp = (h2.calc_H2_cp(t1, p_air) + cp_air * w_air_h2)/(1 + w_air_h2)
        
        # calculate the next temperature step
        t1 += dh/cp * rel_fac
        
        condition_bool = abs(dh) > tolerance
    return t1

def calculate_enthalpy_change(w_air_h2_0, w_h2o_air_0, t0, p, Q):
    """
    Calculate the specific heat that is transferred in the heat exchanger 
    following the parallel combustion chamber if the exhaust gases are cooled 
    to 25˚C/298.15 K

    Parameters
    ----------
    w_air_h2_0 : float
        mass ratio of dry air [before combustion] to hydrogen (-)
    w_h2o_air_0 : float
        specific humidity of inlet bleed air (-)
    t0 : float
        Initial temperature in combustion chamber (K)
    p : float
        static pressure of bleed air (Pa)
    Q : float
        Heat demand (W)

    Returns
    -------
    qm_h2 : float
        hydrogen mass flow (kg/s)

    """
    # calculate the specific enthalpy difference of hydrogen to the reference point
    dh_h2_0 = h2.calc_H2_enthalpy(t0, p) - h2.calc_H2_enthalpy(t_ref, p)
    
    # calculate the specific enthalpy differenc of bleed air to the reference point
    air0 = psp.MoistAir()
    airref0 = psp.MoistAir()
    air0.TW = t0, w_h2o_air_0
    airref0.TW = t_ref, w_h2o_air_0
    dh_air_0 = (air0.specific_enthalpy-airref0.specific_enthalpy)*1e3
    
    # calculate the specific enthalpy difference of the mixture
    dh0 = (dh_h2_0+w_air_h2_0 * dh_air_0)/(1+w_air_h2_0)
    
    # calculate the [mixture] specific heat transferred in the heat exchanger
    q = dh0 + lhv_h2/(1+w_air_h2_0)
    # calculate the [hydrogen] specific heat transferred in the heat exchanger
    q_h2 = q*(1+w_air_h2_0)
    # calculate the mass of hydrogen needed to supply the heat demand
    qm_h2 = Q/q_h2

    return qm_h2

if __name__ == "__main__":
    print(calc_parallel_combustion(1/10, 250, 344, 1.6e6, 1e5, 150e3))