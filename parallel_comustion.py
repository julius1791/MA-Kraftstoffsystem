# -*- coding: utf-8 -*-

# heat capacities
cp_o2=919.63           # J/kgK
cp_n2=1041.3          # J/kgK
cp_h2o=1911.8         # J/kgK
cp_h2=14858         # J/kgK

# mass numbers
Mr_o2=31.9988           # g/mol
Mr_h2=2.01588	        # g/mol
Mr_h2o=18.0153          # g/mol

# mass ratio of oxygen in air
w_o2_air = 0.231781

# lower heating value of hydrogen
lhv_h2 = 119.96e6       # J/kg
t_ref = 298.15          # K


# stoichometric mass ratio of hydrogen and oxygen
mr_o2_h2 = Mr_o2/(2*Mr_h2)

# ratio of water produced
mr_h2o_h2 = Mr_h2o/Mr_h2

def parallel_combustion(Q, t_h2, t_air=272.63, t_hx=400, phi=0.3):
    """
    Calculate the mass flow of hydrogen to a bleed air supplied 
    parallel combustion chamber required to supply the set amount of heat  

    Parameters
    ----------
    phi : float
        fuel - air equivalence ratio (-)
    t_h2 : float
        static temperature of hydrogen prior to injection (K)
    t_air : float
        stagnation temperature of inlet bleed air (K)
    t_hx : float
        static temperature of the cooled exhaust gas (K)
    Q : float
        Heat demand (W)

    Returns
    -------
    qm_h2 : float
        hydrogen mass flow (kg/s)

    """
    if phi > 1:
        raise ValueError("Unexpected Phi value > 1")
        
    # ratio unburnt oxygen
    mr_bo2_h2 = mr_o2_h2*(1/phi-1)
    
    # ratio of nitrogen
    mr_n2_h2 = mr_o2_h2*(1-w_o2_air)/(phi*w_o2_air)
    
    # energy balance
    dh2_n2 = mr_n2_h2*cp_n2*(t_air-t_hx)
    dh2_o2b = mr_bo2_h2*cp_o2*(t_air-t_hx)
    dh2_h2 = lhv_h2+cp_h2*(t_h2-t_ref)
    dh2_o2 = mr_o2_h2*cp_o2*(t_air-t_ref)
    dh2_h2o = mr_h2o_h2*cp_h2o*(t_hx-t_ref)

    m_h2 = Q/(dh2_n2+dh2_o2+dh2_o2b+dh2_h2-dh2_h2o)
    
    return m_h2

if __name__ == "__main__":
    print(parallel_combustion(450e3, 400, phi=0.3))
