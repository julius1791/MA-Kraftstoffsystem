# -*- coding: utf-8 -*-

# cp of air
cp_air = 1004 # J/kg
# cp of hydrogen combustion products at far = 0.01
cp_b = 1104 # J/kg
# fuel to air ratio
far = 0.01

# lower heating value of hydrogen
lhv_h2 = 119.96e6       # J/kg

def parallel_combustion(Q, t_air=272.63, t_hx=400):
    """
    Calculate the mass flow of hydrogen to a bleed air supplied 
    parallel combustion chamber required to supply the set amount of heat  

    Parameters
    ----------
    t_air : float
        stagnation temperature of inlet bleed air (K)
    t_hx : float
        static temperature of the cooled exhaust gas (K)
    Q : float
        Heat demand (W)

    Returns
    -------
    m_h2 : float
        hydrogen mass flow (kg/s)
    m_z : float
        bleed air mass flow (kg/s)  

    """
    # simplified formula for hydrogen demand
    m_h2 = Q/(lhv_h2+cp_air*t_air/far-(1+1/far)*cp_b*t_hx)
    
    # calculate bleed air demand
    m_z = m_h2/far
    
    return m_h2, m_z

if __name__ == "__main__":
    print(parallel_combustion(450e3, t_hx=300))
