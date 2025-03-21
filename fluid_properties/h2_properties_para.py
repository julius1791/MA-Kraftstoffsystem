import numpy as np
import math

# A - Helmholtzsche freie Energie A = f(rho, T) [J/kg]
# alpha - reduzierte Energie alpha = A/RT = f(delta, tau)
# delta: reduzierte Dichte delta=rho/rho_c
# tau: reduzierte reziproke Temperatur tau=T_c/T

# alpha = alpha_0(tau, deltha) + alpha_r(tau, deltha)

##################################################

# Kritische Werte
rho_c = 31.076     # [kg/m^3] Kritische Dichte für H2
T_c   = 32.938     # [K]        ; Kritische Temperatur für H2
p_c   = 1.2858e6   # [Pa]       ; Kritischer Druck für H2

R     = 4124.2     # [J/(kg*K)] ; Spezifische Gaskonstante für H2

# Koeffizienten für Normalwasserstoff
para_hydrogen = {
    "a_k": {
        1: -1.4485891134,
        2: 1.884521239,
        3: 4.30256,
        4: 13.0289,
        5: -47.7365,
        6: 50.0013,
        7: -18.6261,
        8: 0.993973,  # Kein Wert gegeben
        9: 0.536078,  # Kein Wert gegeben
    },
    "b_k": {
        1: 0.0,  # Kein Wert gegeben
        2: 0.0,  # Kein Wert gegeben
        3: -499.0/T_c,
        4: -826.5/T_c,
        5: -970.8/T_c,
        6: -1166.2/T_c,
        7: -1341.4/T_c,
        8: -5395.0/T_c,  # Kein Wert gegeben
        9: -10185.0/T_c  # Kein Wert gegeben
    },
    "N_i": {
        1: -7.33375,
        2: 0.01,
        3: 2.60375,
        4: 4.66279,
        5: 0.682390,
        6: -1.47078,
        7: 0.135801,
        8: -1.05327,
        9: 0.328239,
        10: -0.0577833,
        11: 0.0449743,
        12: 0.0703464,
        13: -0.0401766,
        14: 0.119510
    },
    "t_i": {
        1: 0.6855,
        2: 1,
        3: 1,
        4: 0.489,
        5: 0.774,
        6: 1.133,
        7: 1.386,
        8: 1.619,
        9: 1.162,
        10: 3.96,
        11: 5.276 ,
        12: 0.99,
        13: 6.791,
        14: 3.19 
    },
    "d_i": {
        1: 1,
        2: 4,
        3: 1,
        4: 1,
        5: 2,
        6: 2,
        7: 3,
        8: 1,
        9: 3,
        10: 2,
        11: 1,
        12: 3,
        13: 1,
        14: 1
    },
    "p_i": {
        8: 1,
        9: 1
    },
    "phi_i": {
        10: -1.7437,
        11: -0.5516,
        12: -0.0634,
        13: -2.1341,
        14: -1.777
    },
    "beta_i": {
        10: -0.194,
        11: -0.2019,
        12: -0.0301,
        13: -0.2383,
        14: -0.3253
    },
    "gamma_i": {
        10: 0.8048,
        11: 1.5248,
        12: 0.6648,
        13: 0.6832,
        14: 1.493
    },
    "D_i": {
        10: 1.5487,
        11: 0.1785,
        12: 1.28,
        13: 0.6319,
        14: 1.7104
    },
    "b_i_v": {
        1: 1.81758329e-7,
        2: 0.683106758,
        3: 9.32706091e-13,
        4: 1.48078541,
        5: 1.27555239
    },
    "c_i_l": {
        1: 2.34498695e-3,
        2: 0.764814482,
        3: 1.39412767e-8,
        4: 0.851102621,
        5: 1.11721850
    },
    "N_i_s": {
        1: -4.87767,
        2: 1.03359,
        3: 0.82668,
        4: -0.129412
    },
    "k_i" : {
        1: 1,
        2: 1.5,
        3: 2.65,
        4: 7.4
    }
}


############################################

# Berechnung von alpha_0
def alpha_0(delta, tau, para_hydrogen):

    k_sum_3_9 = 0
    for k in range(3, 10):
        if para_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            k_sum_3_9 += para_hydrogen["a_k"].get(k, 0) * np.log(1 - np.exp(
                para_hydrogen["b_k"][k] * tau))

    alpha_0 = (np.log(delta) + 1.5 * np.log(tau) + para_hydrogen["a_k"][1] +
               para_hydrogen["a_k"][2] * tau + k_sum_3_9)

    # print(f"alpha_0: {alpha_0:.4f}")

    return alpha_0


# Berechnung von alpha_r
def alpha_r(delta, tau, para_hydrogen):

    i_sum_1_7 = 0
    for i in range(1, 8):
        i_sum_1_7 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * tau**para_hydrogen["t_i"][i])

    i_sum_8_9 = 0
    for i in range(8, 10):
        i_sum_8_9 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * tau**para_hydrogen["t_i"][i] * np.exp(
                -delta**para_hydrogen["p_i"][i]))

    i_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        i_sum_10_14 += (N_i * delta**d_i * tau**t_i * np.exp(phi_i * (
            delta - D_i)**2 + beta_i * (tau - gamma_i)**2))

    alpha_r = i_sum_1_7 + i_sum_8_9 + i_sum_10_14

    # print(f"alpha_r: {alpha_r:.4f}")

    return alpha_r


# alpha


def alpha(alpha_0, alpha_r):

    alpha = alpha_0 + alpha_r

    # print(f"alpha: {alpha:.4f}")

    return alpha

##############################################
# ########   Partielle Ableitungen ######### #
##############################################

######################################################################
# ###################### Erste Ableitungen ########################  #
######################################################################
# Für Druck alpha_r nach delta
# Für Enthalpie alpha_r nach delta; alpha_r nach tau; alpha_0 nach tau
######################################################################

# 1. alpha_0 nach tau
def alpha_0_tau_calc(delta, tau, para_hydrogen):
    exp_sum_3_9 = 0
    for k in range(3, 10):
        if para_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            a_k = para_hydrogen["a_k"].get(k, 0)
            b_k = para_hydrogen["b_k"].get(k, 0)
            exp_sum_3_9 += a_k * b_k * (-np.exp(b_k * tau) / (
                1 - np.exp(b_k * tau)))

    alpha_0_tau = 1.5/tau + para_hydrogen["a_k"][2] + exp_sum_3_9

    # print(f"Ableitung von alpha_0 nach tau: {alpha_0_tau:.3f}")

    return alpha_0_tau

# 2. alpha_r nach tau
def alpha_r_tau_calc(delta, tau, para_hydrogen):

    artau_sum_1_7 = 0
    for i in range(1, 8):
        artau_sum_1_7 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * para_hydrogen["t_i"][i] * tau**(para_hydrogen[
                "t_i"][i]-1))

    artau_sum_8_9 = 0
    for i in range(8, 10):
        artau_sum_8_9 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * para_hydrogen["t_i"][i] * tau**(para_hydrogen[
                "t_i"][i]-1) * np.exp(-delta**para_hydrogen["p_i"][i]))

    artau_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        artau_sum_10_14 += (N_i*delta**d_i * np.exp(phi_i*(
            delta - D_i)**2 + beta_i*(tau - gamma_i)**2) * (
                t_i * tau**(t_i - 1) + 2*beta_i*tau**t_i*(tau - gamma_i)))

    alpha_r_tau = artau_sum_1_7 + artau_sum_8_9 + artau_sum_10_14

    # print(f"Ableitung von alpha_r nach tau: {alpha_r_tau:.3f}")

    return alpha_r_tau


# 3. alpha_r nach delta
def alpha_r_delta_calc(delta, tau, para_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (para_hydrogen["N_i"][i] *
                            delta**(para_hydrogen["d_i"][i]-1) *
                            para_hydrogen["d_i"][i] *
                            tau**para_hydrogen["t_i"][i])

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        p_i = para_hydrogen["p_i"][i]
        t_i = para_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (N_i * tau**t_i * np.exp(-delta**p_i) *
                            (d_i*delta**(d_i-1) - p_i*delta**(d_i + p_i - 1)))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        ardelta_sum_10_14 += (N_i * tau**t_i * np.exp(
            phi_i*(delta-D_i)**2 + beta_i*(tau-gamma_i)**2) * (
                d_i*delta**(d_i-1) + 2*delta**d_i * phi_i*(delta-D_i)
            ))

    alpha_r_delta = ardelta_sum_1_7 + ardelta_sum_8_9 + ardelta_sum_10_14

    # print(f"Ableitung von alpha_r nach delta: {alpha_r_delta:.3f}")

    return alpha_r_delta

#######################################################################################
# ############################### Zweite Ableitungen ################################ #
#######################################################################################
# Für c_v alpha_0 nach tau nach tau und alpha_r nach tau nach tau
# Für c_p alpa_r nach delta; alpha_r nach delta nach tau; alpha_r nach delta nach delta
#######################################################################################

# 1. alpha_0_tau nach tau
def alpha_0_tau_tau_calc(delta, tau, para_hydrogen):
    exp_sum_3_9 = 0
    for k in range(3, 10):
        if para_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            a_k = para_hydrogen["a_k"].get(k, 0)
            b_k = para_hydrogen["b_k"].get(k, 0)
            exp_sum_3_9 += a_k * b_k**2 * (-np.exp(b_k * tau) / (1 - np.exp(
                b_k * tau))**2)

    alpha_0_tau_tau = -1.5/tau**2 + exp_sum_3_9

    # print(f"Ableitung von alpha_0 nach tau: {alpha_0_tau:.3f}")

    return alpha_0_tau_tau

# 2. alpha_r_tau nach tau
def alpha_r_tau_tau_calc(delta, tau, para_hydrogen):

    artau_sum_1_7 = 0
    for i in range(1, 8):
        artau_sum_1_7 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * para_hydrogen["t_i"][i] * tau**(para_hydrogen[
                "t_i"][i]-2) * (para_hydrogen["t_i"][i] - 1))

    artau_sum_8_9 = 0
    for i in range(8, 10):
        artau_sum_8_9 += (para_hydrogen["N_i"][i] * delta**para_hydrogen[
            "d_i"][i] * para_hydrogen["t_i"][i] * tau**(para_hydrogen[
                "t_i"][i]-2) * np.exp(-delta**para_hydrogen["p_i"][i]) * (
                    para_hydrogen["t_i"][i] - 1))

    artau_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        artau_sum_10_14 += (N_i*delta**d_i * np.exp(phi_i*(
            delta - D_i)**2 + beta_i*(tau - gamma_i)**2) * (
                2*beta_i*(tau - gamma_i) * (
                    t_i*tau**(t_i-1) + 2*beta_i*tau**t_i * (
                        tau - gamma_i)) + t_i*tau**(t_i - 2)*(
                            t_i - 1) + 2*beta_i*(
                                t_i*tau**(t_i - 1)*(
                                    tau-gamma_i) + tau**t_i)))

    alpha_r_tau_tau = artau_sum_1_7 + artau_sum_8_9 + artau_sum_10_14

    # print(f"Ableitung von alpha_r_tau nach tau: {alpha_r_tau_tau:.3f}")

    return alpha_r_tau_tau


# 3. alpha_r_delta nach tau
def alpha_r_delta_tau_calc(delta, tau, para_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (para_hydrogen["N_i"][i] *
                            delta**(para_hydrogen["d_i"][i] - 1) *
                            para_hydrogen["d_i"][i] *
                            tau**(para_hydrogen["t_i"][i] - 1) * (
                                para_hydrogen["t_i"][i]))

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        p_i = para_hydrogen["p_i"][i]
        t_i = para_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (N_i * t_i * tau**(t_i - 1) * np.exp(
            -delta**p_i) * (d_i*delta**(d_i-1) - p_i*delta**(
                d_i + p_i - 1)))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        ardelta_sum_10_14 += (N_i * np.exp(phi_i*(delta-D_i)**2 + beta_i*(
            tau-gamma_i)**2) * (d_i*delta**(d_i-1) + 2*delta**d_i * phi_i*(
                delta-D_i) * (t_i * tau**(t_i - 1) + 2*beta_i*(
                    tau - gamma_i) * tau**t_i)))

    alpha_r_delta_tau = ardelta_sum_1_7 + ardelta_sum_8_9 + ardelta_sum_10_14

    # print(f"Ableitung von alpha_r nach delta: {alpha_r_delta_tau:.3f}")

    return alpha_r_delta_tau

# 4. alpha_r_delta nach delta
def alpha_r_delta_delta_calc(delta, tau, para_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (para_hydrogen["N_i"][i] *
                            delta**(para_hydrogen["d_i"][i]-2) *
                            para_hydrogen["d_i"][i] *
                            tau**para_hydrogen["t_i"][i] * (
                                para_hydrogen["d_i"][i] - 1))

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        p_i = para_hydrogen["p_i"][i]
        t_i = para_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (
            N_i * tau**t_i * np.exp(-delta**p_i) *
            (-p_i * delta**(p_i - 1)*(d_i*delta**(
                d_i-1) - p_i*delta**(d_i + p_i - 1)) + (
                    d_i*(d_i - 1) * delta**(d_i - 2) - p_i*(
                        d_i + p_i - 1) * delta**(
                            d_i + p_i - 2))))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = para_hydrogen["N_i"][i]
        d_i = para_hydrogen["d_i"][i]
        t_i = para_hydrogen["t_i"][i]
        phi_i = para_hydrogen["phi_i"][i]
        D_i = para_hydrogen["D_i"][i]
        beta_i = para_hydrogen["beta_i"][i]
        gamma_i = para_hydrogen["gamma_i"][i]
        ardelta_sum_10_14 += (N_i * tau**t_i * np.exp(
            phi_i*(delta-D_i)**2 + beta_i*(tau-gamma_i)**2) * (
                d_i * (d_i - 1) * delta**(d_i - 2) + 2*d_i*delta**(
                    d_i - 1) * phi_i * (delta - D_i) + (
                        2*phi_i*delta**d_i) + 2*phi_i*(
                            delta - D_i) * (d_i * delta**(
                                d_i-1) + 2*delta**d_i*phi_i*(
                                    delta - D_i))))

    alpha_r_delta_delta = (
        ardelta_sum_1_7 + ardelta_sum_8_9 + ardelta_sum_10_14)

    # print(f"Ableitung von alpha_r nach delta: {alpha_r_delta:.3f}")

    return alpha_r_delta_delta

########################################################################

###########################################################
# ###### Helmholtzsche Modelle für Druck und Dichte ##### #
###########################################################

def initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen):
    # Vapour pressure nach Leachmann (2009)
    N_i = 0
    k_i = 0
    N_eta_sum = 0
    for i in range(1,5):
        N_i = para_hydrogen["N_i_s"][i]
        k_i = para_hydrogen["k_i"][i]
        N_eta_sum += N_i * (abs((1 - (T_ZP / T_c))) ** k_i)

    p_v = p_c * math.exp((T_c / T_ZP) * N_eta_sum)

    if p_ZP < p_v:
        # Zustandspunkt untehalb der Verdampfungslinie
        # Gasförmiger Zustand
        # Ideale Gasgleichung:
        rho_init = p_ZP/(R*T_ZP)  # kg/m^3 Startwert für die Dichte
    else:
        # Zustandsgleichung oberhalb der Verdampfungslinie
        # Flüssiger Zustand
        rho_init = 70.9  # kg/m^3 Startwert für die Dichte
    
    return rho_init


def berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                    max_iteration, differenz):

    for iteration in range(max_iteration):
        # Reduzierte Dichte und Temperatur
        delta = rho_i / rho_c
        tau = T_c / T_ZP

        # Berechnung des Drucks mit geratenem rho
        alpha_r_delta = alpha_r_delta_calc(delta, tau, para_hydrogen)
        p_raten = rho_i * R * T_ZP * (1 + delta * alpha_r_delta)

        # Zentrale finite Differenz für dp/drho
        ################################################
        delta_rho = rho_i / 1000
        rho_delta_plus = rho_i + delta_rho
        rho_delta_minus = rho_i - delta_rho

        p_delta_plus = rho_delta_plus * R * T_ZP * (
            1 + (rho_delta_plus / rho_c) * alpha_r_delta_calc(
                rho_delta_plus/rho_c, tau, para_hydrogen))

        p_delta_minus = rho_delta_minus * R * T_ZP * (
            1 + (rho_delta_minus / rho_c) * alpha_r_delta_calc(
                rho_delta_minus / rho_c, tau, para_hydrogen))

        dp_nach_dr = (p_delta_plus - p_delta_minus) / (2 * delta_rho)

        # Neue Schätzung von rho
        rho_i_plus = rho_i + (p_ZP - p_raten) / (dp_nach_dr)

        # Prüfen, ob Konvergenz erreicht ist
        if abs(rho_i_plus - rho_i) < differenz:
            # Bei Konvergenz zurückgeben
            # print(f"Konvergenz nach {iteration + 1} Iterationen erreicht.")
            return rho_i_plus

        # Falls keine Konvergenz, rho aktualisieren
        rho_i = rho_i_plus

    # Falls nach max_iteration keine Konvergenz erreicht wurde
    raise ValueError("Newton-Raphson-Verfahren konvergierte nicht.")


#############################################
# ########## Enthalpie h [J/kg] ########### #
#############################################

def calc_H2_enthalpy(T_ZP, p_ZP):
    
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    a_0_t = alpha_0_tau_calc(delta_neu, tau, para_hydrogen)
    a_r_t = alpha_r_tau_calc(delta_neu, tau, para_hydrogen)
    a_r_d = alpha_r_delta_calc(delta_neu, tau, para_hydrogen)

    h = R * T_ZP * (1 + tau * (a_0_t + a_r_t) + delta_neu * a_r_d) - 522e3

    return h

#############################################
# ########## Entropie s [J/kgK] ########### #
#############################################

def calc_H2_entropy(T_ZP, p_ZP):
    
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    a_0_t = alpha_0_tau_calc(delta_neu, tau, para_hydrogen)
    a_r_t = alpha_r_tau_calc(delta_neu, tau, para_hydrogen)
    a_0 = alpha_0(delta_neu, tau, para_hydrogen)
    a_r = alpha_r(delta_neu, tau, para_hydrogen)

    s = R * (tau * (a_0_t + a_r_t) - a_0 - a_r)

    return s

#############################################
# ########## Dichte rho [kg/m3] ########### #
#############################################

def calc_H2_density(T_ZP, p_ZP):
    
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                      max_iteration, differenz)

    return rho

###################################################
# ## c_v  (isochore Wärmekapazität) [J/(kg*K)] ## #
###################################################

def calc_H2_cv(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    a_0_t_t = alpha_0_tau_tau_calc(delta_neu, tau, para_hydrogen)
    a_r_t_t = alpha_r_tau_tau_calc(delta_neu, tau, para_hydrogen)

    c_v = - R * tau**2 * (a_0_t_t + a_r_t_t)

    return c_v

###################################################
# ## c_p  (isobare Wärmekapazität)  [J/(kg*K)] ## #
###################################################

def calc_H2_cp(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, para_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, para_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    c_v = calc_H2_cv(T_ZP, p_ZP)
    
    a_r_d_t = alpha_r_delta_tau_calc(delta_neu, tau, para_hydrogen)
    a_r_d_d = alpha_r_delta_delta_calc(delta_neu, tau, para_hydrogen)
    a_r_d = alpha_r_delta_calc(delta_neu, tau, para_hydrogen)

    Zaehler = (1 + delta_neu * a_r_d - delta_neu * tau * a_r_d_t)**2
    Nenner = 1 + 2 * delta_neu * a_r_d + delta_neu**2 * a_r_d_d
    c_p = c_v + R * Zaehler / Nenner

    return c_p

######################################################################

# Test run

# T_ZP = float(
#     input("Temperatur in Zustandspunkt T_ZP (K): "))
# p_ZP = float(
#     input("Druck in Zustandspunkt p_ZP (Pa): "))

# T_ZP = 250 #[K]
# p_ZP = 1e6 #[Pa]

# enthalpy_H2 = calculate_H2_enthalpy(T_ZP, p_ZP)
# cv_H2 = calculate_H2_cv(T_ZP, p_ZP)
# cp_H2 = calculate_H2_cp(T_ZP, p_ZP)
# eta_H2 = calculate_H2_eta(T_ZP, p_ZP)
# lambda_calc_H2 = calculate_H2_lambda(T_ZP, p_ZP)

# print(f"Enthalpie: {enthalpy_H2:.6f} J/kg")
# print(f"Isochore Wärmekapazität: {cv_H2:.6f} J/(kg*K)")
# print(f"Isobare Wärmekapazität: {cp_H2:.6f} J/(kg*K)")
# print(f"Viskosität: {eta_H2:.3f} [µPa*s]")
# print(f"Wärmeleitfähigkeit: {lambda_calc_H2:.6f} [W/(m*K)]")