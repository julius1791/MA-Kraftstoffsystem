import numpy as np
import math

# A - Helmholtzsche freie Energie A = f(rho, T) [J/kg]
# alpha - reduzierte Energie alpha = A/RT = f(delta, tau)
# delta: reduzierte Dichte delta=rho/rho_c
# tau: reduzierte reziproke Temperatur tau=T_c/T

# alpha = alpha_0(tau, deltha) + alpha_r(tau, deltha)

##################################################

# Kritische Werte
rho_c = 31.262     # [kg/m^3] Kritische Dichte für H2
T_c   = 33.145     # [K]        ; Kritische Temperatur für H2
p_c   = 1.2964e6   # [Pa]       ; Kritischer Druck für H2

R     = 4124.2     # [J/(kg*K)] ; Spezifische Gaskonstante für H2

# Koeffizienten für Normalwasserstoff
normal_hydrogen = {
    "a_k": {
        1: -1.4579856475,
        2: 1.888076782,
        3: 1.616,
        4: -0.4117,
        5: -0.792,
        6: 0.758,
        7: 1.217,
        8: 0.0,  # Kein Wert gegeben
        9: 0.0,  # Kein Wert gegeben
    },
    "b_k": {
        1: 0.0,  # Kein Wert gegeben
        2: 0.0,  # Kein Wert gegeben
        3: -16.0205159149,
        4: -22.6580178006,
        5: -60.0090511389,
        6: -74.9434303817,
        7: -206.9392065168,
        8: 0.0,  # Kein Wert gegeben
        9: 0.0  # Kein Wert gegeben
    },
    "N_i": {
        1: -6.93643,
        2: 0.01,
        3: 2.1101,
        4: 4.52059,
        5: 0.732564,
        6: -1.34086,
        7: 0.130985,
        8: -0.777414,
        9: 0.351944,
        10: -0.0211716,
        11: 0.0226312,
        12: 0.032187,
        13: -0.0231752,
        14: 0.0557346
    },
    "t_i": {
        1: 0.6844,
        2: 1,
        3: 0.989,
        4: 0.489,
        5: 0.803,
        6: 1.1444,
        7: 1.409,
        8: 1.754,
        9: 1.311,
        10: 4.187,
        11: 5.646,
        12: 0.791,
        13: 7.249,
        14: 2.986
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
        10: -1.685,
        11: -0.489,
        12: -0.103,
        13: -2.506,
        14: -1.607
    },
    "beta_i": {
        10: -0.171,
        11: -0.2245,
        12: -0.1304,
        13: -0.2785,
        14: -0.3967
    },
    "gamma_i": {
        10: 0.7164,
        11: 1.3444,
        12: 1.4517,
        13: 0.7204,
        14: 1.5445
    },
    "D_i": {
        10: 1.506,
        11: 0.156,
        12: 1.736,
        13: 0.67,
        14: 1.662
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
    "A_1_i": {
        0: -3.40976e-1,
        1: 4.58820,
        2: -1.45080,
        3: 3.26394e-1,
        4: 3.16939e-3,
        5: 1.90592e-4,
        6: -1.13900e-6
    },
    "A_2_i": {
        0: 1.38497e2,
        1: -2.21878e1,
        2: 4.57151,
        3: 1.0,
        4: 0.0,  # Kein Wert gegeben
        5: 0.0,  # Kein Wert gegeben
        6: 0.0   # Kein Wert gegeben
    },
    "B_1_i": {
        1: 3.63081e-2,
        2: -2.07629e-2,
        3: 3.14810e-2,
        4: -1.43097e-2,
        5: 1.74980e-3
    },
    "B_2_i": {
        1: 1.83370e-3,
        2: -8.86716e-3,
        3: 1.58260e-2,
        4: -1.06283e-2,
        5: 2.80673e-3
    },
    "a_i": {
        0: 2.09630e-1,
        1: -4.55274e-1,
        2: 1.43602e-1,
        3: -3.35325e-2,
        4: 2.76981e-3
    },
    "b_i": {
        0: -0.1870,
        1: 2.4871,
        2: 3.7151,
        3: -11.0972,
        4: 9.0965,
        5: -3.8292,
        6: 0.5166
    },
    "N_i_s": {
        1: -4.89789,
        2: 0.988558,
        3: 0.349689,
        4: 0.499356
    },
    "k_i" : {
        1: 1,
        2: 1.5,
        3: 2,
        4: 2.85
    }
}

# Für Wärmeleitfähigkeit
C_1 = 6.24e-4
C_2 = -2.58e-7
C_3 = 0.837

# Für Viskosität
molar_mass   = 2.01588       # [g/mol]
sigma        = 0.297         # [nm]
epsilon_k_B  = 30.41         # [K]
rho_sc       = 90.909090909  # [kg/m^3]
# N_A          = 6.022137e23 # Avogrado-Zahl
c_1_v        = 6.43449673
c_2_v        = 4.56334068e-2
c_3_v        = 2.32797868e-1
c_4_v        = 9.58326120e-1
c_5_v        = 1.27941189e-1
c_6_v        = 3.63576595e-1

############################################

# Berechnung von alpha_0
def alpha_0(delta, tau, normal_hydrogen):

    k_sum_3_9 = 0
    for k in range(3, 10):
        if normal_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            k_sum_3_9 += normal_hydrogen["a_k"].get(k, 0) * np.log(1 - np.exp(
                normal_hydrogen["b_k"][k] * tau))

    alpha_0 = (np.log(delta) + 1.5 * np.log(tau) + normal_hydrogen["a_k"][1] +
               normal_hydrogen["a_k"][2] * tau + k_sum_3_9)

    # print(f"alpha_0: {alpha_0:.4f}")

    return alpha_0


# Berechnung von alpha_r
def alpha_r(delta, tau, normal_hydrogen):

    i_sum_1_7 = 0
    for i in range(1, 8):
        i_sum_1_7 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * tau**normal_hydrogen["t_i"][i])

    i_sum_8_9 = 0
    for i in range(8, 10):
        i_sum_8_9 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * tau**normal_hydrogen["t_i"][i] * np.exp(
                -delta**normal_hydrogen["p_i"][i]))

    i_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
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
def alpha_0_tau_calc(delta, tau, normal_hydrogen):
    exp_sum_3_9 = 0
    for k in range(3, 10):
        if normal_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            a_k = normal_hydrogen["a_k"].get(k, 0)
            b_k = normal_hydrogen["b_k"].get(k, 0)
            exp_sum_3_9 += a_k * b_k * (-np.exp(b_k * tau) / (
                1 - np.exp(b_k * tau)))

    alpha_0_tau = 1.5/tau + normal_hydrogen["a_k"][2] + exp_sum_3_9

    # print(f"Ableitung von alpha_0 nach tau: {alpha_0_tau:.3f}")

    return alpha_0_tau

# 2. alpha_r nach tau
def alpha_r_tau_calc(delta, tau, normal_hydrogen):

    artau_sum_1_7 = 0
    for i in range(1, 8):
        artau_sum_1_7 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * normal_hydrogen["t_i"][i] * tau**(normal_hydrogen[
                "t_i"][i]-1))

    artau_sum_8_9 = 0
    for i in range(8, 10):
        artau_sum_8_9 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * normal_hydrogen["t_i"][i] * tau**(normal_hydrogen[
                "t_i"][i]-1) * np.exp(-delta**normal_hydrogen["p_i"][i]))

    artau_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
        artau_sum_10_14 += (N_i*delta**d_i * np.exp(phi_i*(
            delta - D_i)**2 + beta_i*(tau - gamma_i)**2) * (
                t_i * tau**(t_i - 1) + 2*beta_i*tau**t_i*(tau - gamma_i)))

    alpha_r_tau = artau_sum_1_7 + artau_sum_8_9 + artau_sum_10_14

    # print(f"Ableitung von alpha_r nach tau: {alpha_r_tau:.3f}")

    return alpha_r_tau


# 3. alpha_r nach delta
def alpha_r_delta_calc(delta, tau, normal_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (normal_hydrogen["N_i"][i] *
                            delta**(normal_hydrogen["d_i"][i]-1) *
                            normal_hydrogen["d_i"][i] *
                            tau**normal_hydrogen["t_i"][i])

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        p_i = normal_hydrogen["p_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (N_i * tau**t_i * np.exp(-delta**p_i) *
                            (d_i*delta**(d_i-1) - p_i*delta**(d_i + p_i - 1)))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
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
def alpha_0_tau_tau_calc(delta, tau, normal_hydrogen):
    exp_sum_3_9 = 0
    for k in range(3, 10):
        if normal_hydrogen["b_k"].get(k, 0) != 0:  # b_k vorhanden sein muss
            a_k = normal_hydrogen["a_k"].get(k, 0)
            b_k = normal_hydrogen["b_k"].get(k, 0)
            exp_sum_3_9 += a_k * b_k**2 * (-np.exp(b_k * tau) / (1 - np.exp(
                b_k * tau))**2)

    alpha_0_tau_tau = -1.5/tau**2 + exp_sum_3_9

    # print(f"Ableitung von alpha_0 nach tau: {alpha_0_tau:.3f}")

    return alpha_0_tau_tau

# 2. alpha_r_tau nach tau
def alpha_r_tau_tau_calc(delta, tau, normal_hydrogen):

    artau_sum_1_7 = 0
    for i in range(1, 8):
        artau_sum_1_7 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * normal_hydrogen["t_i"][i] * tau**(normal_hydrogen[
                "t_i"][i]-2) * (normal_hydrogen["t_i"][i] - 1))

    artau_sum_8_9 = 0
    for i in range(8, 10):
        artau_sum_8_9 += (normal_hydrogen["N_i"][i] * delta**normal_hydrogen[
            "d_i"][i] * normal_hydrogen["t_i"][i] * tau**(normal_hydrogen[
                "t_i"][i]-2) * np.exp(-delta**normal_hydrogen["p_i"][i]) * (
                    normal_hydrogen["t_i"][i] - 1))

    artau_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
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
def alpha_r_delta_tau_calc(delta, tau, normal_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (normal_hydrogen["N_i"][i] *
                            delta**(normal_hydrogen["d_i"][i] - 1) *
                            normal_hydrogen["d_i"][i] *
                            tau**(normal_hydrogen["t_i"][i] - 1) * (
                                normal_hydrogen["t_i"][i]))

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        p_i = normal_hydrogen["p_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (N_i * t_i * tau**(t_i - 1) * np.exp(
            -delta**p_i) * (d_i*delta**(d_i-1) - p_i*delta**(
                d_i + p_i - 1)))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
        ardelta_sum_10_14 += (N_i * np.exp(phi_i*(delta-D_i)**2 + beta_i*(
            tau-gamma_i)**2) * (d_i*delta**(d_i-1) + 2*delta**d_i * phi_i*(
                delta-D_i) * (t_i * tau**(t_i - 1) + 2*beta_i*(
                    tau - gamma_i) * tau**t_i)))

    alpha_r_delta_tau = ardelta_sum_1_7 + ardelta_sum_8_9 + ardelta_sum_10_14

    # print(f"Ableitung von alpha_r nach delta: {alpha_r_delta_tau:.3f}")

    return alpha_r_delta_tau

# 4. alpha_r_delta nach delta
def alpha_r_delta_delta_calc(delta, tau, normal_hydrogen):

    ardelta_sum_1_7 = 0
    for i in range(1, 8):
        ardelta_sum_1_7 += (normal_hydrogen["N_i"][i] *
                            delta**(normal_hydrogen["d_i"][i]-2) *
                            normal_hydrogen["d_i"][i] *
                            tau**normal_hydrogen["t_i"][i] * (
                                normal_hydrogen["d_i"][i] - 1))

    ardelta_sum_8_9 = 0
    for i in range(8, 10):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        p_i = normal_hydrogen["p_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        ardelta_sum_8_9 += (
            N_i * tau**t_i * np.exp(-delta**p_i) *
            (-p_i * delta**(p_i - 1)*(d_i*delta**(
                d_i-1) - p_i*delta**(d_i + p_i - 1)) + (
                    d_i*(d_i - 1) * delta**(d_i - 2) - p_i*(
                        d_i + p_i - 1) * delta**(
                            d_i + p_i - 2))))

    ardelta_sum_10_14 = 0
    for i in range(10, 15):
        N_i = normal_hydrogen["N_i"][i]
        d_i = normal_hydrogen["d_i"][i]
        t_i = normal_hydrogen["t_i"][i]
        phi_i = normal_hydrogen["phi_i"][i]
        D_i = normal_hydrogen["D_i"][i]
        beta_i = normal_hydrogen["beta_i"][i]
        gamma_i = normal_hydrogen["gamma_i"][i]
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

def initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen):
    # Vapour pressure nach Leachmann (2009)
    N_i = 0
    k_i = 0
    N_eta_sum = 0
    for i in range(1,5):
        N_i = normal_hydrogen["N_i_s"][i]
        k_i = normal_hydrogen["k_i"][i]
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


def berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                    max_iteration, differenz):

    for iteration in range(max_iteration):
        # Reduzierte Dichte und Temperatur
        delta = rho_i / rho_c
        tau = T_c / T_ZP

        # Berechnung des Drucks mit geratenem rho
        alpha_r_delta = alpha_r_delta_calc(delta, tau, normal_hydrogen)
        p_raten = rho_i * R * T_ZP * (1 + delta * alpha_r_delta)

        # Zentrale finite Differenz für dp/drho
        ################################################
        delta_rho = rho_i / 1000
        rho_delta_plus = rho_i + delta_rho
        rho_delta_minus = rho_i - delta_rho

        p_delta_plus = rho_delta_plus * R * T_ZP * (
            1 + (rho_delta_plus / rho_c) * alpha_r_delta_calc(
                rho_delta_plus/rho_c, tau, normal_hydrogen))

        p_delta_minus = rho_delta_minus * R * T_ZP * (
            1 + (rho_delta_minus / rho_c) * alpha_r_delta_calc(
                rho_delta_minus / rho_c, tau, normal_hydrogen))

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

def calculate_H2_enthalpy(T_ZP, p_ZP):
    
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    a_0_t = alpha_0_tau_calc(delta_neu, tau, normal_hydrogen)
    a_r_t = alpha_r_tau_calc(delta_neu, tau, normal_hydrogen)
    a_r_d = alpha_r_delta_calc(delta_neu, tau, normal_hydrogen)

    h = R * T_ZP * (1 + tau * (a_0_t + a_r_t) + delta_neu * a_r_d)

    return h

###################################################
# ## c_v  (isochore Wärmekapazität) [J/(kg*K)] ## #
# #################################################

def calculate_H2_cv(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    a_0_t_t = alpha_0_tau_tau_calc(delta_neu, tau, normal_hydrogen)
    a_r_t_t = alpha_r_tau_tau_calc(delta_neu, tau, normal_hydrogen)

    c_v = - R * tau**2 * (a_0_t_t + a_r_t_t)

    return c_v

###################################################
# ## c_p  (isobare Wärmekapazität)  [J/(kg*K)] ## #
###################################################

def calculate_H2_cp(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                      max_iteration, differenz)
    
    tau = T_c / T_ZP
    delta_neu = rho / rho_c

    c_v = calculate_H2_cv(T_ZP, p_ZP)
    
    a_r_d_t = alpha_r_delta_tau_calc(delta_neu, tau, normal_hydrogen)
    a_r_d_d = alpha_r_delta_delta_calc(delta_neu, tau, normal_hydrogen)
    a_r_d = alpha_r_delta_calc(delta_neu, tau, normal_hydrogen)

    Zaehler = (1 + delta_neu * a_r_d - delta_neu * tau * a_r_d_t)**2
    Nenner = 1 + 2 * delta_neu * a_r_d + delta_neu**2 * a_r_d_d
    c_p = c_v + R * Zaehler / Nenner

    return c_p

############################################
# #        Viskosität   [µPa*s]         # #
############################################

def calculate_H2_eta(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                      max_iteration, differenz)

    # Korrelation nach Muzny et al. (2013), von NIST genutzt werden und für alle Aggregatzustände gültig
    # Zero-density limit of the viscosity
    ln_S_stern = 0
    S_stern = 0
    eta_0 = 0
    for i in range(0,5):
        a_i = normal_hydrogen["a_i"][i]
        ln_S_stern += a_i * ((math.log((T_ZP / epsilon_k_B))) ** i)

    S_stern = math.exp(ln_S_stern)

    eta_0 = (0.021357 * ((molar_mass * T_ZP) ** 0.5))/(S_stern * (sigma ** 2))

    # Initial-density coefficient of the viscosity
    b_sum_0_6 = 0
    B_stern = 0
    B_eta = 0
    eta_1 = 0
    for i in range(0,7):
        b_i = normal_hydrogen["b_i"][i]
        b_sum_0_6 += b_i * ((T_ZP / epsilon_k_B) ** (-i))

    B_stern = b_sum_0_6
    B_eta = B_stern * (sigma ** 3)

    eta_1 = B_eta * eta_0

    # Final functional form of the viscosity
    rho_r = rho / rho_sc
    T_r = T_ZP / T_c
    eta = eta_0 + eta_1 * rho + c_1_v * (rho_r ** 2) * math.exp(c_2_v * T_r + c_3_v / T_r + (c_4_v * (rho_r ** 2)) / (c_5_v + T_r) + c_6_v * (rho_r ** 6))

    return eta

#####################################################
# ###      Wärmeleitfähigkeit λ [W/(m*K)]       ### #
#####################################################

def calculate_H2_lambda(T_ZP, p_ZP):
    rho_i = initialize_rho(T_ZP, p_ZP, p_c, T_c, normal_hydrogen)

    differenz = 1e-9  # Konvergenzkriterium: (rho_i+1 - rho_i) < 10^-9
    max_iteration = 50
    # Max_Iterationen: Konvergenzkriterium
    # (Begrenzung falls keine Konvergenz erreicht wird)

    rho = berechne_dichte(rho_i, T_ZP, p_ZP, rho_c, T_c, R, normal_hydrogen,
                      max_iteration, differenz)
    
    # Korrelation nach Assael et al. (2011), von NIST genutzt und für alle Aggregatzustände gültig
    # Contribution to the thermal conductivity in the dilute-gas limit
    A_1_sum_0_6 = 0
    A_2_sum_0_6 = 0
    lambda_0 = 0
    for i in range(0, 7):
        A_1 = normal_hydrogen["A_1_i"][i]
        A_1_sum_0_6 += A_1 * (T_ZP / T_c) ** i
        if normal_hydrogen["A_2_i"].get(i, 0) != 0:  # A_2_i vorhanden sein muss
            A_2 = normal_hydrogen["A_2_i"].get(i, 0)
            A_2_sum_0_6 += A_2 * (T_ZP / T_c) ** i

    lambda_0 = A_1_sum_0_6 / A_2_sum_0_6

    # Excess thermal conductivity (contribution of all other effects at elevated densities including many-body collisions, molecular-velocity correlations, collisional transfer)
    B_sum_1_5 = 0
    delta_lambda = 0
    for i in range(1,6):
        B_1 = normal_hydrogen["B_1_i"][i]
        B_2 = normal_hydrogen["B_2_i"][i]
        B_sum_1_5 += (B_1 + B_2 * (T_ZP/T_c)) * (rho/rho_c) ** i

    delta_lambda = B_sum_1_5

    # Empirical critical enhancement (for state points relatively distant (about 10 - 15 K) from the critical point)
    delta_lambda_c = C_1 / (C_2 + abs((T_ZP / T_c) - 1)) * math.exp(-1 * (C_3 * ((rho / rho_c) - 1)) ** 2)

    # Final functional form of the thermal conductivity as the sum of previously calculated components
    lambda_H2 = lambda_0 + delta_lambda + delta_lambda_c

    return lambda_H2

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