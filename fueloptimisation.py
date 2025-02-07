import matplotlib.pyplot as plt
import fuelflow
import pygad

def station_update(mass_flow, temperature, pressure, enthalpy, jetaflow):
    mass_flow.append(jetaflow.qm)
    temperature.append(jetaflow.t)
    pressure.append(jetaflow.p)
    enthalpy.append(jetaflow.calc_h())
    return mass_flow, temperature, pressure, enthalpy

def reference(qm_cb, qm_r, qm_t, p_lpfp):
    stations = ["LPP Eintritt", "LPP Austritt", "Rezirkulation Austritt", "FOHE Austritt", "HPP Austritt", "Brennkammer Eintritt", "FOHE-IDG Austritt", "FRTT", "Rezirkulation"]
    mass_flow = list()
    temperature = list()
    pressure = list()
    enthalpy = list()
    t0 = 250
    p0 = 0.4e5
    eta_lpfp = 0.83
    Q_fohe = 200000
    p_hpfp = 3e6
    eta_hpfp = 0.88
    Q_idg = 5500
    dp_inj = 2e5
    # init flow LPP Eintritt
    jetaflow = fuelflow.JetaFlow(qm_cb+qm_t, t0, p0)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate LPP Austritt
    P_lpp, _ = jetaflow.pump_hydraulic(p_lpfp, eta_lpfp)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate Rezirkulation Austritt
    recirc = fuelflow.JetaFlow(qm_r, 300, jetaflow.p)
    recirc.raise_to_h(h_r)
    _ = jetaflow.mix_flows(recirc)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate FOHE Austritt
    _ = jetaflow.heat_fixed_power(Q_fohe)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate HPP Austritt
    P_hpp, _ = jetaflow.pump_hydraulic(p_hpfp, eta_hpfp)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate Brennkammer Eintritt
    cb_ff = jetaflow.split_flows(qm_cb)
    cb_ff.reduce_pressure(cb_ff.p-dp_inj)
    jetaflow.reduce_pressure(jetaflow.p - 2e6)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, cb_ff)
    # calculate FOHE-IDG Austritt
    _ = jetaflow.heat_fixed_power(Q_idg)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    # calculate FRTT
    frtt = jetaflow.split_flows(qm_t)
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, frtt)
    # calculate Rezirkulation
    h_ra = jetaflow.calc_h()
    mass_flow, temperature, pressure, enthalpy = station_update(mass_flow, temperature, pressure, enthalpy, jetaflow)
    return stations, mass_flow, temperature, pressure, enthalpy, h_ra,  P_lpp, P_hpp


def fitness_func(ga_instance, solution, solution_idx):
    qm_cb, qm_r, qm_t, p_lpfp, h_r = solution
    try:
        _, mass_flow, temperature, pressure, enthalpy, h_ra,  P_lpp, P_hpp = reference(qm_cb, qm_r, qm_t, p_lpfp, h_r)
        fitness = -(P_lpp + P_hpp)
    except:
        fitness = -9999999
    return fitness

num_generations = 50
num_parents_mating = 4

sol_per_pop = 8
num_genes = 5

init_range_low = -2
init_range_high = 5

parent_selection_type = "sss"
keep_parents = 1

crossover_type = "single_point"

mutation_type = "random"
mutation_percent_genes = 60

def on_gen(ga_instance):
    print("Generation : ", ga_instance.generations_completed)
    print("Fitness of the best solution :", ga_instance.best_solution()[1])

ga_instance = pygad.GA(
    num_generations=num_generations,
    num_parents_mating=num_parents_mating,
    fitness_func=fitness_func,
    sol_per_pop=sol_per_pop,
    num_genes=num_genes,
    init_range_low=init_range_low,
    init_range_high=init_range_high,
    on_generation=on_gen,
    parent_selection_type=parent_selection_type,
    keep_parents=keep_parents,
    crossover_type=crossover_type,
    mutation_type=mutation_type,
    mutation_percent_genes=mutation_percent_genes)

ga_instance.run()
ga_instance.plot_fitness()

temp = fuelflow.JetaFlow(1, 400, 1e5)
stations, mass_flow, temperature, pressure, enthalpy, h_ra,  P_lpp, P_hpp = reference(0.3, 0.4, 0.3, 5e5, temp.calc_h())

for i in range(len(stations)):
    print(stations[i], str(mass_flow[i]), str(temperature[i]), str(pressure[i]), str(enthalpy[i]))
