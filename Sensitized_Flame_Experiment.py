# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:52:04 2020

@author: Kodo Bear
"""

import os
import shutil
import copy
import errno
import pickle
import time, sys
import numpy as np
import cantera as ct
import datetime

import flames
import common_functions as cf

ct.suppress_thermo_warnings() #Suppress cantera warnings!


def run_flame_simulation(mech, arrtype, pres, temp, fue, oxi, dilu, mix_params,
                         safi, par, Mingrid, Mul_soret, Loglevel):
    """
    Takes information from initializer and runs necessary functions to perform
    a one-dimensional simulation. Simulation results will be saved if booleans
    are set to True.

    Parameters
    ----------
    mech : str
        A .cti mechanism file containing all reaction and species information.
    arrtype : str
        Defines the scale that conditions are in. Either linear or logarithmic
    pres : list
        A list of pressure conditions to test over [initial, final, number of points].
    temp : list
        A list of temperature conditions to test over [initial, final, number of points].
    fue : str or list
        As a string the variable represents a single species of fuel being used.
        As a list the variable represents multicomponent fuel species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    oxi : str or list
        As a string the variable represents a single species of oxidizer being used.
        As a list the variable represents multicomponent oxidizer species
        followed by the percentage to the total oxidizer [Component1, % of total, ...]
    dilu : str or list
        As a string the variable represents a single species of diluent being used.
        As a list the variable represents multicomponent diluent species
        followed by the percentage to the total diluent [Component1, % of total, ...]
    mix_params : list
        A list of the two mixture parameters and mixtrue type used in creating
        a mixture.
    safi : boolean
        If true simulation conditions and ranking results will be saved.
    par : bool
        If true, run simulations in paralle
    Mingrid: int
        Number of points to be solved in the simulation
    Mul_soret : boolean
        Multicomponent diffuction and Soret effect are calculated if true
    Loglevel : int
        A number from 1 to 10. The larger the number the more information that
        is printed during the simulation to the user.

    Returns
    -------
    None.

    """
    mechan = cf.model_folder(mech)
    condi = initialization(mechan, arrtype, pres, temp, fue, oxi, dilu,
                           mix_params, Mingrid, Mul_soret, Loglevel)
    paralist = cf.case_maker(condi)
    simtime = run_simulations(condi, paralist, par)
    flame_info, flame_info_unfiltered, sim_info = load_filter_flame_info(condi)
    siminfo = [simtime, sim_info[0], simtime + sim_info[1]]
    if safi:
        file_saving(condi, flame_info, paralist, siminfo)



def initialization(mechanism, array_type, Press, Temperature, fuel, oxidizer,
                   diluent, mix_params, mingrid, mul_soret, loglevel):
    """


    Parameters
    ----------
    mechanism : str
        A .cti mechanism file containing all reaction and species information.
    array_type : str
        Defines the scale that conditions are in. Either linear or logarithmic
    Press : list
        A list of pressure conditions to test over [initial, final, number of points].
    Temperature : list
        A list of temperature conditions to test over [initial, final, number of points].
    fuel : str or list
        As a string the variable represents a single species of fuel being used.
        As a list the variable represents multicomponent fuel species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    oxidizer : str or list
        As a string the variable represents a single species of oxidizer being used.
        As a list the variable represents multicomponent oxidizer species
        followed by the percentage to the total oxidizer [Component1, % of total, ...]
    diluent : str or list
        As a string the variable represents a single species of diluent being used.
        As a list the variable represents multicomponent diluent species
        followed by the percentage to the total diluent [Component1, % of total, ...]
    mix_params : list
        A list of the two mixture parameters and mixtrue type used in creating
        a mixture.
    mingrid : int
        Number of points to be solved in the simulation
    mul_soret : boolean
        Multicomponent diffuction and Soret effect are calculated if true
    loglev : int
        A number from 1 to 10. The larger the number the more information that
        is printed during the simulation to the user.

    Returns
    -------
    conditions : dict
        Simulation information organized into a dictionary with the following
        structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture': [fuel, diluent, oxidizer],
             'Flame': [mingrid, mul_soret, loglevel],
             'Files': [mechanism, save_path]}

    """
    #Working directory
    parent_dir = 'Flame_Sensitivity_Results'
    try:
        os.makedirs(parent_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y_%m_%d %H.%M.%S Flame_Speed_Sens")
    directory = dt_string
    save_path = os.path.join(parent_dir, directory)
    os.makedirs(save_path)

    os.makedirs(os.path.join(save_path, 'all_flame_sims'))

    shutil.copyfile('input.yaml', os.path.join(save_path, 'input.yaml'))

    oxidizer = cf.normalize_mixture(oxidizer)
    fuel = cf.normalize_mixture(fuel)
    diluent = cf.normalize_mixture(diluent)

    conditions = {'Parameters': [Press, Temperature, mix_params, array_type],
                  'Mixture': [fuel, diluent, oxidizer],
                  'Flame': [mingrid, mul_soret, loglevel],
                  'Files': [mechanism, save_path]}
    return conditions


def run_simulations(conditions, paramlist, par):
    """
    Runs the simulation through all functions inlcuding the simulation
    calculations and the rankings.

    Parameters
    ----------
    conditions : dict
        Simulation information organized into a dictionary with the following
        structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture': [fuel, diluent, oxidizer],
             'Flame': [mingrid, mul_soret, loglevel],
             'Files': [mechanism, flame_temp]}
    paramlist : list
        Simulation case information the following structure:
        [[Pressure, Temperature, Mixture], ...]
    par : bool
        If true, run simulations in paralle

    Returns
    -------
    flame_info_filtered : list
        A filtered list where duplicate reactions in f_sens are summed together
    flame_info_unfiltered : list
        Unfiltered list of dictionary information iof the results of the
        simulation and conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    sim_time : float
        The time it took to run the simulations.

    """
    tic = time.time()
    print('Initial number of cases: '+format(len(paramlist)))
    print('\nStart of simulations...')
    sim_start  = time.time()
    cf.parallelize(paramlist, conditions, flame_sens, par)
    sim_end    = time.time()
    sim_time   = sim_end - sim_start
    print('End of simulations')
    print('Simulations took '+format(sim_time, '0.5f')+' seconds.')
    return sim_time

def load_filter_flame_info(conditions):
    """
    Load and filter the individual flame simulations.
    Deal with duplicate reactions.

    Parameters
    ----------
    conditions : dict
        Simulation information organized into a dictionary with the following
        structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture': [fuel, diluent, oxidizer],
             'Flame': [mingrid, mul_soret, loglevel],
             'Files': [mechanism, flame_temp]}

    Returns
    -------
    flame_info_filtered : list
        list of flame_info dictionaries after filtering
    flame_info_unfiltered : list
        list of unfiltered flame_info dictionaries
    sim_info : list
        [# of converged simulations, time for filtering]

    """
    print('Loading and filtering flame information...')
    flame_info = collect_flame_info(conditions['Files'][1])
    filter_start = time.time()
    converged    = 0
    for x in flame_info:
        if x['Flame'][0] is None:
            continue
        else:
            converged += 1
    flame_info_unfiltered = copy.deepcopy(flame_info)

    chem = conditions['Files'][0]
    gas = ct.Solution(chem)
    duplicate_rxns = cf.duplicate_reactions(gas)
    flame_info_filtered = flame_info_filter(flame_info, duplicate_rxns)
    filter_end  = time.time()
    filter_time = filter_end - filter_start
    print('Loading and Filtering took '+format(filter_time, '0.5f')+ ' seconds.')
    print('\nNumber of cases converged:' +str(converged))
    sim_info = [converged, filter_time]
    return flame_info_filtered, flame_info_unfiltered, sim_info


def flame_sens(p, T, mix, cond):
    """
    Runs a single flame simulation.

    Parameters
    ----------
    P : list
        A list of pressure conditions to test over [initial, final, number of points].
    T : list
        A list of temperature conditions to test over [initial, final, number of points].
    mix : dict
        Dictionary where keys are mixture components, values are mole fracs.
    cond : dict
        Simulation information organized into a dictionary with the following
        structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture': [fuel, diluent, oxidizer],
             'Flame': [mingrid, mul_soret, loglevel],
             'Files': [mechanism, flame_temp]}

    Returns
    -------
    flame_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}

    """
    at            = cond['Parameters'][3]
    chem          = cond['Files'][0]
    tempfile      = cond['Files'][1]
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    mg            = cond['Flame'][0]
    ms            = cond['Flame'][1]
    logl          = cond['Flame'][2]

    # More summary parameters to save later
    Fuel = 'unused'
    Oxidizer = 'unused'
    Fue_Percent = cf.mixture_percentage(Fuel_name, mix)
    Oxi_Percent = cf.mixture_percentage(Oxidizer_name, mix)
    Dil_Percent = cf.mixture_percentage(Diluent_name, mix)
    mixture = flames.Mixture(mix, chem)  # Create mixture object
    phi = mixture.phi

    f = flames.Flame(mixture, p, T, tempfile, chemfile=chem)
    f.run(mingrid=mg, loglevel=logl, mult_soret=ms)
    if f.flame_result is None:
        flame_info = {'Flame': [None, 'Flame did not converge'],
                      'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent,
                                     at]}
    else:
        f.sensitivity()
        f_sens = f.sens  # Rxn sensitivities
        Su = f.flame_result.velocity[0]  # Flame speed at the front
        flame_T = f.flame_result.T[-1]  # Flame temperature at end
        flame_rho = f.flame_result.density_mass[0]  # Flame density at the front
        flame_info = {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
                      'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent,
                                     at]}
    now = datetime.datetime.now()
    filename = now.strftime("simulation %H_%M_%S.pkl")
    with open(os.path.join(tempfile, 'all_flame_sims', filename), 'wb') as f:
        pickle.dump(flame_info, f)
    return None


def flame_info_filter(flame_information, duplicate_reactions):
    """
    Takes the sensitivty of duplicate reactions and summs them togther to get
    a total sensitivity of the reaction.

    Parameters
    ----------
    flame_information : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    duplicate_reactions : dict
        Dictionary containing duplicate reactions of the following format:
        dup_rxns = {'Reaction Equation 1': [Rxn Number_a, Rxn Number_b]
                            :
                            :
                    'Reaction Equation N': [Rxn Number_a, Rxn Number_b]}

    Returns
    -------
    flame_information : dict
        A filtered list where duplicate reactions in f_sens are summed together.

    """
    for f in flame_information:
        if f['Flame'][0] is None:
            continue
        else:
            for d in duplicate_reactions:
                sum = 0
                for n in duplicate_reactions[d]:
                     sum += f['Flame'][0][n][1]
                     f['Flame'][0][n][1] = 0
                f['Flame'][0][duplicate_reactions[d][0]][1] = sum
    return flame_information


def file_saving(cond, fla_inf, p_list, s_info):
    """
    Information from the simulation is pickled and saved in a generated folder
    with a text document filled with details of the simulation performed.

    Parameters
    ----------
    cond : dict
        Simulation information organized into a dictionary with the following
        structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture': [fuel, diluent, oxidizer],
             'Flame': [mingrid, mul_soret, loglevel],
             'Files': [mechanism, flame_temp]}
    fla_inf : list
        Information of the results of the simulation as well as simulation
        conditions. List of dictionaries with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    p_list : list
        Simulation case information the following structure:
        [[Pressure, Temperature, Mixture], ...]
    s_info : list
        A list of information of from performing the simulation including
        the time it took to run the simulations, the number of cases that
        converged, and the duration of the function.

    Returns
    -------
    None.

    """
    save_time_start = time.time()
    p             = cond['Parameters'][0]
    T             = cond['Parameters'][1]
    mix_params    = cond['Parameters'][2]
    chem          = cond['Files'][0]
    save_path = cond['Files'][1]
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    mg            = cond['Flame'][0]
    ms            = cond['Flame'][1]
    s_time        = datetime.timedelta(seconds=s_info[0])
    conv          = s_info[1]
    dur           = datetime.timedelta(seconds=s_info[2])

    figures_dir = 'Flame_Sensitivity_Plots'
    figure_path = os.path.join(save_path, figures_dir)
    os.makedirs(figure_path)
    print('Directory Created')

    print('\nCreating text file...')
    #Text Description
    filename  = 'Case Description.txt'
    filename  = os.path.join(save_path, filename)
    f         = open(filename, "w")

    text_description = ("This file provides simulation information.\n"
                        "The following information are the parameters "
                        "and cases simulated\n\n" +
                        cf.parameters_string('1D', p, T, mix_params, chem,
                                               Fuel_name, Oxidizer_name,
                                               Diluent_name) +
                        "\n\n======Flame Simulation Information======"
                        "\nMingrid = " + format(mg) +
                        "\nMult_Soret = " + format(ms) +
                        "\n========================================\n"
                        "\n=============Time/Converged=============\n"
                        "Total Cases: " + format(len(p_list)) +
                        "\nSim time: " + str(s_time) +
                        "\nCases Converged: " + str(conv) +
                        "\nRun time: " + str(dur) +
                        "\n========================================")

    f.write(text_description)
    f.close()
    print('Text file created')

    print('\nStart file saving...')
    file = os.path.join(save_path, 'Flame Information.pkl')
    with open(file, 'wb') as f:
        pickle.dump(fla_inf, f)
    shutil.rmtree(os.path.join(save_path, 'all_flame_sims')) # clean up individual flame files
    save_time_end = time.time()
    save_time = save_time_end - save_time_start
    print('Total File Save Time: '+format(save_time, '0.5f')+' seconds.\n')
    print('End file saving')


if __name__ == "__main__":
    #Set experiment parameters
    mech_name = 'Li_model_modified_trioxane.cti' #Mechanism file
    Mechanism = cf.model_folder(mech_name)

    #Parameters for main loop
    # P and Phi are used in all cases.
    # Whether FTD or OtD is used depends on mixture type.
    # Each list should follow patter [First Point, Last Point, # of Points]
    # Array_type's (log, lin)
    #  log creates a logspace array of parameters
    #  lin creates a linspace array of parameters
    Array_type = 'lin'
    Press      = [0.5, 1, 2]  #Pressure [atm]
    E_Ratio    = [0.25, 1.2, 2] #Equivalence ratio
    FO_to_D    = [0.05, 0.95, 2] #Amount of Fuel/Oxidizer to Diluent........This isn't used?
    F_to_D     = [0.75, 0.95, 2] #Fuel/(Fuel + Diluent)
    O_to_D     = [0.05, 0.95, 2] #Oxidizer/(Oxidizer + Diluent)

    #Initial temperature of unburned mixture
    Tint = 373 #Temperature [K]

    #Parameters for mixture (Fuel, Oxidizer, Diluent)
    Fuel     = 'C3H6O3' #chemical formula of fuel
    # Fuel  = ['CH4', .50 , 'CH3OH', .50]
    Oxidizer = 'O2' #chemical formula of oxidizer
    # Oxidizer = ['O2', .35 , 'NO2', .65]
    Diluent  = 'N2' #chemical formula of diluent
    Air      = 'O2:1, N2:3.76' #chemical components for air as an oxidizer

    #Flame Conditions
    Mingrid   = 200
    Mul_soret = False
    Loglevel  = 0

    #True/False statements
    Save_files = True # If true, save files for plotting script

    #Mixture Type (Debug, Custom, Oxi_Dil, Fue_Dil)
    #Provide one of the four types of mixtures into
    # variable mixture_type as a string
    #  Debug is used for simulating a single flame to check code
    #  Custom is under construction
    #  Oxi_Dil creates a mixture where the Diluent is a ratio of the Oxidizer used
    #  Fue_Dil creates a mixture where the Diluent is a ratio of the Fuel used
    Mixture_type = 'Oxi_Dil'

    #Run Code
    run_flame_simulation(Mechanism, Array_type, Press, E_Ratio, F_to_D,
                         O_to_D, Tint, Fuel, Oxidizer, Diluent, Air, Mingrid,
                         Mul_soret, Loglevel, Mixture_type, Save_files)

def collect_flame_info(path):
    """
    Collect flame information from individual pickled files or summary file.

    Parameters
    ----------
    path : string
        Path to the main folder for a group of simulations, usually formatted
        as Y_M_D H.M.S Flame_Speed_Sens

    Returns
    -------
    flame_info : list
        list of flame information dictionaries.

    """
    all_sim_path = os.path.join(path, 'all_flame_sims')
    if os.path.exists(all_sim_path):  # Load each individual simulation
        flame_info = []
        for file in os.listdir(all_sim_path):
            with open(os.path.join(all_sim_path, file), 'rb') as f:
                flame_info.append(pickle.load(f))
    else:
        with open(os.path.join(path, 'Flame Information.pkl'), 'rb') as f:
            flame_info = pickle.load(f)

    return flame_info