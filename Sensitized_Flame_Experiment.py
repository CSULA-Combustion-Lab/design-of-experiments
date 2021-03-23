# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:52:04 2020

@author: Kodo Bear
"""

import os
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
                         mgrid, msoret, loglev, safi):
    """


    Parameters
    ----------
    mech : TYPE
        DESCRIPTION.
    arrtype : TYPE
        DESCRIPTION.
    pres : TYPE
        DESCRIPTION.
    temp : TYPE
        DESCRIPTION.
    fue : TYPE
        DESCRIPTION.
    oxi : TYPE
        DESCRIPTION.
    dilu : TYPE
        DESCRIPTION.
    mix_params : TYPE
        DESCRIPTION.
    mgrid : TYPE
        DESCRIPTION.
    msoret : TYPE
        DESCRIPTION.
    loglev : TYPE
        DESCRIPTION.
    safi : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    mechan = cf.model_folder(mech)
    condi = initialization(mechan, arrtype, pres, temp, fue, oxi, dilu,
                           mix_params, mgrid, msoret, loglev)
    paralist = cf.case_maker(condi)
    flame_info, flame_info_unfiltered, siminfo = run_simulations(condi,
                                                                  paralist)
    if safi:
        file_saving(condi, flame_info, paralist, siminfo)



def initialization(mechanism, array_type, Press, Temperature, fuel, oxidizer,
                   diluent, mix_params, mingrid, mul_soret, loglevel):
    """


    Parameters
    ----------
    mechanism : TYPE
        DESCRIPTION.
    array_type : TYPE
        DESCRIPTION.
    Press : TYPE
        DESCRIPTION.
    Temperature : TYPE
        DESCRIPTION.
    fuel : TYPE
        DESCRIPTION.
    oxidizer : TYPE
        DESCRIPTION.
    diluent : TYPE
        DESCRIPTION.
    mix_params : TYPE
        DESCRIPTION.
    mingrid : TYPE
        DESCRIPTION.
    mul_soret : TYPE
        DESCRIPTION.
    loglevel : TYPE
        DESCRIPTION.

    Returns
    -------
    conditions : TYPE
        DESCRIPTION.

    """
    #Working directory
    flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

    if type(oxidizer) is list:
        multioxidizer = True
    elif type(oxidizer) is str:
        multioxidizer = False

    if type(fuel) is list:
        multifuel = True
    elif type(fuel) is str:
        multifuel = False

    #Multifuel mixture percentage of total fuel
    # fuel is a list of fuels and their percentages in the total fuel
    #  odds are fuel name as a string
    #  evens are percetage of previous fuel name in total fuel
    # percentages should sum to 1 or script will not run
    if multifuel:
        check = 0
        for c in range(1, len(fuel), 2):
            check += fuel[c]
        if not np.isclose(check, 1):
            print('Error in Multifuel.'+
                  'Sum of individual fuel percentage must add up to 1!')
            sys.exit()

    #Multioxidizer mixture percentage of total oxidizer
    # fuel is a list of fuels and their percentages in the total fuel
    #  odds are fuel name as a string
    #  evens are percetage of previous fuel name in total fuel
    # percentages should sum to 1 or script will not run
    if multioxidizer:
        check    = 0
        for c in range(1, len(oxidizer), 2):
            check += oxidizer[c]
        if not np.isclose(check, 1):
            print('Error in Multioxidizer.'+
                  'Sum of individual fuel percentage must add up to 1!')
            sys.exit()

    conditions = {'Parameters': [Press, Temperature, mix_params, array_type],
                  'Mixture': [fuel, diluent, oxidizer],
                  'Flame': [mingrid, mul_soret, loglevel],
                  'Files': [mechanism, flame_temp],
                  'T/F': [multifuel, multioxidizer]}
    return conditions


def run_simulations(conditions, paramlist):
    """


    Parameters
    ----------
    conditions : TYPE
        DESCRIPTION.
    paramlist : TYPE
        DESCRIPTION.

    Returns
    -------
    flame_info_filtered : TYPE
        DESCRIPTION.
    flame_info_unfiltered : TYPE
        DESCRIPTION.
    sim_info : TYPE
        DESCRIPTION.

    """
    tic = time.time()
    chem = conditions['Files'][0]
    gas = ct.Solution(chem)

    print('Initial number of cases: '+format(len(paramlist)))
    print('\nStart of simulations...')
    sim_start  = time.time()
    flame_info = cf.parallelize(paramlist, conditions, flame_sens)
    sim_end    = time.time()
    sim_time   = sim_end - sim_start
    print('End of simulations')
    print('Simulations took '+format(sim_time, '0.5f')+' seconds.')
    print('\nStart flame information filtering...')
    filter_start = time.time()
    converged    = 0
    for x in flame_info:
        if x['Flame'][0] is None:
            continue
        else:
            converged += 1
    flame_info_unfiltered = copy.deepcopy(flame_info)
    duplicate_rxns = cf.duplicate_reactions(gas)
    flame_info_filtered = flame_info_filter(flame_info, duplicate_rxns)
    filter_end  = time.time()
    filter_time = filter_end - filter_start
    print('End of filtering')
    print('Filtering took '+format(filter_time, '0.5f')+ ' seconds.')
    print('\nNumber of cases converged:' +str(converged))
    toc      = time.time() #Main Code Time End
    duration = toc-tic
    print('Total time '+format(duration, '0.5f')+' seconds.\n')
    sim_info = [sim_time, converged, duration]
    return flame_info_filtered, flame_info_unfiltered, sim_info


def flame_sens(p, T, mix, cond):
    """
    Run one flame simulation.

    Parameters
    ----------
    p : float
        Pressure in atmospheres.
    T : float
        Temperature in Kelvin.
    mix : dict
        Dictionary where keys are mixture components, values are mole fracs.
    cond : dict
        Dictionary with detailed information on the condition.

    Returns
    -------
    flame_info : TYPE
        DESCRIPTION.

    """
    """[Fill in information]"""
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
    return flame_info


def flame_info_filter(flame_information, duplicate_reactions):
    """


    Parameters
    ----------
    flame_information : TYPE
        DESCRIPTION.
    duplicate_reactions : TYPE
        DESCRIPTION.

    Returns
    -------
    flame_information : TYPE
        DESCRIPTION.

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


    Parameters
    ----------
    cond : TYPE
        DESCRIPTION.
    fla_inf : TYPE
        DESCRIPTION.
    p_list : TYPE
        DESCRIPTION.
    s_info : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    save_time_start = time.time()
    p             = cond['Parameters'][0]
    T             = cond['Parameters'][1]
    mix_params    = cond['Parameters'][2]
    chem          = cond['Files'][0]
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    multif        = cond['T/F'][0]
    multio        = cond['T/F'][1]
    mg            = cond['Flame'][0]
    ms            = cond['Flame'][1]
    s_time        = datetime.timedelta(seconds=s_info[0])
    conv          = s_info[1]
    dur           = datetime.timedelta(seconds=s_info[2])

    #Save Path/Parent Directory
    parent_dir = 'Flame_Sensitivity_Results'
    try:
        os.makedirs(parent_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #Create Directory Name
    print('Creating Directory...')
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y_%m_%d %H.%M.%S Flame_Speed_Sens")
    directory = dt_string
    save_path = os.path.join(parent_dir, directory)
    os.makedirs(save_path)

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
                        cf.parameters_string(p, T, mix_params, chem,
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
    save_time_end = time.time()
    save_time = save_time_end - save_time_start
    print('Total File Save Time: '+format(save_time, '0.5f')+' seconds.\n')
    print('End file saving')


if __name__ == "__main__":
    ###########################Initializing####################################
    """[Fill in information]"""
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