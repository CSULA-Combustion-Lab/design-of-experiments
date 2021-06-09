#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 13:29:22 2018

@author: boss
"""

""" This is just a change to the discription to test
pushing a change in the code.
 """

import os
import math
import time
import errno
import pickle
import cantera
import shutil
import numpy as np
from tqdm import tqdm
from datetime import datetime
import common_functions as cf

cantera.suppress_thermo_warnings()

def run_0D_simulation(mech, arrtype, pres, temp, fue, oxi, dilu, m_type,
                      mix_params, safi, SpecificSpecies, Starttime, Endtime,
                      Delta_T, threshold, Save_time):
    """
    Takes information from initializer and runs necessary functions to perform
    a zero-dimensional simulation. Simulation results will be saved if booleans
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
    m_type : str
        A string providing the mixture type being used to create a mixture.
    mix_params : list
        A list of the two mixture parameters and mixtrue type used in creating
        a mixture.
    safi : bool
        If true simulation conditions and ranking results will be saved.
    SpecificSpecies : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated.
    Starttime : float
        Default set to zero in case the sim.time starts before 0.
    Endtime : float
        Sets the end time for the sim.time to be performed.
    Delta_T : float
        A float to check that the temperature did not change above the given
        limit. If the temperature change does exceed Delta_T then the simulation
        assumes a flame was created and the case is thrown out of the results.
    threshold : float
        A limit set by the experimental device that measures the mole fractions.
        If a species mole fraction is set below this limit the species mole
        fraction is set 0.
    Save_time : bool
        If true mole fraction sensitivities in time will be saved.

    Returns
    -------
    None.

    """

    #start time
    tic = time.time()
    mechan = cf.model_folder(mech)
    pack, paralist = initialization(mechan, arrtype, pres, temp,
                                    fue, oxi, dilu, m_type, mix_params,
                                    SpecificSpecies, Starttime, Endtime,
                                    Delta_T, float(threshold))
    zerod_info, siminfo, s_lists = run_simulations(pack, paralist)
    toc       = time.time()
    #end time
    duration  = toc-tic
    siminfo.append(duration)
    print('Number of cases ran', len(paralist))
    print('Time it took to run all simulations and rank them was '
          +format(duration, '0.5f')+' seconds.\n')

    if safi:
        file_saving(pack, s_lists, zerod_info, paralist, siminfo, Save_time)

def initialization(mechan, arrtype, pres, temp, fue, oxi, dilu, m_type,
                   mix_params, s_species, stime, etime, dT, ppm):
    """
    Creates a list and dictionary of simulation information that will be used
    the simulation calculations and ranking.

    Parameters
    ----------
    mechan : str
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
    m_type : str
        A string providing the mixture type being used to create a mixture.
    mix_params : list
        A list of the two mixture parameters and mixtrue type used in creating
    s_species : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated.
    Starttime : float
        Default set to zero in case the sim.time starts before 0.
    etime : float
        Sets the end time for the sim.time to be performed.
    dT : float
        A float to check that the temperature did not change above the given
        limit. If the temperature change does exceed dT then the simulation
        assumes a flame was created and the case is thrown out of the results.
    ppm : float
        A limit set by the experimental device that measures the mole fractions.
        If a species mole fraction is set below this limit the species mole
        fraction is set 0.

    Returns
    -------
    Packed : dict
        A dictionary of the following simulation information with the
        following structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture':[fuel_name, diluent_name, oxidizer_name],
             'ZeroD': [SpecificSpecies, dup_reactions],
             'Time_Info': [starttime, endtime, SCORE3_TIME],
             'Files': [mechanism],
             'Limits': [delta_T, ppm]}
    Parameter_List : list
        Simulation case information the following structure:
        [[Pressure, Temperature, Mixture], ...]

    """
    #create gas object
    gas                   = cantera.Solution(mechan)
    dup_reactions         = cf.duplicate_reactions(gas)
    Reactions             = range(0,len(gas.reaction_equations()))
    names                 = gas.species_names
    SpecificSpecieNumbers = [names.index(speci) for speci in s_species]
    SCORE3_TIME = stime + (etime - stime)*0.75 # time to evaluate score3

    Packed = {'Parameters': [pres, temp, mix_params, arrtype],
              'Mixture':[fue, dilu, oxi, m_type],
              'ZeroD': [s_species, dup_reactions,
                        Reactions, SpecificSpecieNumbers],
              'Time_Info': [stime, etime, SCORE3_TIME],
              'Files': [mechan],
              'Limits': [dT, ppm]}
    Parameter_List = cf.case_maker(Packed)
    return Packed, Parameter_List

def run_simulations(pack, plist):
    """
    Runs the simulation through all functions inlcuding the simulation
    calculations and the rankings.

    Parameters
    ----------
    pack : dict
        A dictionary of the following simulation information with the
        following structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture':[fuel_name, diluent_name, oxidizer_name],
             'ZeroD': [SpecificSpecies, dup_reactions],
             'Time_Info': [starttime, endtime, SCORE3_TIME],
             'Files': [mechanism],
             'Limits': [delta_T, ppm]}
    plist : list
        Simulation case information the following structure:
        [[Pressure, Temperature, Mixture], ...]

    Returns
    -------
    sim_package : list
        A list of simulation results for each case with the following
        sturcture:
            [t_T_P_AMol, t_SMol, t_AllSpecieSens, All_time_Sens_test,
            All_tTP_AMo_test, temp, mix, pressure]
    sim_info : list
        Simulation time and ranking time are recorded as floats.
    score_lists : list
        List of ranked simulation results with the following structure:
            [Parameters, score_T_P_MaxSensavg, scoretimes, rank_score,
             rank_plot, All_Ranks, All_Ranks_Params, score2_T_P_MaxSens,
             score2times, score2_Max_sens_rxn,
             score2_Params_MaxSens_Name_Params, score3_T_P_MaxSens,
             score3times, score3_Max_sens_rxn,
             score3_Params_MaxSens_Name_Params, int_strength,
             MoleFraction, All_time_Sens, All_tTP_AMo]

    """

    dT  = pack['Limits'][0]
    ppm = pack['Limits'][1]
    sspecies = pack['ZeroD'][0]
    rxns     = pack['ZeroD'][2]
    sspecnum = pack['ZeroD'][3]
    gas = cantera.Solution(pack['Files'][0])
    score3time = pack['Time_Info'][2]
    fuel     = pack['Mixture'][0]
    diluent  = pack['Mixture'][1]
    oxidizer = pack['Mixture'][2]

    #Lists
    Parameters      = []
    #first scoring criteria
    score_T_P_MaxSensavg = []
    scoretimes           = []
    rank_score           = []
    rank_plot            = []
    All_Ranks            = []
    All_Ranks_Params     = []
    #second scoring criteria
    score2_T_P_MaxSens  = []
    score2times         = []
    score2_Max_sens_rxn = []
    score2_Params_MaxSens_Name_Params = []
    #third scoring criteria
    score3_T_P_MaxSens  = []
    score3times         = []
    score3_Max_sens_rxn = []
    score3_Params_MaxSens_Name_Params = []
    # Fourth scoring criteria
    int_strength = []
    #mole/massfraction function
    MoleFraction  = []
    All_time_Sens = []
    All_tTP_AMo   = []

    print('Start of Simulations')
    sim_start   = time.time()
    sim_package = cf.parallelize(plist, pack, reac_sens_parallel)
    sim_end     = time.time()
    sim_time    = sim_end - sim_start
    print('End of Simulations')
    print('It took '+format(sim_time, '0.5f')+' seconds to simulate all cases')

    print('Start Ranking Process')
    rank_start = time.time()
    for sims in tqdm(sim_package):
        [t_T_P_AMol, t_SMol, t_AllSpecieSens,
          All_time_Sens_test, All_tTP_AMo_test,
          temp, mix, pressure] = sims
        temp_condition = T_limit(t_T_P_AMol, temp, dT)
        if temp_condition == 1:
            # if debug:
            #     print('T limit ' + DEBUG_FMT.format(temp, pressure, mix))
            del All_tTP_AMo_test[-1]
            del All_time_Sens_test[-1]
            continue
        else:
            All_tTP_AMo.append(All_tTP_AMo_test)
            All_time_Sens.append(All_time_Sens_test)
            gas.X = (mix)
            condition_info = {
                'temperature': temp, 'pressure': pressure*101.325,
                'mixture': mix, 'phi':gas.get_equivalence_ratio(),
                'fuel':cf.mixture_percentage(fuel,mix),
                'oxidizer':cf.mixture_percentage(oxidizer,mix),
                'diluent':cf.mixture_percentage(diluent,mix)}
            Parameters.append(condition_info)

        #initialize parameters
        SpecificSpecieSens = [] #list of sens of interest per rxn
        all_ranks          = []
        senstime           = np.array([x[0] for x in t_AllSpecieSens]) #times for all specie sensitivities
        sensitivities      = np.absolute(np.array([x[1] for x in t_AllSpecieSens])) #sensitivities
        molfrac_conditions = mole_fractions(t_SMol, pressure, temp, mix, ppm,
                                            MoleFraction, sspecies)

        specific_sens(sspecies, rxns, t_SMol, sspecnum, molfrac_conditions,
                      sensitivities, SpecificSpecieSens)

        sensitivity_score(SpecificSpecieSens, sspecies, temp, pressure, rxns,
                          senstime, score_T_P_MaxSensavg, scoretimes)

        sensitivity_score2(SpecificSpecieSens, sspecies, temp, pressure, mix,
                           rxns, senstime, gas, score2_T_P_MaxSens,
                           score2times, score2_Max_sens_rxn,
                           score2_Params_MaxSens_Name_Params,
                           fuel, oxidizer, diluent)

        sensitivity_score3(sspecies, mix, temp, pressure, score3time,
                           t_AllSpecieSens, SpecificSpecieSens, rxns, gas,
                           score3_Max_sens_rxn,
                           score3_Params_MaxSens_Name_Params,
                           fuel, oxidizer, diluent)

        rank_all(SpecificSpecieSens, temp, pressure, mix, rxns, senstime,
                 sspecies, all_ranks, All_Ranks, All_Ranks_Params)

        IS_norm = integrated(t_AllSpecieSens, sspecnum, len(rxns),
                              molfrac_conditions)
        int_strength.append([condition_info, IS_norm])

    score_lists = [Parameters, score_T_P_MaxSensavg, scoretimes, rank_score,
                    rank_plot, All_Ranks, All_Ranks_Params, score2_T_P_MaxSens,
                    score2times, score2_Max_sens_rxn,
                    score2_Params_MaxSens_Name_Params, score3_T_P_MaxSens,
                    score3times, score3_Max_sens_rxn,
                    score3_Params_MaxSens_Name_Params, int_strength,
                    MoleFraction, All_time_Sens, All_tTP_AMo]

    rank_end  = time.time()
    rank_time = rank_end - rank_start
    print('End Ranking Process')
    print('It took '+format(rank_time, '0.5f')+' seconds to rank all cases')
    sim_info = [sim_time, rank_time]
    return sim_package, sim_info, score_lists


def reac_sens_parallel(pressure, temp, mix, pack):
    """ Function which saves simulation history.

    This function saves the different parameters through each time step
    per simulation. Time, Temperature, Pressure, Mole Fractions,
    Sensitivities. It reduces the amount of information through the
    calling of the 'reduce_resolution' function. Lastly, it seperates all
    information into different lists to facilitate other calculations.

    Dummy Variables
    ---------------
    row : list
        List which collects all information from current time step.
    rows1 : list
        List which appends all information from 'row' to pass on to
        'reduce resolution' function
    rows : list
        List which contains the reduced amount of points returned from
        the 'reduce_resolution' function

    Appends
    --------
    t_T_P_AMol : list
        List which contains time, Temperature, Pressure, all species mole
        fractions per simulation.
    t_SMol : list
        List which contains time, and 'SpecificSpecies' mole fractions
        per simulation.
    t_AllSpeciesSens : list
        List which contains time, and all sensitivities of species with
        respect to reactions available in GRI-Mech 3.0, an optimized
        mechanism designed to model natural gas combustion, including
        NO formation and reburn chemistry. List is unique per simulation.
    All_time_Sens : list
        List which saves time and all species sensitivities with respect
        to reactions available in GRI-Mech 3.0 mechanism file. This list
        is a list of lists, in which each list corresponds to a different
        simulation or set of initial conditions.
    All_tTP_AMo : list
        This is another list which saves time, Temperature, Pressure, and
        all species mole fractions. This list also is a list of lists in
        which each individual list corresponds to a different simulation
        or set of initial conditions.

    """

    specificspecies = pack['ZeroD'][0]
    mech            = pack['Files'][0]
    dup_reactions   = pack['ZeroD'][1]
    starttime       = pack['Time_Info'][0]
    endtime         = pack['Time_Info'][1]
    SCORE3_TIME     = pack['Time_Info'][2]

    gas1 = cantera.Solution(mech)
    gas1.TPX = temp, pressure*101325, mix
    rxns     = range(0,len(gas1.reaction_equations()))

    reac = cantera.IdealGasConstPressureReactor(gas1, name=None, energy='on')
    sim  = cantera.ReactorNet([reac])

    for i in rxns:
        reac.add_sensitivity_reaction(i)

    row                = [None]*6
    rows1              = []
    rows               = []

    while sim.time < endtime and sim.time >= starttime:
        row[0]   = sim.step() #save time
        row[1:3] = gas1.TP    #save T, P
        row[3]   = gas1[gas1.species_names].X #save mole fraction value of all species
        row[4]   = gas1[specificspecies].X    #save molesfrac value of specific species
        ssensitivities_uf = sim.sensitivities()  #sensitivities
        ssensitivities = sens_dup_filter(ssensitivities_uf[2:], dup_reactions)
        row[5]   = ssensitivities #first two rows are mass and enthalpy or temperature, not sensitivities
        rows1.append([x for x in row])
    rows = reduce_resolution(rows1,100, SCORE3_TIME) #isnt rows 1 always going to have a length of 1?
    pts  = len(rows)
    t_T_P_AMol      = [rows[i][0:4] for i in range(0,pts)]
    t_SMol          = [[rows[i][0],rows[i][4]] for i in range(0,pts)]
    t_AllSpecieSens = [[rows[i][0],rows[i][5]] for i in range(0,pts)]
    All_time_Sens_test = t_AllSpecieSens
    All_tTP_AMo_test = t_T_P_AMol

    package = [t_T_P_AMol, t_SMol, t_AllSpecieSens, All_time_Sens_test,
               All_tTP_AMo_test, temp, mix, pressure]

    return package


def reduce_resolution(mylist, maxlength, time_of_interest):
    """ Reduces the number of elements in a list if the list length is greater
    than the specificed maximum. The function saves every nth term from the list
    to a new list, where n is calculated to be an optimal sequential number to
    evenly represent the original list. The resulting number of elements will be
    less than or equal to the speicified maximum.

    Parameters
    ----------
    mylist : list
        List for which reduced resolution is desired. 1st entry of each row
        must be time in seconds.
    maxlength : int
        Integer greater than zero which sets the maximum number of elements
        desired in the list.
    time_of_interest : float
        Include the first row after this time.

    Returns
    ---------
    reduced : list
        List with a reduced resolution less than or equal to the maxlength.

    """

    reduced = []
    length  = len(mylist)
    if length > maxlength:
        time = np.array([x[0] for x in mylist])
        ind = np.argwhere(time > time_of_interest)[0][0]
        nth     = math.ceil(length/maxlength)
        reduced = mylist[:ind:nth]
        if reduced[-1] != mylist[ind]:
            reduced.append(mylist[ind])
        reduced.extend(mylist[ind+1::nth])
    else:
        reduced = mylist
    return reduced


def T_limit(t_T_P_AMol,temp,delta_T):
    """
    Determines if the temperature change per time-step is greater than the
    maximum allowable temperautre change

    Parameters
    ----------
    t_T_P_AMol : list
        Temperatures of the previous time-step.
    temp : float
        Temperature of the current time-step.
    delta_T : int or float
        Maximum allowable tempearture change.

    Returns
    -------
    temp_condition : int
        A integer of either 1 or 0. If 0 the temperature change did not surpass
        the limit. If 1 the temperatue did surpass the limit.

    """

    max_temp = t_T_P_AMol[-1][1] #temp of the last time step for condition
    if max_temp-temp > delta_T:
        temp_condition = 1
    else:
        temp_condition = 0
    return temp_condition


def mole_fractions(t_SMol, pressure, temp, mix,
                   ppm, MoleFraction, specificspecies):
    """Function which zeros mole fractions values.

    This function is checking if the mole fractions from the t_SMol
    list are above a predefined ppm (part per million). In this case,
    if the mole fractions are not above one ppm then the values get
    replaced with zeros. This will facilitate future analysis and
    ranking of the sensitivities to reactions of interest.

    Parameters
    ----------
    SpecificSpecies : list
        List created by user to identify Species of interest.
    molfrac_time : list
        List which contains all of the time steps for the simulation
        being considered.
    molfrac : list
        List which contains all of the mole fractions relating to each
        time step in the 'molfrac_time' list for 'SpecificSpecies.'
    ppm : float
        Floating point number which will define how large the mole
        fractions should be, in ppm.

    Dummy Variables
    ---------------
    row : list
        Dummy variable used to replace actual mole fraction values with
        zero if they are not above the predefined 'ppm' value

    Appends
    -------
    molfrac_conditions : list
        List which contains current simulation mixture, Temperature,
        Pressure, mole fractions, and time steps. This list contains
        the updated mole fraction values after checking if they are
        above the predefined ppm value. If they are above the ppm value,
        they are left alone, if they are not they are replaced with 0.
    MoleFraction : list
        This is a list of lists. Each list is identical to the
        information in 'molfrac_conditions' and each list corresponds
        to a different simulation or initial conditions.

    """
    molefrac_time = np.array([x[0] for x in t_SMol])
    molfrac       = np.absolute(np.array([x[1] for x in t_SMol]))  #specific specie moles frac

    molfrac_conditions      = [None]*(len(specificspecies)+2)
    molfrac_conditions[0:2] = [mix, temp, pressure]

    for i in range(0,len(specificspecies)):
        row = np.zeros((len(molfrac),2))
        for j in range(0,len(molfrac)):
            if molfrac[j,i] >= ppm:
                row[j,0] = molfrac[j,i]
                row[j,1] = molefrac_time[j]
            else:
                row[j,0] = 0
                row[j,1] = molefrac_time[j]
        molfrac_conditions[i+3] = row
    MoleFraction.append([x for x in molfrac_conditions])
    return molfrac_conditions


def specific_sens(specificspecies, rxns, t_SMol, SpecificSpecieNumbers,
                  molfrac_conditions, sensitivities, SpecificSpecieSens):
    """ Function which zeroes sensitivity values.

    This function iterates through all of the 'SpecificSpecie'
    'molefrac_conditions', if the mole fractions are above the predefined
    ppm value within the 'mole_fractions' function then the sensitivities
    get saved, if not the sensitivities are replaced with value of 0.

    Parameters
    ----------
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural
        gas combustion, including NO formation and reburn chemistry.
    SpecificSpecies : list
        List created by user to identify Species of interest.
    t_SMol : list
        List of lists containing time as the first element, and the
        corresponding mole fractions, of the chosen SpecificSpecies, as a
        list for the second element.
    SpecificSpecieNumbers : list
        List containing the index number of species in 'SpecificSpecies'
        list with respect to the used mechanism file. (GRI - Mech 3.0).
    molfrac_conditions : list
        Please refer to 'mole_fractions' function.
    sensitivities : list
        List containing sensitivity values for all species in GRI - Mech
        mechanism file with respect to each reaction within the same
        mechanism file.

    Dummy Variables
    ---------------
    row : list
        Dummy variable used for iteration purposes and for appending
        SpecificSpecieSens.
    molfractions : list
        List which grabs the mole fractions and time steps from the
        molfrac_conditions list.
    MolFractions : list
        List which grabs the mole fraction values from molfractions list
        defined above.
    MolFractionsTime : list
        List which grabs the time step values from the molfractions list
        defined above.

    Appends
    -------
    SpecificSpecieSens: list
        List which contains the sensitivities of 'SpecificSpecies.' This
        list is dependent on the 'mole_fractions' function, if the mole
        fractions are not above the predefined ppm value then the
        sensitivites get replaced with zero, otherwise they are
        not changed. This list is unique and depends on initial
        conditions or mixture of current simulation. The size of this
        list should be len(SpecificSpecies)*len(rxns).
        Format of this list: It contains one list for each time step, so
        SpecificSpecieSens[0] and SpecificSpecieSens[1] are different tiems.
        These lists are then ordered as sensitivity of: [species1 to rxn1,
        species2 to rxn 1, species3 to rxn 1.... species1 to rxn2, species2
        to rxn2....]

    """
    row = [None]*len(specificspecies)*len(rxns)
    molefrac_time = np.array([x[0] for x in t_SMol])
    for i in range(0,len(molefrac_time)):
        for k in range(0,len(rxns)):
            for j in range(0,len(SpecificSpecieNumbers)):
                molfractions = molfrac_conditions[3+j]
                MolFractions = molfractions[:,0]
#                MolFractionsTime = molfractions[:,1]
                if MolFractions[i] == 0:
                    row[j+k*len(specificspecies)] = 0
                elif MolFractions[i] != 0:
                    row[j+k*len(specificspecies)] = sensitivities[i,SpecificSpecieNumbers[j],k]
        SpecificSpecieSens.append([x for x in row])
    return SpecificSpecieSens


def sensitivity_score(SpecificSpecieSens, specificspecies, temp, pressure,
                      rxns, senstime, score_T_P_MaxSensavg, scoretimes):
    """Ratio of sum of all maximum reaction values to number of reactions.

    This function obtains the maximum sensitivity from a complete simulation
    per specie, per reaction, and takes the average of all maximum values
    with respect to the number of reactions within the used mechanism file.

    Parameters
    ----------
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural
        gas combustion, including NO formation and reburn chemistry.
    SpecificSpecies : list
        List created by user to identify Species of interest.

    Dummy Variables
    ---------------
    dataa : ndarray
        N dimensional array which contains a list of lists. Each list
        corresponds to a different time step within a specific simulation.
        Within each list, there are sensitivity values for each
        'SpecificSpecies' for each reaction within the used mechanism file.
        The size of each list should be len(SpecificSpecies)*len(rxns)
    row : list
        List which contains Temperature, Pressure and the average of the
        found maximum sensitivities.
    row1 : list
        List which contains the average time at which all maximum
        sensitivities were found.
    rxn_maxs : ndarray
        N - dimensional array used to collect maximum sensitivities per
        reaction per species throughout all time steps in current
        simulation results being analyzed.
    rxn_t : ndarray
        N - dimensional array used to store the time at which the maximum
        sensitivities from 'rxn_maxs' occur in the simulation results.

    Appends
    -------
    score_T_P_MaxSensavg : list
       A list in which elements 0 and 1 corresponds to Temperature,
       Pressure, and elements 2-6 correspond to the ratio of the sum
       of all maximum sensitivites to the number of reactions in the
       mechanism file being used.
    scoretimes : list
        A list which relates time to elements 2-len(SpecificSpecies) of
        score_T_P_MaxSensavg.

    """
    dataa = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    row   = [None]*(len(specificspecies)+2)
    row1  = [None]*(len(specificspecies))

    row[0] = temp             #K
    row[1] = pressure*101.325 #kPa

    rxn_maxs = np.zeros(len(rxns))
    rxn_t    = np.zeros(len(rxns))

    for i in range(0,len(specificspecies)):
        for j in range(0,len(rxns)):
            rxn_maxs[j] = max(dataa[:,i+j*len(specificspecies)])
            rxn_t[j]=senstime[np.argmax(dataa[:,i+j*len(specificspecies)])]
        row[i+2] = sum(rxn_maxs)/len(rxns) #avg of max sensitivity per reactions
        row1[i]  = sum(rxn_t)/len(rxns)
    score_T_P_MaxSensavg.append([x for x in row]) #T,P,avg between maximum sens per reaction
    scoretimes.append([x for x in row1]) #what does this populate?


def sensitivity_score2(SpecificSpecieSens, specificspecies, temp, pressure,
                       mix, rxns, senstime, gas, score2_T_P_MaxSens,
                       score2times, score2_Max_sens_rxn,
                       score2_Params_MaxSens_Name_Params,
                       fuel, oxidizer, diluent):
    """
    Find the reaction with the highest absolute value of sensitivity for each species of interest.
    The maximum sensitivity is found for all time - not at a specific time step.
    This score will find spikes, but may not find a reaction that is sensitive for a longer time.

    Parameters
    ----------
    SpecificSpecieSens : list
        List which contains the sensitivities of 'SpecificSpecies.' This
        list is dependent on the 'mole_fractions' function, if the mole
        fractions are not above the predefined ppm value then the
        sensitivites get replaced with zero, otherwise they are
        not changed. This list is unique and depends on initial
        conditions or mixture of current simulation. The size of this
        list should be len(SpecificSpecies)*len(rxns).
        Format of this list: It contains one list for each time step, so
        SpecificSpecieSens[0] and SpecificSpecieSens[1] are different tiems.
        These lists are then ordered as sensitivity of: [species1 to rxn1,
        species2 to rxn 1, species3 to rxn 1.... species1 to rxn2, species2
        to rxn2....]
    specificspecies : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated.
    temp : float
        Temperature in Kelvin of the current case.
    pressure : float
        Pressure in atmospheres of the current case.
    mix : list
        Molefraction of species in the mixture.
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural
        gas combustion, including NO formation and reburn chemistry.
    senstime : array
        Times for all specie sensitivities.
    gas : object
        Cantera generated gas object created using user provided mechanism.
    score2_T_P_MaxSens : list
        Thermodynamic properties and maximum sensitivity.
    score2times : list
        Times of maximum sensitivity occurances.
    score2_Max_sens_rxn : list
        Simulation information appended as well as the reaction with the
        maximum sensitivity
    score2_Params_MaxSens_Name_Params : list
        Simulation information stored as a dictionary as well as a list of
        reaction number and equation that was most sensitive.
    fuel : str or list
        As a string the variable represents a single species of fuel being used.
        As a list the variable represents multicomponent fuel species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    oxidizer : str or list
        As a string the variable represents a single species of oxidizer being used.
        As a list the variable represents multicomponent oxidizer species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    diluent : str or list
        As a string the variable represents a single species of diluent being used.
        As a list the variable represents multicomponent diluent species
        followed by the percentage to the total diluent [Component1, % of total, ...]

    Returns
    -------
    None.

    """
    dataa2  = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    ss      = [None]*(len(specificspecies)+2)  #ss-sensititvity score
    ss_time = [None]*(len(specificspecies))
    max_rxn = [None]*(len(specificspecies)+3)

    ss[0] = temp
    ss[1] = pressure*101.325 #in kPa for the plot
    max_rxn[0] = mix
    max_rxn[1] = temp
    max_rxn[2] = pressure*101.325
    rxn_name   = [None]*(len(specificspecies)+1)
    rxn_name[0] = {'temperature': temp, 'pressure': pressure*101.325,
                   'mixture': mix, 'phi':gas.get_equivalence_ratio(),
                   'fuel':cf.mixture_percentage(fuel,mix),
                   'oxidizer':cf.mixture_percentage(oxidizer,mix),
                   'diluent':cf.mixture_percentage(diluent,mix)}

    max_sens = np.zeros(len(rxns))
    rxn_t    = np.zeros(len(rxns))

    for i in range(0,len(specificspecies)):
        for j in range(0,len(rxns)):
            max_sens[j] = max(dataa2[:,i+j*len(specificspecies)])
            rxn_t[j] = senstime[np.argmax(dataa2[:,i+j*len(specificspecies)])]
        ss[i+2]      = max(max_sens)              #maximum sensitivity
        ss_time[i]   = rxn_t[np.argmax(max_sens)] #time of max sensitivity
        max_rxn[i+3] = rxns[np.argmax(max_sens)]  #rxn with max sensitivity
        rxn_name[i+1] = gas.reaction_equations([rxns[np.argmax(max_sens)]])

    score2_T_P_MaxSens.append([x for x in ss])
    score2times.append([x for x in ss_time])
    score2_Max_sens_rxn.append([x for x in max_rxn])
    score2_Params_MaxSens_Name_Params.append([x for x in rxn_name])


def sensitivity_score3(specificspecies, mix, temp, pressure, SCORE3_TIME,
                       t_AllSpecieSens, SpecificSpecieSens, rxns, gas,
                       score3_Max_sens_rxn, score3_Params_MaxSens_Name_Params,
                       fuel, oxidizer, diluent):
    """
    Determines what reactions was most sensitivty and appends the simulation
    case information, reaction information, and simulation results

    Parameters
    ----------
    specificspecies : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated.
    mix : list
        Molefraction of species in the mixture.
    temp : float
        Temperature in Kelvin of the current case.
    pressure : float
        Pressure in atmospheres of the current case.
    SCORE3_TIME : float
        A specified time near the end of the simulation case to capture
        information about the reaction about.
    t_AllSpecieSens : list
        All species sensitivty information including thermodynamic information,
        molefractions, and time-steps.
    SpecificSpecieSens : list
        List which contains the sensitivities of 'SpecificSpecies.' This
        list is dependent on the 'mole_fractions' function, if the mole
        fractions are not above the predefined ppm value then the
        sensitivites get replaced with zero, otherwise they are
        not changed. This list is unique and depends on initial
        conditions or mixture of current simulation. The size of this
        list should be len(SpecificSpecies)*len(rxns).
        Format of this list: It contains one list for each time step, so
        SpecificSpecieSens[0] and SpecificSpecieSens[1] are different tiems.
        These lists are then ordered as sensitivity of: [species1 to rxn1,
        species2 to rxn 1, species3 to rxn 1.... species1 to rxn2, species2
        to rxn2....]
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural
        gas combustion, including NO formation and reburn chemistry.
    gas : object
        Cantera generated gas object created using user provided mechanism.
    score3_Max_sens_rxn : list
        Simulation information appended as well as the reaction with the
        maximum sensitivity
    score3_Params_MaxSens_Name_Params : list
        Simulation information stored as a dictionary as well as a list of
        reaction number and equation that was most sensitive.
    fuel : str or list
        As a string the variable represents a single species of fuel being used.
        As a list the variable represents multicomponent fuel species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    oxidizer : str or list
        As a string the variable represents a single species of oxidizer being used.
        As a list the variable represents multicomponent oxidizer species
        followed by the percentage to the total fuel [Component1, % of total, ...]
    diluent : str or list
        As a string the variable represents a single species of diluent being used.
        As a list the variable represents multicomponent diluent species
        followed by the percentage to the total diluent [Component1, % of total, ...]

    Returns
    -------
    None.

    """
    num_spec = len(specificspecies)
    max_rxn = [None]*(num_spec+4)
    max_rxn[0]  = mix
    max_rxn[1]  = temp
    max_rxn[2]  = pressure*101.325
    max_rxn[3]  = gas.get_equivalence_ratio()
    rxn_name    = [None]*(num_spec+1)
    rxn_name[0] = {'temperature': temp, 'pressure': pressure*101.325,
                   'mixture': mix, 'score3time': SCORE3_TIME,
                   'phi':gas.get_equivalence_ratio(),
                   'fuel':cf.mixture_percentage(fuel,mix),
                   'oxidizer':cf.mixture_percentage(oxidizer,mix),
                   'diluent':cf.mixture_percentage(diluent,mix)}

    time = np.array([x[0] for x in t_AllSpecieSens])
    ind = np.argwhere(time > SCORE3_TIME)[0][0]
    sens_at_t = SpecificSpecieSens[ind]

    for i in range(num_spec):
        unshuffled = sens_at_t[i::num_spec]
        max_sens_ind = np.argmax(unshuffled)
        if unshuffled[max_sens_ind] == 0:  # Below ppm limit
            max_rxn[i+4] = None
            rxn_name[i+1] = None
        else:
            max_rxn[i+4] = rxns[max_sens_ind]  #rxn with max sensitivity
            rxn_name[i+1] = [rxns[max_sens_ind],
                             gas.reaction_equations([rxns[max_sens_ind]])]

    score3_Max_sens_rxn.append([x for x in max_rxn])
    score3_Params_MaxSens_Name_Params.append([x for x in rxn_name])


def rank_all(SpecificSpecieSens, temp, pressure, mix, rxns, senstime,
             specificspecies, all_ranks, All_Ranks, All_Ranks_Params):
    """
    Creates a list of reaction rankings for all timesteps and all species.
    Needs work. Understanding output rank_mat: The lists cycle though the
    species list, and after every len(species) its a new time_step.

    Parameters
    ----------
    SpecificSpecieSens : list
        List which contains the sensitivities of 'SpecificSpecies.' This
        list is dependent on the 'mole_fractions' function, if the mole
        fractions are not above the predefined ppm value then the
        sensitivites get replaced with zero, otherwise they are
        not changed. This list is unique and depends on initial
        conditions or mixture of current simulation. The size of this
        list should be len(SpecificSpecies)*len(rxns).
        Format of this list: It contains one list for each time step, so
        SpecificSpecieSens[0] and SpecificSpecieSens[1] are different tiems.
        These lists are then ordered as sensitivity of: [species1 to rxn1,
        species2 to rxn 1, species3 to rxn 1.... species1 to rxn2, species2
        to rxn2....]
    temp : float
        Temperature in Kelvin of the current case.
    pressure : float
        Pressure in atmospheres of the current case.
    mix : list
        Molefraction of species in the mixture.
    rxns : ndarray
        Array containing the number of available Chemical Reactions based
        on GRI-Mech 3.0, an optimized mechanism designed to model natural
        gas combustion, including NO formation and reburn chemistry.
    senstime : array
        Times for all specie sensitivities.
    specificspecies : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated.
    all_ranks : list
        An empty list filled with ranking information from ranking_per_step
    All_Ranks : list
        A list of list of all_ranks information per case.
    All_Ranks_Params : list
        A list of dictinaries of simulation information per case.

    Returns
    -------
    None.

    """
    data3 = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    all_ranks_params    = [None]*1
    all_ranks_params[0] = {'temperature': temp, 'pressure': pressure*101.325,
                         'mixture': mix}
    step  = [None]*len(rxns)
    for i in range(0,len(senstime)):
        for j in range(0,len(specificspecies)):
            for k in range(0,len(rxns)):
                step[k] = data3[i,j+k*len(specificspecies)] #creates list of specie sensitivities per rxns, per time, per specie
            ranking = ranking_per_step(step)
            all_ranks.append([x for x in ranking])
            all_ranks_params.append([x for x in ranking])

    All_Ranks.append([x for x in all_ranks])
    All_Ranks_Params.append([x for x in all_ranks_params])


def ranking_per_step(sval):
    """
    Assigns a rank to each reaction sensitivity per time-step and species

    Parameters
    ----------
    sval : list
        A list of specie sensitivites per reactions.

    Returns
    -------
    step_rank : list
        Rank given to each reaction sensitivity per time-step and species.
    """

    sval  = np.array(sval)   #sval = step from above per time, per specie
    order = np.argsort(sval) #orders sval from smallest to largeest by listing indices which sval elements should be arranged in
    ranks = ranks_possible(sval)
    step_rank = [None]*len(sval)
    for j in range (0,len(sval)):
        step_rank[order[j]] = ranks[j]
    return step_rank


def ranks_possible(sval):
    """
    Determines how many ranks there will be for a given time-step, since
    sensitivities can tie if they are equal

    Parameters
    ----------
    sval : list
        A list of specie sensitivites per reactions..

    Returns
    -------
    ranks : list
        A list of possible ranks.

    """
    sval  = np.sort(sval)        #orders elements from smallest to largest
    sval  = np.flip(sval,axis=0) #re-orders from largest to smallest
    ranks = [None]*len(sval)
    if sval[0] == 0:
        ranks[0] = len(sval)
    else:
        ranks[0] = 1
    for i in range (1,len(sval)):
        if sval[i] == 0:
            ranks[i] = len(sval)
        else:
            if sval[i] == sval[i-1]:
                ranks[i] = ranks[i-1]
            else:
                ranks[i] = ranks[i-1]+1
    ranks = ranks[::-1] #Flips the list
    return ranks


def integrated(t_sens, spec_nums, num_rxns, mole_frac):
    """ Calculated normalized integrated strength for certain species.

    Parameters
    ----------
    t_sens : list
        Structure is a list of lists, where each list represents one time step.
        Within each list, the first item is the time [s], and the second item
        is an array of sensitivities where each row is a species and each
        column is a reaction.
    spec_nums : list
        The species numbers of interest (0-indexed)
    mole_frac : list
        Created in mole_fractions() and returned as molfrac_conditions. This
        list indicates where the species of interest mole fraction is below the
        ppm limit.

    Returns
    -------
    IS_norm : list
        This is a list of lists. Each list is the result for one species.
        Within each species is a list of normalized integrated strength for each
        reaction in order. In the future, it may be beneficial to also return
        IS, the non-normalized integrated strength.
    """
    t = 0
    num_spec = len(spec_nums)
    IS = num_spec * [np.zeros(num_rxns)]
    mole_frac_rearranged = np.hstack((mole_frac[3][:, 1][:, None],
                                      mole_frac[3][:, 0][:, None]))
    if num_spec > 1:
        # Check that the time-basis is the same for all
        for i in range(1, num_spec):
            assert np.all(np.isclose(mole_frac[3][:, 1], mole_frac[3+i][:, 1]))
            mole_frac_rearranged = np.hstack((mole_frac_rearranged,
                                              mole_frac[3+i][:, 0][:, None]))

    # mole_frac_rearranged is now a big array where the first column is time,
    # and the other columns are mole fractions of each species of interest.
    for t_index in range(len(t_sens)):  # For every time step
        dt = t_sens[t_index][0] - t
        t = t_sens[t_index][0]
        assert np.isclose(t, mole_frac_rearranged[t_index, 0])
        if dt == 0:
            continue
        for i in range(num_spec):  # For each species
            if mole_frac_rearranged[t_index, i+1] != 0:  # check that species is detectable
                spec = spec_nums[i]  # Get the species index
                sens = t_sens[t_index][1][spec, :]  # Sensitivity of this species to all rxns at this time
                IS[i] = IS[i] + sens * dt

    # Normalize it
    IS_norm = []
    for spec in IS:
        top5 = spec[np.argsort(np.abs(spec))[-5:]]
        if all([x == 0 for x in spec]):
            IS_norm.append(None)
        else:
            IS_norm.append(spec / np.mean(np.abs(top5)))

    return IS_norm


def sens_dup_filter(sens_information, duplicate_reactions):
    """
    Sums the sensitivites of duplicate reactions.

    Parameters
    ----------
    sens_information : List
        List of sensitivites per reaction.
    duplicate_reactions : Dictionary
        Dictionary of duplicate reactions from the mechanism file.

    Returns
    -------
    sens_information : list
        A reformated list of parameter sens_information which represents the
        sensitivity analysis of all reactions. Any duplicate reactions are
        summed together into one reaciton.

    """
    for s in sens_information:
        for d in duplicate_reactions.keys():
                    sumt = 0
                    for n in duplicate_reactions[d]:
                         sumt += s[n]
                         s[n] = 0
                    s[duplicate_reactions[d][0]] = sumt
    return sens_information


def file_saving(pack, slists, zerod_info, plist, siminfo, sati):
    """
    Information from the simulation is pickled and saved in a generated folder
    with a text document filled with details of the simulation performed.

    Parameters
    ----------
    pack : dict
        A dictionary of the following simulation information with the
        following structure:
            {'Parameters': [Press, Temperature, mix_params, array_type],
             'Mixture':[fuel_name, diluent_name, oxidizer_name],
             'ZeroD': [SpecificSpecies, dup_reactions],
             'Time_Info': [starttime, endtime, SCORE3_TIME],
             'Files': [mechanism],
             'Limits': [delta_T, ppm]}
    slists : list
        List of ranked simulation results with the following structure:
            [Parameters, score_T_P_MaxSensavg, scoretimes, rank_score,
             rank_plot, All_Ranks, All_Ranks_Params, score2_T_P_MaxSens,
             score2times, score2_Max_sens_rxn,
             score2_Params_MaxSens_Name_Params, score3_T_P_MaxSens,
             score3times, score3_Max_sens_rxn,
             score3_Params_MaxSens_Name_Params, int_strength,
             MoleFraction, All_time_Sens, All_tTP_AMo]
    zerod_info : list
        A list of simulation results for each case with the following
        sturcture:
            [t_T_P_AMol, t_SMol, t_AllSpecieSens, All_time_Sens_test,
            All_tTP_AMo_test, temp, mix, pressure]
    plist : list
        Simulation case information the following structure:
        [[Pressure, Temperature, Mixture], ...].
    siminfo : list
        Simulation time and ranking time are recorded as floats
    sati : boolean
        If true molefraction time sensitivity information is pickeled.
        Note: Time senstivity information list may become to large to pickle.

    Returns
    -------
    None.

    """

    P           = pack['Parameters'][0]
    T           = pack['Parameters'][1]
    m_params    = pack['Parameters'][2]
    a_type      = pack['Parameters'][3]
    sspecies    = pack['ZeroD'][0]
    rxns        = pack['ZeroD'][2]
    mech        = pack['Files'][0]
    SCORE3_TIME = pack['Time_Info'][2]
    dT          = pack['Limits'][0]
    ppm         = pack['Limits'][1]
    f_name      = pack['Mixture'][0]
    d_name      = pack['Mixture'][1]
    o_name      = pack['Mixture'][2]


    print('Creating Directory...')
    #Save Path/Parent Directory
    parent_dir = '0D_Sensitivity_Results'
    try:
        os.makedirs(parent_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    #Create Directory Name
    now = datetime.now()
    dt_string = now.strftime("%Y_%m_%d %H.%M.%S")
    directory = dt_string
    save_path = os.path.join(parent_dir, directory)
    os.makedirs(save_path)

    figures_dir = 'ZeroD_Sensitivity_Plots'
    figure_path = os.path.join(save_path, figures_dir)
    os.makedirs(figure_path)

    print('\nDirectory Created\nCreating Text File...')
    #Text Description
    filename  = 'Case Description.txt'
    filename  = os.path.join(save_path, filename)
    f         = open(filename, "w")
    text_description = ("This file provides a description of the parameters "
            "and specific species of interest investigated in the "
            "Flow Reactor Simulation.\r\n"
            "=======Species of Interest======="
            "\r"+format(sspecies)+"\r"
            "=================================\r\n" +
            cf.parameters_string('0D', P, T, m_params, mech, f_name, o_name, d_name)
            + "\n"
            "==========Number of Cases==========\r"
            "Total Cases Simulated 	     = "+format(len(plist))+"\r"
            "Cases Meeting All Conditions = "+format(len(slists[5]))+"\r"
            "===================================\r\n"
            "=====Conditions=====\r"
            "PPM Limit     = "+format(ppm*1e6)+"\r"
            "delta_T Limit = "+format(dT)+"\r"
            "====================\r\n"
            "===============Times===============\r"
            "Simulation Time = "+format(siminfo[0], '0.2f')+" [seconds]\r"
            "Ranking Time    = "+format(siminfo[1], '0.2f')+" [seconds]\r"
            "Total Time      = "+format(siminfo[2], '0.2f')+" [seconds]\r"
            "===================================")
    f.write(text_description)
    f.close()

    print('\nText File Created\nSaving Files...\n')
    #File saving
    gas_reactions = cantera.Solution(mech).reaction_equations()
    species_rxn_file = [sspecies, [len(rxns)], [gas_reactions],
                        [SCORE3_TIME], [a_type]]

    # This loop can make the saving easier. Just put the thing to be
    # saved and the file name in the save_list
    save_list = [(slists[15], 'Integrated Strength.pkl'),
                 (slists[6], 'All_Ranks.pkl'),
                 (slists[14], 'Max_Sens_Rxn.pkl'),
                 (species_rxn_file, 'Species_Rxn.pkl'),
                 (slists[0], 'Case_Parameters.pkl')]
    for item in save_list:
        file = os.path.join(save_path, item[1])
        with open(file, 'wb') as f:
            pickle.dump(item[0], f)

    shutil.copyfile('input.yaml', os.path.join(save_path, 'input.yaml'))
    if sati:
        Mole_Frac_fname = 'Mole_fractions.pkl'
        file = os.path.join(save_path, Mole_Frac_fname)
        with open(file, "wb") as f:
            pickle.dump(slists[18], f)
        f.close()

        file = os.path.join(save_path, 'All_time_sens.pkl')
        with open(file, "wb") as f:
            pickle.dump(slists[17], f)
        f.close()
        print('Files Saved')

    else:
        print('File Saved')

if __name__ == "__main__":
    ####Set experiment parameters
    mechanism = 'mech-FFCM1_modified.cti' #Mechanism file
    #Parameters for main loop
    Temperature = [600, 2500, 2]    #Temperature [K]
    Press       = [.1, 100, 2]      #Pressure [atm]
    Phi         = [0.00025, 2.5, 3] #Equivalence ratio
    Fuel        = [0.00001, 0.1, 3] #Fuel mole fraction
    Dilper      = [0.00001, 0.1, 4]

    #Parameters for mixture
    fuel_name = 'H2' #chemical formula of fuel
    diluent_name = 'N2' #chemical formula of diluent
    oxidizer_name = 'O2'
    Mixture_type = 'phi_fuel'
    array_type = 'lin'
    Mix_params = (Mixture_type, Phi, Fuel)

    SpecificSpecies = ['OH'] #Species of interest for rxn ranking data
    Starttime = 0     #in the case a reading is too early
    Endtime   = 0.001 #one milisecond
    Delta_T   = 100
    PPM       = 1/1000000 #one ppm

    save_files = True # If true, save files for plotting script
    save_time  = True # If true, also save files for GUI. Not recommended for large runs
    debug      = True  # If true, print lots of information for debugging.

    run_0D_simulation(mechanism, array_type, Press, Temperature,
                      fuel_name, oxidizer_name, diluent_name, Mixture_type,
                      Mix_params, SpecificSpecies, Starttime, Endtime,
                      Delta_T, PPM, save_files, save_time)
