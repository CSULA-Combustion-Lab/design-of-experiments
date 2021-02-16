#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 13:29:22 2018

@author: boss
"""

""" This is just a change to the discription to test
pushing a change in the code.
 """

import cantera
import numpy as np
import math
import time, sys
import itertools as it
import pickle
import os
# import copy
from datetime import datetime
import common_functions as cf

cantera.suppress_thermo_warnings()


def reac_sens_parallel(temp, pressure, phi, fuel, pack):
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

    starttime       = pack['Time'][0]
    endtime         = pack['Time'][1]
    fuel_name       = pack['Names'][0]
    diluent_name    = pack['Names'][1]
    SpecificSpecies = pack['Names'][2]
    mix_type        = pack['Mixture_Info'][0]
    mechanism       = pack['Mixture_Info'][1]
    dup_reactions   = pack['Mixture_Info'][2]
    starttime       = pack['Time_Info'][0]
    endtime         = pack['Time_Info'][1]
    SCORE3_TIME     = pack['Time_Info'][2]
    a               = pack['Parameters'][5]


    # oxygen   = a/phi*fuel
    # diluent  = 1 - oxygen - fuel
    # mix      = {fuel_name:fuel, 'O2':oxygen, diluent_name:diluent}
    gas1     = cantera.Solution(mechanism)
    mix      = mixture_maker(gas1, phi, fuel, a,
                             mix_type, fuel_name, diluent_name)
    gas1.TPX = temp, pressure*101325, mix
    rxns     = range(0,len(gas1.reaction_equations()))
    #Insert for test#
    # dup_reactions = duplicate_reactions(gas1)
    #End Insert#

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
        row[4]   = gas1[SpecificSpecies].X    #save molesfrac value of specific species
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
               All_tTP_AMo_test, temp, mix, pressure, phi, fuel]

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

    max_temp = t_T_P_AMol[-1][1] #temp of the last time step for condition
    if max_temp-temp > delta_T:
        temp_condition = 1
    else:
        temp_condition = 0
    return temp_condition


def mole_fractions():
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

    molfrac_conditions      = [None]*(len(SpecificSpecies)+2)
    molfrac_conditions[0:2] = [mix, temp, pressure]

    for i in range(0,len(SpecificSpecies)):
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


def specific_sens():
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
    row = [None]*len(SpecificSpecies)*len(rxns)
    molefrac_time = np.array([x[0] for x in t_SMol])
    for i in range(0,len(molefrac_time)):
        for k in range(0,len(rxns)):
            for j in range(0,len(SpecificSpecieNumbers)):
                molfractions = molfrac_conditions[3+j]
                MolFractions = molfractions[:,0]
#                MolFractionsTime = molfractions[:,1]
                if MolFractions[i] == 0:
                    row[j+k*len(SpecificSpecies)] = 0
                elif MolFractions[i] != 0:
                    row[j+k*len(SpecificSpecies)] = sensitivities[i,SpecificSpecieNumbers[j],k]
        SpecificSpecieSens.append([x for x in row])
    return SpecificSpecieSens


def sensitivity_score():
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
    row   = [None]*(len(SpecificSpecies)+2)
    row1  = [None]*(len(SpecificSpecies))

    row[0] = temp             #K
    row[1] = pressure*101.325 #kPa

    rxn_maxs = np.zeros(len(rxns))
    rxn_t    = np.zeros(len(rxns))

    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(rxns)):
            rxn_maxs[j] = max(dataa[:,i+j*len(SpecificSpecies)])
            rxn_t[j]=senstime[np.argmax(dataa[:,i+j*len(SpecificSpecies)])]
        row[i+2] = sum(rxn_maxs)/len(rxns) #avg of max sensitivity per reactions
        row1[i]  = sum(rxn_t)/len(rxns)
    score_T_P_MaxSensavg.append([x for x in row]) #T,P,avg between maximum sens per reaction
    scoretimes.append([x for x in row1]) #what does this populate?


def sensitivity_score2():
    """ Find the reaction with the highest absolute value of sensitivity for each species of interest.
    The maximum sensitivity is found for all time - not at a specific time step.
    This score will find spikes, but may not find a reaction that is sensitive for a longer time.
    """
    dataa2  = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    ss      = [None]*(len(SpecificSpecies)+2)  #ss-sensititvity score
    ss_time = [None]*(len(SpecificSpecies))
    max_rxn = [None]*(len(SpecificSpecies)+3)

    ss[0] = temp
    ss[1] = pressure*101.325 #in kPa for the plot
    max_rxn[0] = mix
    max_rxn[1] = temp
    max_rxn[2] = pressure*101.325
    rxn_name   = [None]*(len(SpecificSpecies)+1)
    rxn_name[0] = {'temperature': temp, 'pressure': pressure*101.325,
                   'equivalence': phi, 'fuel': fuel}

    max_sens = np.zeros(len(rxns))
    rxn_t    = np.zeros(len(rxns))

    for i in range(0,len(SpecificSpecies)):
        for j in range(0,len(rxns)):
            max_sens[j] = max(dataa2[:,i+j*len(SpecificSpecies)])
            rxn_t[j] = senstime[np.argmax(dataa2[:,i+j*len(SpecificSpecies)])]
        ss[i+2]      = max(max_sens)              #maximum sensitivity
        ss_time[i]   = rxn_t[np.argmax(max_sens)] #time of max sensitivity
        max_rxn[i+3] = rxns[np.argmax(max_sens)]  #rxn with max sensitivity
        rxn_name[i+1] = gas.reaction_equations([rxns[np.argmax(max_sens)]])

    score2_T_P_MaxSens.append([x for x in ss])
    score2times.append([x for x in ss_time])
    score2_Max_sens_rxn.append([x for x in max_rxn])
    score2_Params_MaxSens_Name_Params.append([x for x in rxn_name])


def sensitivity_score3():
    num_spec = len(SpecificSpecies)
    max_rxn = [None]*(num_spec+3)
    max_rxn[0] = mix
    max_rxn[1] = temp
    max_rxn[2] = pressure*101.325
    rxn_name   = [None]*(num_spec+1)
    rxn_name[0] = {'temperature': temp, 'pressure': pressure*101.325,
                   'equivalence': phi, 'fuel': fuel, 'score3time': SCORE3_TIME}

    time = np.array([x[0] for x in t_AllSpecieSens])
    ind = np.argwhere(time > SCORE3_TIME)[0][0]
    sens_at_t = SpecificSpecieSens[ind]

    for i in range(num_spec):
        unshuffled = sens_at_t[i::num_spec]
        max_sens_ind = np.argmax(unshuffled)
        if unshuffled[max_sens_ind] == 0:  # Below ppm limit
            max_rxn[i+3] = None
            rxn_name[i+1] = None
        else:
            max_rxn[i+3] = rxns[max_sens_ind]  #rxn with max sensitivity
            rxn_name[i+1] = gas.reaction_equations([rxns[max_sens_ind]])

    score3_Max_sens_rxn.append([x for x in max_rxn])
    score3_Params_MaxSens_Name_Params.append([x for x in rxn_name])


def rank_all(SpecificSpecieSens):
    """ Creates a list of reaction rankings for all timesteps and all species.
        Needs work. Understanding output rank_mat: The lists cycle though the
        species list, and after every len(species) its a new time_step."""
    data3 = np.absolute(np.array([x[:] for x in SpecificSpecieSens]))
    all_ranks_params    = [None]*1
    all_ranks_params[0] = {'temperature': temp, 'pressure': pressure*101.325,
                         'equivalence': phi, 'fuel': fuel}
    step  = [None]*len(rxns)
    for i in range(0,len(senstime)):
        for j in range(0,len(SpecificSpecies)):
            for k in range(0,len(rxns)):
                step[k] = data3[i,j+k*len(SpecificSpecies)] #creates list of specie sensitivities per rxns, per time, per specie
            ranking = ranking_per_step(step)
            all_ranks.append([x for x in ranking])
            all_ranks_params.append([x for x in ranking])

    All_Ranks.append([x for x in all_ranks])
    All_Ranks_Params.append([x for x in all_ranks_params])


def ranking_per_step(sval):
    """ Assigns a rank to each reaction sensitivity per time-step and species
    """
    sval  = np.array(sval)   #sval = step from above per time, per specie
    order = np.argsort(sval) #orders sval from smallest to largeest by listing indices which sval elements should be arranged in
    ranks = ranks_possible(sval)
    step_rank = [None]*len(sval)
    for j in range (0,len(sval)):
        step_rank[order[j]] = ranks[j]
    return step_rank


def ranks_possible(sval):
    """ Determines how many ranks there will be for a given time-step, since
    sensitivities can tie if they are equal"""
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
        IS_norm.append(spec / np.mean(np.abs(top5)))

    return IS_norm


def update_progress(progress):
    barLength = 25 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block +
                       "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def flow_reactor_parallel(params, pack, fun):
    """Runs the simulation in Parallel by first determining the largest length
    of the four parameters. If all parameters are of equal length, the code will
    create a pool of temperatures to run. The code must be outside of ipython
    console, if not it will not parallize the simulation"""
    from multiprocessing import cpu_count, Pool
    #Find optimal number of cpus to use
    numcases = len(params) #Number of cases to run
    if cpu_count() == 2 or cpu_count() == 1:
        proc = 1 #Less Powerful Computer
    elif numcases > cpu_count():
        #Number of cases to run on each processor, rounded up
        loops = [np.ceil(numcases/proc) for proc in range(1, cpu_count())]
        # First entry in loops with the minumum number. Add one because
        # of 0-based indexing, add another in case one process is much slower.
        proc = loops.index(min(loops))+2
    else: # More cpus than cases
        proc = numcases

    pool = Pool(processes=proc)

    results = []
    for x in params:
        results.append(pool.apply_async(fun, args=(*x, pack)))
    pool.close()
    pool.join()

    # Get the results
    datadict = dict()
    casenum  = 0
    for p in results:
        try:
            # Assign it to datadict. This is ordered by the time when each
            # simulation starts, not when they end
            datadict[casenum] = p.get()
        except RuntimeError: # I'm not sure what this is
            print('\nUnknown Cody RunTimeError.')
            datadict[casenum] = None
        casenum += 1
    outlist = [datadict[k] for k in datadict.keys()] # Convert back to list
    return outlist


def reduce_cases(paramlist, a, debug=False):
    """ Reduces cases if the diluent is negative."""
    paramlist_condition = []
    parameters_list     = []
    for x in paramlist:
        (temp, pressure, phi, fuel) = x
        oxygen  = a/phi*fuel
        diluent = 1 - oxygen - fuel
        if diluent<0:
            paramlist_condition.append(1)
            if debug:
                print('Diluent limit. ' + DEBUG_FMT.format(*x))
        else:
            paramlist_condition.append(0)

    Cases = range(0,len(paramlist))
    for x in Cases:
        if paramlist_condition[x] == 0:
            parameters_list.append(paramlist[x])
    return parameters_list


def sens_dup_filter(sens_information, duplicate_reactions):
    """
    

    Parameters
    ----------
    sens_information : TYPE
        DESCRIPTION.
    duplicate_reactions : TYPE
        DESCRIPTION.

    Returns
    -------
    sens_information : TYPE
        DESCRIPTION.

    """
    for s in sens_information:
        for d in duplicate_reactions.keys():
                    sumt = 0
                    for n in duplicate_reactions[d]:
                         sumt += s[n]
                         s[n] = 0
                    s[duplicate_reactions[d][0]] = sumt
    return sens_information


def calculate_a(fuel, mech):
    """
    

    Parameters
    ----------
    fuel : TYPE
        DESCRIPTION.
    mech : TYPE
        DESCRIPTION.

    Returns
    -------
    a : TYPE
        DESCRIPTION.

    """
    #fuel C(x)H(y)O(z)
    gas        = cantera.Solution(mech)
    fuel_index = gas.species(gas.species_index(fuel_name)).composition
    if 'C' in fuel_index:
        x = fuel_index['C']
    else:
        x = 0
    if 'H' in fuel_index:
        y = fuel_index['H']
    else:
        y = 0
    if 'O' in fuel_index:
        z = fuel_index['O']
    else:
        z = 0
    a = x+y/4-z/2
    return a


def mixture_maker(gas, phi, fuel, a, mtype, fname, dname):
    if mtype == "Fuel_Phi":
        oxygen   = a/phi*fuel
        diluent  = 1 - oxygen - fuel
        mixture  = {fname:fuel, 'O2':oxygen, dname:diluent}
    elif mtype == "DiluentMix_Fuel":
        fstring = fname+":"+str(fuel)+","+dname+":"+str(1-fuel)
        gas.set_equivalence_ratio(phi, (fstring), ('O2'))
        mixture = gas.mole_fraction_dict()
    elif mtype == "DiluentMix_Oxi":
        fstring = 'O2'+":"+str(fuel)+","+dname+":"+str(1-fuel)
        gas.set_equivalence_ratio(phi, (fname), (fstring))
        mixture = gas.mole_fraction_dict()
    return mixture


if __name__ == "__main__":
    ####Set experiment parameters
    mech_name = 'mech-FFCM1_modified.cti' #Mechanism file
    # mechanism = cf.model_folder(mech_name)
    mechanism = os.path.join('Models',mech_name)
    #Parameters for main loop
    T    = np.linspace(600, 2500, 2)                        #Temperature [K]
    P    = np.logspace(np.log10(.1), np.log10(100), 2)      #Pressure [atm]
    Phi  = np.logspace(np.log10(0.00025), np.log10(2.5), 2) #Equivalence ratio
    Fuel = np.logspace(np.log10(0.00001), np.log10(0.1), 2) #Fuel mole fraction
    Dilper = np.logspace(np.log10(0.00001), np.log10(0.1), 2)

    #Parameters for mixture
    fuel_name = 'H2' #chemical formula of fuel
    diluent_name = 'N2' #chemical formula of diluent
    mixture_type = 'Custom'

    SpecificSpecies = ['OH'] #Species of interest for rxn ranking data
    starttime = 0     #in the case a reading is too early
    endtime   = 0.001 #one milisecond
    delta_T   = 100
    ppm       = 1/1000000 #one ppm
    SCORE3_TIME = starttime + (endtime - starttime)*0.75 # time to evaluate score3

    save_files = True # If true, save files for plotting script
    save_time  = True # If true, also save files for GUI. Not recommended for large runs
    debug = False  # If true, print lots of information for debugging.

    DEBUG_FMT = 'Removing condition: T={:.0f}, P={:.0f}, phi={:.3g}, fuel={:.3g}'
    #start time
    tic = time.time()
    #initializing
    iteration       = 0
    totaliterations = len(T)*len(P)*len(Fuel)*len(Phi)
    paramlist       = list(it.product(T,P,Phi,Fuel))
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
    #create gas object
    gas                   = cantera.Solution(mechanism)
    a                     = calculate_a(fuel_name, mechanism)
    dup_reactions         = cf.duplicate_reactions(gas)
    rxns                  = range(0,len(gas.reaction_equations()))
    names                 = gas.species_names
    SpecificSpecieNumbers = [names.index(speci) for speci in SpecificSpecies]

    #Package key parameters for simulations
    Packed = {'Time': [starttime, endtime],
              'Names': [fuel_name, diluent_name, SpecificSpecies],
              'Parameters': [T, P, Phi, Fuel, Dilper, a],
              'Mixture_Info':[mixture_type, mechanism, dup_reactions],
              'Time_Info': [starttime, endtime, SCORE3_TIME],
              'Limits': [delta_T, ppm]}

    print('Initial Number of Cases: '+format(len(paramlist)))
    sim_start   = time.time()
    if mixture_type == 'Custom':
        params = reduce_cases(paramlist, a, debug)
    else:
        params = paramlist
    print(format(len(params))+' cases have a dilutent greater than -0')
    print('Start of Simulations')
    # sim_package = cf.parallelize(params, Packed, reac_sens_parallel)
    sim_package = flow_reactor_parallel(params, Packed, reac_sens_parallel)
    sim_end     = time.time()
    sim_time    = sim_end - sim_start
    # sim_package_unfiltered = copy.deepcopy(sim_package_uf)
    # sim_package            = sim_info_filter(sim_package_uf, dup_reactions)
    print('End of Simulations')
    print('It took '+format(sim_time, '0.5f')+' seconds to simulate all cases')

    Cases = len(params)
    case  = 0
    print('Start Ranking Process')
    rank_start = time.time()
    update_progress(0)
    for sims in sim_package:
        [t_T_P_AMol, t_SMol, t_AllSpecieSens,
         All_time_Sens_test, All_tTP_AMo_test,
         temp, mix, pressure, phi, fuel] = sims
        temp_condition = T_limit(t_T_P_AMol, temp, delta_T)
        if temp_condition == 1:
            if debug:
                print('T limit ' + DEBUG_FMT.format(temp, pressure, phi, fuel))
            del All_tTP_AMo_test[-1]
            del All_time_Sens_test[-1]
            case +=1
            update_progress(case/Cases)
            continue
        else:
            All_tTP_AMo.append(All_tTP_AMo_test)
            All_time_Sens.append(All_time_Sens_test)
            condition_info = {
                'temperature': temp, 'pressure': pressure*101.325,
                'equivalence': phi, 'fuel': fuel}
            Parameters.append(condition_info)

        #initialize parameters
        SpecificSpecieSens = []   #list of sens of interest per rxn
        all_ranks          = []
        senstime           = np.array([x[0] for x in t_AllSpecieSens])   #times for all specie sensitivities
        sensitivities      = np.absolute(np.array([x[1] for x in t_AllSpecieSens])) #sensitivities
        molfrac_conditions = mole_fractions()
        specific_sens()
        sensitivity_score()
        sensitivity_score2()
        sensitivity_score3()
        rank_all(SpecificSpecieSens)
        IS_norm = integrated(t_AllSpecieSens, SpecificSpecieNumbers, len(rxns),
                             molfrac_conditions)
        int_strength.append([condition_info, IS_norm])
        case += 1
        update_progress(case/Cases)

    rank_end  = time.time()
    rank_time = rank_end - rank_start
    print('End Ranking Process')
    print('It took '+format(rank_time, '0.5f')+' seconds to rank all cases')

    toc       = time.time()
    duration  = toc-tic
    numbr_mix = len(Phi)
    tp_combs  = len(T)*len(P)
    print('Number of cases ran', len(paramlist))
    print('Number of cases that meet all conditions', len(All_Ranks))
    print(numbr_mix, 'mixture with',tp_combs,
          'Temperature/Pressure combinations took '+format(duration, '0.5f')+
          ' seconds.\n')

    if save_files:
        print('Creating Directory...')
        #Creat Directory Name
        now = datetime.now()
        dt_string = now.strftime("%d_%m_%Y %H.%M.%S")
        # dt_string = 'test2'
        directory = dt_string

        #Save Path
        parent_dir = r'Outputs'
        save_path = os.path.join(parent_dir, directory)
        os.makedirs(save_path)

        figures_dir = 'Figures'
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
                "\r"+format(SpecificSpecies)+"\r"
                "=================================\r\n"
                "=====================Parameters==========================\r"
                "Temperatures 	    = "+format(T)+" [K]\r"
                "Pressures    	    = "+format(P)+" [atm]\r"
                "Equivalence Ratios  = "+format(Phi)+" [phi]\r"
                "Fuel Mole Fractions = "+format(Fuel)+" [mole fraction]\r"
                "=========================================================\r\n"
                "==========Number of Cases==========\r"
                "Total Cases Simulated 	     = "+format(len(paramlist))+"\r"
                "Cases Meeting All Conditions = "+format(len(All_Ranks))+"\r"
                "===================================\r\n"
                "=====Conditions=====\r"
                "PPM Limit     = "+format(ppm*1e6)+"\r"
                "delta_T Limit = "+format(delta_T)+"\r"
                "====================\r\n"
                "====================Times====================\r"
                "Simulation Time = "+format(sim_time, '0.2f')+" [seconds]\r"
                "Ranking Time    = "+format(rank_time, '0.2f')+" [seconds]\r"
                "Total Time      = "+format(duration, '0.2f')+" [seconds]\r"
                "=============================================")
        f.write(text_description)
        f.close()

        print('\nText File Created\nSaving Files...\n')
        #File saving
        gas_reactions = cantera.Solution(mechanism).reaction_equations()
        species_rxn_file = [SpecificSpecies, [len(rxns)], [a],
                                             [gas_reactions], [SCORE3_TIME]]

        # This loop can make the saving easier. Just put the thing to be
        # saved and the file name in the save_list
        save_list = [(int_strength, 'Integrated Strength.pkl'),
                     (All_Ranks_Params, 'All_Ranks.pkl'),
                     (score3_Max_sens_rxn, 'Max_Sens_Rxn.pkl'),
                     (species_rxn_file, 'Species_Rxn.pkl'),
                     (Parameters, 'Case_Parameters.pkl')]
        for item in save_list:
            file = os.path.join(save_path, item[1])
            with open(file, 'wb') as f:
                pickle.dump(item[0], f)

        if save_time:
            Mole_Frac_fname = 'Mole_fractions.pkl'
            file = os.path.join(save_path, Mole_Frac_fname)
            with open(file, "wb") as f:
                pickle.dump(All_tTP_AMo, f)
            f.close()

            file = os.path.join(save_path, 'All_time_sens.pkl')
            with open(file, "wb") as f:
                pickle.dump(All_time_Sens, f)
            f.close()
            print('Files Saved')

        else:
            print('File Saved')
