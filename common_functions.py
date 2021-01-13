# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 14:15:27 2021

@author: Kodo Bear
"""

import os
import sys
import errno
import numpy as np
import cantera as ct
import itertools as it
from tqdm import tqdm
from multiprocessing import cpu_count, Pool


def model_folder(mechanism):
    """
    Generates the path given a mechanism file. If the Model folder is not
    all ready created one will be generated. All mechanism files should be in
    the model folder.

    Parameters
    ----------
    mechanism : .cti
        Mechanism file with information on species, reactions, thermodynamic
        properties of a fuel.

    Returns
    -------
    model_path : Path
        Path for the mechanism file found inside the Models folder.

    """
    #Save Path/Parent Directory
    parent_dir = 'Models'
    try:
        os.makedirs(parent_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    model_path = os.path.join(parent_dir, mechanism)
    if not os.path.exists(model_path):
        print("Mechanism file not found in Models folder.")
        sys.exit()

    return model_path


def duplicate_reactions(gas):
    """
    Generates a dictionary of reaction equations and
    their reaction numbers for a given mechanism.

    Parameters
    ----------
    gas : object
        Cantera generated gas object created using user provided mechanism

    Returns
    -------
    dup_rxns : dict
        Dictionary of the following format:
        dup_rxns = {'Reaction Equation 1': [Rxn Number_a, Rxn Number_b]
                            :
                            :
                    'Reaction Equation N': [Rxn Number_a, Rxn Number_b]}
    """
    dup_rxns = {}
    eqns = gas.reaction_equations()
    for i in range(len(eqns)):
        if gas.reaction(i).duplicate:
            rxn_eqn = gas.reaction_equation(i)
            if rxn_eqn not in dup_rxns.keys():
                dup_rxns[rxn_eqn] = [i]
            elif rxn_eqn in dup_rxns.keys():
                dup_rxns[rxn_eqn].append(i)
            else:
                print('Something went wrong!')
    return dup_rxns


def case_maker(cond):
    """
    Generate parameters for each simulation.

    Parameters
    ----------
    cond : dict
        Dictionary of the following format:
        for 1D
        conditions = {'Parameters': [P, Phi, FtD, Tint, OtD, array_type],
                      'Mixture': [fuel, diluent, oxidizer, mixture_type],
                      'Flame': [mingrid, mul_soret, loglevel],
                      'Files': [mechanism, flame_temp],
                      'T/F': [multifuel, multioxidizer]}
        for 0D
        conditions =
    Returns
    -------
    paramlist : list
        List that is used to define the initial state for simulations.
        Format is [[Pressure, Temperature, Mixture dictionary], [...], ...]

    """
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    mix_type = cond['Mixture'][3]
    p        = cond['Parameters'][0]
    phi      = cond['Parameters'][1]
    T        = cond['Parameters'][3]
    # fod      = cond['Parameters'][2]
    otd      = cond['Parameters'][4]
    ftd      = cond['Parameters'][2]
    chem     = cond['Files'][0]

    gas = ct.Solution(chem)

    msg = ('Fuel or oxidizer input format is incorrect. It should be a string '
           'with the name of the component, or a list of the format '
           '[component1, quantity1, component2, quantity2...]')
    paramlist = []

    if mix_type == 'Debug':
        print('Debug Loop Enabled')
    elif mix_type == 'Custom':
        print('Custom Loop Enabled')
        #Under Construction
        for i in p:
            for k in otd:
                paramlist.append((i, ftd, k))
        #Under Construction

    # TODO: The quantity that the code is referring to as "X to diluent ratio"
    # is actually "X/(X+diluent)". This is a different quantity, with a
    # different meaning. "X to diluent ratio" would be "X/diluent". We should
    # either change the definition in this function, or change the labels on
    # many figures and other outputs.
    elif mix_type == 'Oxi_Dil':
        print('Oxidizer to Diluent Loop Enabled')
        P_T_phi_otd = it.product(p, T, phi, otd)

        for pressure, temperature, equiv, ox_to_dil in P_T_phi_otd:
            if ox_to_dil > 1:
                continue  # Impossible mixture

            if type(Fuel_name) is str:
                Fuel = Fuel_name
            elif type(Fuel_name) is list:  # multi_fuel
                Fuel = ''
                for fl in range(0,len(Fuel_name),2):
                    Fuel += Fuel_name[fl]+':'+str(Fuel_name[fl+1])+' '
            else:
                raise ValueError(msg)

            if type(Oxidizer_name) is str:
                Oxidizer = Oxidizer_name + ':'+str(ox_to_dil) + ' '
            elif type(Oxidizer_name) is list:  # multi_ox
                Oxidizer = ''
                for ox in range(0, len(Oxidizer_name), 2):
                    Oxidizer += Oxidizer_name[ox] + ':' + str(Oxidizer_name[ox+1]*ox_to_dil) + ' '
            else:
                raise ValueError(msg)
            Oxidizer += Diluent_name + ':' + str(1 - ox_to_dil)

            gas.set_equivalence_ratio(equiv, (Fuel), (Oxidizer))
            paramlist.append([pressure, temperature, gas.mole_fraction_dict()])

    elif mix_type == 'Fue_Dil':
        print('Fuel to Diluent Loop Enabled')

        P_T_phi_ftd = it.product(p, T, phi, ftd)

        for pressure, temperature, equiv, f_to_dil in P_T_phi_ftd:
            if f_to_dil > 1:
                continue  # Impossible mixture

            if type(Fuel_name) is str:
                Fuel = Fuel_name + ':' + str(f_to_dil) + ' '
            elif type(Fuel_name) is list:  # multi_fuel
                Fuel = ''
                for fl in range(0,len(Fuel_name),2):
                    Fuel += Fuel_name[fl]+':'+str(Fuel_name[fl+1]*f_to_dil)+' '
            else:
                raise ValueError(msg)
            Fuel += Diluent_name + ':' + str(1 - f_to_dil)

            if type(Oxidizer_name) is str:
                Oxidizer = Oxidizer_name
            elif type(Oxidizer_name) is list:  # multi_ox
                Oxidizer = ''
                for ox in range(0, len(Oxidizer_name), 2):
                    Oxidizer += Oxidizer_name[ox] + ':' + str(Oxidizer_name[ox+1])
            else:
                raise ValueError(msg)

            gas.set_equivalence_ratio(equiv, (Fuel), (Oxidizer))
            paramlist.append([pressure, temperature, gas.mole_fraction_dict()])


    # TODO: Following this example, add other mixture types.
    # Rewrite 0D conditions dictionary so that is it formated similar to 1D


    # elif mix_type == 'Oxi_Dil' or mix_type == 'Fue_Dil':
    #     if mix_type == 'Oxi_Dil':
    #         print('Oxidizer to Diluent Loop Enabled')
    #     elif mix_type == 'Fue_Dil':
    #         print('Fuel to Diluent Loop Enabled')
    #     totaliterations = len(p)*len(phi)*len(fod)
    #     paramlist       = list(it.product(p,phi,fod))
    else:
        print('Error creating mixtures. Check mixture_type variable.')
        sys.exit()

    return paramlist


def parallelize(param, cond, fun):
    """
    Parrallelize all cases found in param using information found from cond
    and calculated using function defined by fun.

    Parameters
    ----------
    param : List
        DESCRIPTION.
    cond : Dictionry
        DESCRIPTION.
    fun : Function
        Name of the function being used per simulation

    Returns
    -------
    outlist : Dictionary
        DESCRIPTION.

    """
    #Find optimal number of cpus to use
    numcases = len(param) #Number of cases to run
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
    for x in param:
        results.append(pool.apply_async(fun, args=(*x, cond)))
    pool.close()
    # pool.join()  # This waits until all all processes have finished. It's commented so that the progress bar works.

    # Get the results
    datadict = dict()
    casenum  = 0
    for p in tqdm(results):
        try:
            # Assign it to datadict. This is ordered by the time when each
            # simulation starts, not when they end
            datadict[casenum] = p.get()
        except RuntimeError: # I'm not sure what this is
            print('\nUnknown RunTimeError.')
            datadict[casenum] = {'Error': [None, 'RunTimeError']}
        casenum += 1
    outlist = [datadict[k] for k in datadict.keys()] # Convert back to list
    return outlist