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
    """Generate the path given a mechanism file.

    If the Model folder is not already created, one will be generated.
    All mechanism files should be in the model folder.

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
    # Save Path/Parent Directory
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
    """Find and report duplicate reactions in a model.

    Parameters
    ----------
    gas : object
        Cantera generated gas object created using user provided mechanism

    Returns
    -------
    dup_rxns : dict
        Dictionary containing duplicate reactions of the following format:
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
    """Generate mixture and thermodynamic parameters for each simulation.

    Parameters
    ----------
    cond : dict
        Dictionary of the following format:
        for 1D
        conditions = {'Parameters': [Press, Temperature, mix_params, array_type],
                      'Mixture': [fuel, diluent, oxidizer],
                      'Flame': [mingrid, mul_soret, loglevel],
                      'Files': [mechanism, flame_temp],
                      'T/F': [multifuel, multioxidizer]}
        for 0D
        conditions = {'Parameters': [Press, Temperature, mix_params, array_type],
                      'Mixture':[fuel_name, diluent_name, oxidizer_name],
                      'ZeroD': [SpecificSpecies, dup_reactions],
                      'Time_Info': [starttime, endtime, SCORE3_TIME],
                      'Files': [mechanism],
                      'Limits': [delta_T, ppm]}

    Returns
    -------
    paramlist : list
        List that is used to define the initial state for simulations.
        Format is [[Pressure, Temperature, Mixture dictionary], [...], ...]

    """
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    Press        = cond['Parameters'][0]
    Temperature  = cond['Parameters'][1]
    mix_params   = cond['Parameters'][2]
    array_type   = cond['Parameters'][3]
    chem         = cond['Files'][0]

    # Create ranges for the parameters
    def logspace(start, stop, number):
        """Custom, simplified logspace function."""
        return np.logspace(np.log10(start), np.log10(stop), number)

    if array_type == 'log':
        function = logspace
    elif array_type == 'lin':
        function = np.linspace
    else:
        print('Error! Check array_type variable for invalid string input')
        sys.exit()

    # Should we have some parameters, like phi, that are always "lin" no matter what?
    P    = function(*Press)
    T    = function(*Temperature)
    param1 = function(*mix_params[1])
    param2 = function(*mix_params[2])
    mix_type = mix_params[0]

    # Deal with multi-fuel and multi-oxidizer
    msg = ('Fuel or oxidizer input format is incorrect. It should be a string '
           'with the name of the component, or a list of the format '
           '[component1, quantity1, component2, quantity2...]')
    if type(Fuel_name) is str:
        Fuel = {Fuel_name: 1}
    elif type(Fuel_name) is list:  # multi_fuel
        Fuel = {}
        for fl in range(0, len(Fuel_name), 2):
            Fuel[Fuel_name[fl]] = Fuel_name[fl+1]
    else:
        raise ValueError(msg)

    if type(Oxidizer_name) is str:
        Oxidizer = {Oxidizer_name: 1}
    elif type(Oxidizer_name) is list:  # multi_ox
        Oxidizer = {}
        for ox in range(0, len(Oxidizer_name), 2):
            Oxidizer[Oxidizer_name[ox]] = Oxidizer_name[ox+1]
    else:
        raise ValueError(msg)

    # Create mixtures
    gas = ct.Solution(chem)

    mixlist = []
    mix_loop = it.product(param1, param2)  # Two parameters for looping.

    # TODO: The quantity that the code is referring to as "X to diluent ratio"
    # is very confusing.
    # Initializer calls it "Diluent_Percentage" which I would think of as
    # "diluent / (X+diluent).
    # The rest of the code (and outputs) call it "X to diluent ratio", which
    # I would understand to mean "X/diluent".
    # But, the quantity is actually "X/(X+diluent)", which I would call
    # "X% in X mixture," like "O2 % in oxidizer." I've started fixing this.
    if mix_type == 'phi_oxi/dil':
        for equiv, ox_to_dil in mix_loop:
            if ox_to_dil > 1:
                continue  # Impossible mixture

            # Mix the oxidizer with diluent
            diluted_ox = {k: v * ox_to_dil for k, v in Oxidizer.items()}
            diluted_ox[Diluent_name] = 1 - ox_to_dil

            gas.set_equivalence_ratio(equiv, Fuel, diluted_ox)
            mixlist.append(gas.mole_fraction_dict())

    elif mix_type == 'phi_fuel/dil':
        for equiv, f_to_dil in mix_loop:
            if f_to_dil > 1:
                continue  # Impossible mixture

            # Mix the fuel with diluent
            diluted_f = {k: v * f_to_dil for k, v in Fuel.items()}
            diluted_f[Diluent_name] = 1 - f_to_dil

            gas.set_equivalence_ratio(equiv, diluted_f, Oxidizer)
            mixlist.append(gas.mole_fraction_dict())

    # elif mix_type == 'phi_fuel':
    #     for equiv, fuel_frac in mix_loop:
    #         if fuel_frac > 1:
    #             continue  # Impossible mixture

    #         gas.set_equivalence_ratio(equiv, Fuel, Oxidizer)  # Without diluent
    #         initial_mix = gas.mole_fraction_dict()
    #         fuel_total = sum([initial_mix[k] for k in Fuel])
    #         if fuel_total < fuel_frac:
    #             continue  # Cannot create mixture at this phi + fuel
    #         mixture = {k: v*fuel_frac/fuel_total for k, v in initial_mix.items()}
    #         mixture[Diluent_name] = 1 - fuel_frac / fuel_total
    #         mixlist.append(mixture)
            
    elif mix_type in('phi_fuel', 'phi_oxi'):
        for equiv, var_frac in mix_loop:
            if var_frac > 1:
                continue  # Impossible mixture

            gas.set_equivalence_ratio(equiv, Fuel, Oxidizer)  # Without diluent
            initial_mix = gas.mole_fraction_dict()
            if mix_type == 'phi_fuel':
                Variable = Fuel
            elif mix_type == 'phi_oxi':
                Variable = Oxidizer
            var_total = sum([initial_mix[k] for k in Variable])
            if var_total < var_frac:
                continue  # Cannot create mixture at this phi + fuel or oxidizer
            mixture = {k: v*var_frac/var_total for k, v in initial_mix.items()}
            mixture[Diluent_name] = 1 - var_frac / var_total
            mixlist.append(mixture)

    elif mix_type == 'oxi_fuel':
        for oxi_frac, fuel_frac in mix_loop:
            if fuel_frac + oxi_frac > 1:
                continue  # Impossible mixture
            reduced_fuel = {k: v * fuel_frac for k, v in Fuel.items()}
            reduced_ox = {k: v * oxi_frac for k, v in Oxidizer.items()}
            mixture = {**reduced_fuel, **reduced_ox,
                       Diluent_name: 1 - fuel_frac - oxi_frac}
            mixlist.append(mixture)
            
    elif mix_type in ('fuel_dil', 'oxi_dil'):
        for var1_frac, dil_frac in mix_loop:
            if var1_frac + dil_frac > 1:
                continue # Impossible mixture
            elif mix_type == 'fuel_dil':
                Variable1 = Fuel
                Variable2 = Oxidizer
            elif mix_type == 'oxi_dil':
                Variable1 = Oxidizer
                Variable2 = Fuel
            var2_frac = 1 - var1_frac - dil_frac
            reduced_var1 = {k: v*var1_frac for k, v in Variable1.items()}
            reduced_var2 = {k: v*var2_frac for k, v in Variable2.items()}
            mixture = {**reduced_var1, **reduced_var2, Diluent_name: dil_frac}
            mixlist.append(mixture)
            
    elif mix_type == 'phi_dil':
        for equiv, dil_frac in mix_loop:
            if dil_frac > 1:
                continue # Impossible mixture
            gas.set_equivalence_ratio(equiv, Fuel, Oxidizer)  # Without diluent
            initial_mix = gas.mole_fraction_dict()
            fueloxi_frac = 1 - dil_frac
            fuel_frac = sum([initial_mix[k] for k in Fuel])
            oxi_frac  = sum([initial_mix[k] for k in Oxidizer])
            reduced_fuel = {k: v*fueloxi_frac*fuel_frac for k, v in Fuel.items()}
            reduced_oxi = {k: v*fueloxi_frac*oxi_frac for k, v in Oxidizer.items()}
            mixture = {**reduced_fuel, **reduced_oxi, Diluent_name: dil_frac}
            mixlist.append(mixture)

    else:
        print('Error creating mixtures. Check mixture_type variable.')
        sys.exit()

    return list(it.product(P, T, mixlist))


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


def parameters_string(P, T, mix_params, chem, fuel, oxidizer, diluent):
    """Return string of useful information.

    Parameters
    ----------
    P : list
        [initial, final, # of points]
    T : list
        [initial, final, # of points]
    mix_params : tuple
        (type, param1, param2)
    chem : str
        Chemistry file
    fuel, oxidizer, diluent : dict, str, or list

    Returns
    -------
    string to be printed

    """
    mt = mix_params[0]
    labels = ['Mixture Type']
    if mt == 'Oxi_Dil':
        labels.extend(['Equivalence Ratio', 'O2 fraction in oxidizer'])
    elif mt == 'Fue_Dil':
        labels.extend(['Equivalence Ratio', 'Fuel fraction in fuel mix'])
    elif mt == 'phi_fuel':
        labels.extend(['Equivalence Ratio', 'Fuel fraction in mixture'])
    elif mt == 'Ox_Fuel':
        labels.extend(['Oxidizer mole fraction', 'Fuel mole fraction'])
    else:
        labels.extend(['Parameter 1', 'Parameter 2'])
    mixture_text = '\n'.join(['\t\t{}: {}'.format(k, v) for
                              k, v in zip(labels, mix_params)])
    string = ("========================Parameters========================" +
              "\nMechanism: " + chem + "\nFuel: " + str(fuel) +
              "\nOxidizer: " + str(oxidizer) + "\nDiluent: " +
              str(diluent) +
              "\n[Initial, Final, # of points]\n\tTemperature: " +
              format(T) + " [Kelvin]\n\tPressure Range: " + format(P) +
              " [atm]\n\tMixture Parameters:\n" + mixture_text +
              "\n==========================================================")
    return string


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
    gas        = ct.Solution(mech)
    fuel_index = gas.species(gas.species_index(fuel)).composition
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
    
def mixture_percentage(components, mix):
    """


    Parameters
    ----------
    components : TYPE
        DESCRIPTION.
    mix : TYPE
        DESCRIPTION.

    Raises
    ------
    TypeError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if type(components) is str:  # Single Component
        try:
            return mix[components]
        except KeyError:
            return 0.0
    elif type(components) is list:
        # Format is [component1, quantity1, component2, quantity2, ...]
        Percentage = 0
        for n in range(0, len(components), 2):
            try:
                Percentage += mix[components[n]]
            except KeyError:  # Component isn't in mixture
                pass
        return Percentage
    else:
        raise TypeError