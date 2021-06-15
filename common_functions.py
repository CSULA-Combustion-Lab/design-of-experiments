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


def normalize_mixture(mix):
    """
    Convert a mixture into a mixture dictionary, normalize mole fractions.

    Parameters
    ----------
    mixture_list : list, str, or dict
        [species1, molefrac1, species2, molefrac2....]
        'species'
        {species1: molefrac1, species2: molefrac2...}

    Returns
    -------
    mixture_dict : dict
        {species1: molefrac1, species2: molefrac2....}

    """
    if type(mix) is str:
        mixture_dict = {mix: 1.0}

    elif type(mix) is list:
        total = sum(mix[1::2])
        mixture_dict = {}
        for i in range(0, len(mix), 2):
            mixture_dict[mix[i]] = mix[i+1] / total

    elif type(mix) is dict:
        mixture_dict = {}
        total = sum(mix.values())
        for k, v in mix.items():
            mixture_dict[k] = v / total

    else:
        msg = ('Fuel, diluent, or oxidizer input format is incorrect.'
               ' It should be a string with the name of the component, a'
               ' list of the format '
               '[component1, quantity1, component2, quantity2...], or a dict '
               '{comp1: quant1, comp2: quant2}')
        raise TypeError(msg)

    return mixture_dict


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
    Fuel         = cond['Mixture'][0]
    Diluent      = cond['Mixture'][1]
    Oxidizer     = cond['Mixture'][2]
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

    # Check that there isn't any overlap between fuel, oxidizer, and diluent
    if len(Fuel.keys() & Diluent.keys()) > 0:
        raise ValueError('Fuel and Diluent cannot contain the same species')
    if len(Fuel.keys() & Oxidizer.keys()) > 0:
        raise ValueError('Fuel and Oxidizer cannot contain the same species')
    if len(Oxidizer.keys() & Diluent.keys()) > 0:
        raise ValueError('Oxidizer and Diluent cannot contain the same species')

    # Create mixtures
    gas = ct.Solution(chem)

    mixlist = []
    mix_loop = it.product(param1, param2)  # Two parameters for looping.

    if mix_type == 'phi_oxi/dil':
        for equiv, ox_to_dil in mix_loop:
            if ox_to_dil > 1:
                continue  # Impossible mixture

            # Mix the oxidizer with diluent
            reduced_ox = {k: v * ox_to_dil for k, v in Oxidizer.items()}
            reduced_dil = {k: v * (1 - ox_to_dil) for k, v in Diluent.items()}
            diluted_ox = {**reduced_ox, **reduced_dil}

            gas.set_equivalence_ratio(equiv, Fuel, diluted_ox)
            mixlist.append(gas.mole_fraction_dict())

    elif mix_type == 'phi_fuel/dil':
        for equiv, f_to_dil in mix_loop:
            if f_to_dil > 1:
                continue  # Impossible mixture

            # Mix the fuel with diluent
            reduced_f = {k: v * f_to_dil for k, v in Fuel.items()}
            reduced_dil = {k: v * (1 - f_to_dil) for k, v in Diluent.items()}
            diluted_f = {**reduced_f, **reduced_dil}

            gas.set_equivalence_ratio(equiv, diluted_f, Oxidizer)
            mixlist.append(gas.mole_fraction_dict())

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
            undil_mixture = {k: v*var_frac/var_total for k, v in initial_mix.items()}
            dil_comp = {k: v * (1 - var_frac / var_total) for k, v in Diluent.items()}
            mixture = {**undil_mixture, **dil_comp}
            mixlist.append(mixture)

    elif mix_type == 'oxi_fuel':
        for oxi_frac, fuel_frac in mix_loop:
            if fuel_frac + oxi_frac > 1:
                continue  # Impossible mixture
            reduced_fuel = {k: v * fuel_frac for k, v in Fuel.items()}
            reduced_ox = {k: v * oxi_frac for k, v in Oxidizer.items()}
            reduced_dil = {k: v * (1 - fuel_frac - oxi_frac) for k, v in Diluent.items()}
            mixture = {**reduced_fuel, **reduced_ox, **reduced_dil}
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
            reduced_dil = {k: v*dil_frac for k, v in Diluent.items()}
            mixture = {**reduced_var1, **reduced_var2, **reduced_dil}
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
            reduced_dil = {k: v*dil_frac for k, v in Diluent.items()}
            mixture = {**reduced_fuel, **reduced_oxi, **reduced_dil}
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
        Simulation case information with the following structure:
        [[Pressure, Temperature, Mixture], ...]
    cond : dict
        A dictionary of the simulation information specific to the type of
        simulation being performed.
    fun : Function
        Name of the function being used per simulation

    Returns
    -------
    outlist : List
        Results of the function are orderd in a list decending in the time
        that each simulation started. Results information depend on the
        function.

    """
    #Find optimal number of cpus to use
    numcases = len(param) #Number of cases to run
    if cpu_count() == 2 or cpu_count() == 1:
        proc = 1 #Less Powerful Computer
    elif numcases > cpu_count():
        #Number of cases to run on each processor, rounded up
        loops = [np.ceil(numcases/proc) for proc in range(1, cpu_count())]
        # First entry in loops with the minumum number. Add one because
        # of 0-based indexing
        # In the past, I added another in case one process is much slower, but
        # This runs the risk of using all cpu's
        proc = loops.index(min(loops))+1
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


def parameters_string(p_type, P, T, mix_params, chem, fuel, oxidizer, diluent):
    """Return string of useful information.

    Parameters
    ----------
    p_type : str
        Problem type
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
              "\nProblem type: " + p_type +
              "\nMechanism: " + chem + "\nFuel: " + str(fuel) +
              "\nOxidizer: " + str(oxidizer) + "\nDiluent: " +
              str(diluent) +
              "\n[Initial, Final, # of points]\n\tTemperature: " +
              format(T) + " [Kelvin]\n\tPressure Range: " + format(P) +
              " [atm]\n\tMixture Parameters:\n" + mixture_text +
              "\n==========================================================")
    return string


# def calculate_a(fuel, mech):
#     """
#     Calculates the stoichiometric ratio for a given mixture of fuel.

#     Parameters
#     ----------
#     fuel : str or list
#         As a string the variable represents a single species of fuel being used.
#         As a list the variable represents multicomponent fuel species
#         followed by the percentage to the total fuel [Component1, % of total, ...]
#     mech : str
#         A .cti mechanism file containing all reaction and species information.

#     Returns
#     -------
#     a : float
#         The stoichiometric ratio of the given fuel mixture

#     """
#     #fuel C(x)H(y)O(z)
#     gas        = ct.Solution(mech)
#     fuel_index = gas.species(gas.species_index(fuel)).composition
#     if 'C' in fuel_index:
#         x = fuel_index['C']
#     else:
#         x = 0
#     if 'H' in fuel_index:
#         y = fuel_index['H']
#     else:
#         y = 0
#     if 'O' in fuel_index:
#         z = fuel_index['O']
#     else:
#         z = 0
#     a = x+y/4-z/2
#     return a

def update_progress(progress):
    """
    A progress bar used in for loops to show haw many iterations have completed
    and how many are left.

    Parameters
    ----------
    progress : int or float
        A value of the current iteration being performed or just finished.

    Returns
    -------
    None.

    """
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
    components : str, list, or dict
        As a string the variable represents a single species of a component.
        As a list the variable represents multiple species in componet .
        As a dictionary species are listed as the key with no given value.
    mix : TYPE
        DESCRIPTION.

    Raises
    ------
    KeyError
        If components are missing in the mixture a 0.0 is returned if the
        components is a string and a pass is used if the components is a list
    TypeError
        Component type is not listed. Function cannot be used with this type.

    Returns
    -------
    Percentage : float
        The sum percentage of requested components in the mixture.

    """
    if type(components) is str:  # Single Component
        try:
            Percentage =  mix[components]
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
    elif type(components) is dict:
        Percentage = sum((mix[key] for key in components))
    else:
        raise TypeError

    return Percentage