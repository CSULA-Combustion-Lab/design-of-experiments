# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:32:37 2020

@author: Kodo Bear
"""
import Sensitized_0D_Experiment as zeroD
import Sensitized_Flame_Experiment as oneD
import common_functions as cf
import Sensitized_Flame_Plotting as flame_plot
import Sensitized_0D_Plotting as zeroD_plot
# import BurnerSimulation
import sys
import yaml


def package(opts):
    """
    Repackage mixture options

    Parameters
    ----------
    mixture_opts : dict
        Mixture options from the input YAML file

    Returns
    -------
    mix_params : list
        (Mixture type, variable1, variable2)

    """
    Mixture_type = opts['Mixture_type']

    # Package the mixture parameters. Each mixture type requires two variables
    if Mixture_type in ('phi_oxi/dil', 'phi_fuel/dil'):
        mix_params = (Mixture_type, opts['Equivalence'],
                      opts['fraction_in_oxidizer_or_fuel'])
    elif Mixture_type == 'phi_fuel':
        mix_params = (Mixture_type, opts['Equivalence'],
                      opts['fuel_fraction'])
    elif Mixture_type == 'oxi_fuel':
        mix_params = (Mixture_type, opts['oxidizer_fraction'],
                      opts['fuel_fraction'])
    else:
        raise ValueError('Mixture_type = {} is not supported'.format(Mixture_type))

    return mix_params


if __name__ == "__main__":

    # Open and parse the input file
    inputs = list(yaml.safe_load_all(open('input.yaml', 'r')))
    sim_inputs = inputs[0]
    plot_inputs = inputs[1]

    mix_options = sim_inputs['Mixture_options']
    mix_params = package(mix_options)

    # Unpackage other items
    Pressure = mix_options['Pressure']
    Temperature = mix_options['Temperature']
    Fuel = mix_options['Fuel']
    Oxidizer = mix_options['Oxidizer']
    Diluent = mix_options['Diluent']
    Array_type = mix_options['Array_type']
    Mixture_type = mix_options['Mixture_type']

    par = sim_inputs['Parallel']
    Mechanism = sim_inputs['Mechanism']
    Save_files = sim_inputs['Save_files']
    problem_type = sim_inputs['Simulation_Type']

    print(cf.parameters_string(problem_type, Pressure, Temperature, mix_params,
                               Mechanism, Fuel, Oxidizer, Diluent))
    if problem_type == '0D':
        zeroD.run_0D_simulation(Mechanism, Array_type, Pressure, Temperature,
                                Fuel, Oxidizer, Diluent, Mixture_type,
                                mix_params, Save_files, par,
                                **sim_inputs['ZeroD_options'])
    elif problem_type == '1D':
        oneD.run_flame_simulation(Mechanism, Array_type, Pressure, Temperature,
                                  Fuel, Oxidizer, Diluent, mix_params,
                                  Save_files, par, **sim_inputs['Flame_options'])
    else:
        print('Error! Simulation_Type String Does Not Match!' +
              '\nMake sure string matches one of two options!')
        sys.exit()

    if sim_inputs['Plot_on_completion']:
        plot_inputs['Folder_name'] = 2
        if problem_type == '0D':
            zeroD_plot.main(plot_inputs)
        elif problem_type == '1D':
            flame_plot.main(**plot_inputs)
