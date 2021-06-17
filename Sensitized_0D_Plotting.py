# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:46:05 2020

@author: Kodo Bear
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import yaml
import matplotlib.patches as mpl_patches
dirname = os.path.normpath(os.path.dirname(__file__))
plt.style.use(os.path.join(dirname, 'CSULA_Combustion.mplstyle'))


def rxn_plots(max_rxns, species, Array_Type, Plot_path, varnames):
    """
    Plots and saves figures of the reactions that were most sensitive per
    simulation case per independant variable.

    Parameters
    ----------
    max_rxns : list
        Simulation information appended as well as the reaction with the
        maximum sensitivity
    species : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated
    varnames : list
        A list of the four independent variables to plot.

    Returns
    -------
    None.

    """
    for specie in species[0]:
        Temperature  = []
        Pressure     = []
        Equivalence  = []
        Fuel         = []
        Rxns         = []
        Oxygen = []
        num_ticks = 10
        for ents in max_rxns:
            if ents[1] is None:
                continue
            else:
                Rxns.append(ents[1][0])
            Temperature.append(ents[0]['temperature'])
            Pressure.append(ents[0]['pressure'])
            Equivalence.append(ents[0]['phi'])
            Fuel.append(ents[0]['fuel'])
            Oxygen.append(ents[0]['oxidizer'])

        # info is a dictionary with an entry for each figure. The entries are
        # [x label, x points, name]
        info = {'P': ['Pressure [kPa]', Pressure, 'Pressure'],
                'Phi': [r'Equivalence Ratio [$\phi$]', Equivalence,
                        'Equivalence'],
                'F': ['Fuel Mole fraction', Fuel, 'Fuel'],
                'Oxi': ['Oxidizer Mole Fraction', Oxygen, 'Oxygen'],
                'T': ['Temperature [K]', Temperature, 'Temperature']}

        for var in varnames:
            if var not in info.keys():
                fmt = ('WARNING: {} is an invalid variable name for rxn_plots.'
                       ' It must be in {}. Using [P, Phi, F, T]')
                print(fmt.format(var, list(info.keys())))
                varnames = ['P', 'Phi', 'F', 'T']

        fig, axes = plt.subplots(nrows=2, ncols=2)
        # Four subplots of rxn versus independent variables
        for var, ax, string in zip(varnames, axes.flatten(),
                                   ['(a)', '(b)', '(c)', '(d)']):
            if Array_Type == 'log':
                ax.set_xscale('log')
            ax.plot(info[var][1], Rxns)
            ax.annotate(xy=(0.90, 0.90), text=string, xycoords='axes fraction')
            ax.set_xlabel(info[var][0])
            ax.set_ylabel('Reaction #')
            ax.set_ylim([0, int(max(Rxns))+1])
            if max(Rxns) <= num_ticks:
                ax.set_yticks(np.arange(0, max(Rxns)+1, step=2))
            else:
                ax.set_yticks(np.arange(0, max(Rxns),
                                        step=int(max(Rxns)/num_ticks)))
            ax.grid(True)

        fig.suptitle('Time Instance: ' + format(species[3][0]) +
                     ' s    Species: ' + specie)
        fig.savefig(Plot_path + '\\[' + format(specie) +
                    '] Most Sensitive Reactions.png')
        plt.close(fig)

        # Individual plots of Rxn versus independent variable
        for var in varnames:
            plt.figure("Sensitivity of "+specie+" to "+info[var][2])
            if Array_Type == 'log':
                plt.xscale('log')
            plt.plot(info[var][1], Rxns)
            plt.xlabel(info[var][0])
            plt.ylabel('Most Sensitive Reaction Number for '+specie)
            plt.ylim([0, int(max(Rxns)) + 1])
            if max(Rxns) <= 25:
                plt.yticks(np.arange(0, max(Rxns) + 1, step=1))
            else:
                plt.yticks(np.arange(0, max(Rxns), step=int(max(Rxns)/20)))
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white",
                                             ec="white", lw=0, alpha=0)] * 2
            labels = ['Time: ' + format(species[3][0]) + ' s','Species: ' + specie]
            plt.legend(handles, labels, loc='best', fontsize='small',
                       fancybox=True, framealpha=0.7,
                       handlelength=0, handletextpad=0)
            plt.grid(True)
            plt.savefig(Plot_path + '\\[' + format(specie) +
                        '] ' + info[var][2] + ' Most Sensitive Reaction.png')
            plt.close()


def integrated_strength_plots(int_strength, species, Plot_path, rxn_num,
                              Array_Type, varnames):
    """
    Plots and saves figures of the integrated senstivities of reactions of
    interests against each indpendant varaible.

    Parameters
    ----------
    int_strength : list
        This is a list of lists. Each list is the result for one species.
        Within each species is a list of normalized integrated strength for each
        reaction in order. In the future, it may be beneficial to also return
        IS, the non-normalized integrated strength.
    species : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated
    Plot_path : str
        A string of the save path to the plots folder.
    rxn_num : list
        Reaction of interest to plot the senstivities of.
    varnames : list
        A list of the four independent variables to plot.

    Returns
    -------
    None.

    """
    # info is a dictionary with an entry for each figure. The entries are
    # [x label, key for int_strength]
    info = {'P': ['Pressure [kPa]', 'pressure'],
            'Phi': ['Equivalence Ratio', 'phi'],
            'F': ['Fuel Mole Fraction', 'fuel'],
            'T': ['Temperature [K]', 'temperature'],
            'Oxi': ['Oxidizer Mole Fraction', 'oxidizer']}
    file_format = '[{:}] Normalized Integrated Strength.png'.format
    for i in range(len(species[0])):
        spec = species[0][i]
        fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True,
                                 subplot_kw={'ylabel': r'$\hat{IS}_{' + spec + ', j}$'})
        for var, ax, string in zip(varnames, axes.flatten(),
                                   ['(a)', '(b)', '(c)', '(d)']):
            if var not in info.keys():
                ax.set_xlabel(var)
                key = None
            else:
                ax.set_xlabel(info[var][0])
                key = info[var][1]
            if Array_Type == 'log':
                ax.set_xscale('log')
            ax.grid(True)
            x = []
            y = []
            for condition in int_strength:
                if condition[1][0] is None:
                    continue
                else:
                    if var not in info.keys():
                        # Assume var is an input species
                        x.append(condition[0]['mixture'][var])
                    else:
                        x.append(condition[0][key])
                    norm_IS = condition[1][i]
                    y.append(norm_IS[rxn_num])
            ax.plot(x, y)
            ax.axhline(0, c='gray', ls='--', marker="None")
            ax.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        fig.suptitle(species[2][0][rxn_num])
        fname = file_format(spec)
        fig.savefig(os.path.join(Plot_path, fname))
        plt.close(fig)


def integrated_strength_rxn_plots(int_strength, species, Plot_path,
                                  rxn_num, threshold, Array_Type):
    """
    Plots and saves figures of the integrated senstivities of reactions of
    interests that are above the threshold against each indpendant varaible.

    Parameters
    ----------
    int_strength : list
        This is a list of lists. Each list is the result for one species.
        Within each species is a list of normalized integrated strength for each
        reaction in order. In the future, it may be beneficial to also return
        IS, the non-normalized integrated strength.
    species : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated
    Plot_path : str
        A string of the save path to the plots folder.
    rxn_num : list
        Reaction of interest to plot the senstivities of.
    threshold : float
        Minimum sensitivity strength defined by user

    Returns
    -------
    None.

    """
    # info is a dictionary with an entry for each figure. The entries are
    # [x label, key for int_strength]
    info = {'P': ['Pressure [kPa]', 'pressure'],
            'phi': ['Equivalence Ratio', 'equivalence'],
            'X': ['Fuel Mole Fraction', 'fuel'],
            'T': ['Temperature [K]', 'temperature']}
    file_format = '[{:}] Normalized Integrated Strength Plots.png'.format
    conditions = [ ('T', 'phi'), ('X', 'phi'), ('P', 'phi'),
                  ('T', 'P'), ('T', 'X'), ('P', 'X')]
    anlabel = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)']

    for i in range(len(species[0])):
        spec = species[0][i]
        fig, axes = plt.subplots(nrows=2, ncols=3)
        for condition, ax, string in zip(conditions, axes.flatten(), anlabel):
            ax.set_xlabel(info[condition[0]][0])
            ax.set_ylabel(info[condition[1]][0])
            if Array_Type == 'log':
                ax.set_xscale('log')
                ax.set_yscale('Log')
            ax.grid(True)
            x = []
            y = []
            x_key = info[condition[0]][1]
            y_key = info[condition[1]][1]
            if y_key == 'equivalence':
                y_key = 'phi'
            for cond in int_strength:
                if cond[1][0] is None:
                    continue
                elif abs(cond[1][i][rxn_num]) >= threshold:
                    x.append(cond[0][x_key])
                    y.append(cond[0][y_key])
            ax.plot(x, y)
            ax.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        fig.suptitle(species[2][0][rxn_num]+
                     '    Threshold >='+format(threshold))
        fname = file_format(spec)
        fig.savefig(os.path.join(Plot_path, fname))
        plt.close(fig)


def reaction_of_interest_plot(max_rxns, species, rxn_number, Array_Type, Plot_path):
    """
    Plots and saves figures comparing independant variables that were most
    senstivity.

    Parameters
    ----------
    max_rxns : list
        Simulation information appended as well as the reaction with the
        maximum sensitivity
    species : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated
    rxn_num : list
        Reaction of interest to plot the senstivities of.
    Plot_path : str
        A string of the save path to the plots folder.

    Returns
    -------
    None.

    """
    for specie in species[0]:
        Temperature  = []
        Pressure     = []
        Equivalence  = []
        Fuel         = []
        for ents in max_rxns:
            if ents[1] is None:
                continue
            elif ents[1][0] == rxn_number:
                Temperature.append(ents[0]['temperature'])
                Pressure.append(ents[0]['pressure'])
                Equivalence.append(ents[0]['phi'])
                Fuel.append(ents[0]['fuel'])

        fig, axes = plt.subplots(nrows=2, ncols=3)
        ax = axes.flatten()

        cond_dict = {'T': [Temperature, 'Temperature [K]', 'Temperature'],
                     'P': [Pressure, 'Pressure [kPa]', 'Pressure'],
                     'X': [Fuel, 'Fuel [mole fraction]', 'Fuel'],
                     'phi': [Equivalence, 'Equivalence Ratio [$\phi$]',
                             'Equivalence Ratio']}
        conditions = [('T', 'P'), ('T', 'X'), ('T', 'phi'),
                      ('P', 'X'), ('P', 'phi'), ('X', 'phi')]

        #Six subplots of unique paired independent variables
        for a, condition, string in zip(ax, conditions,
                                        ['(a)', '(b)', '(c)',
                                         '(d)', '(e)', '(f)']):
            x_key = condition[0]
            y_key = condition[1]
            a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
            a.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
            a.set_xlabel(cond_dict[x_key][1])
            a.set_ylabel(cond_dict[y_key][1])
            if Array_Type == 'log':
                a.set_xscale('log')
                a.set_yscale('log')
            a.grid(True)

        fig.suptitle('Reaction Number '+format(rxn_number)+
                      ' '+species[2][0][rxn_number]+
                      '\nTime Instance: '+format(species[3][0])+' s')
        plt.savefig(Plot_path+'\\['+specie+'] Reaction Number '
                    +format(rxn_number)+'.png')
        plt.close(fig)

        #Individual plots for unique paired independant variables
        for cond in conditions:
            x_key = cond[0]
            y_key = cond[1]
            plt.figure()
            plt.plot(cond_dict[x_key][0], cond_dict[y_key][0])
            plt.xlabel(cond_dict[x_key][1])
            plt.ylabel(cond_dict[y_key][1])
            if Array_Type == 'log':
                plt.xscale('log')
                plt.yscale('log')
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white",
                                             ec="white", lw=0, alpha=0)] * 2
            labels = ['Time: '+format(species[3][0]) + ' s','Species: '+specie]
            plt.legend(handles, labels, loc='best', fontsize='small',
                       fancybox=True, framealpha=0.7,
                       handlelength=0, handletextpad=0)
            plt.grid(True)
            plt.savefig(Plot_path+'\\['+specie+'] Reaction Number '
                        +format(rxn_number)+' '+cond_dict[y_key][2]+
                        ' vs '+cond_dict[x_key][2]+'.png')
            plt.close()


def create_csv(int_strength, out_path, rxn_num, spec):
    """
    Creates a csv of integrated sensitivity for requested species with
    thermodynamic properties and mixture information.

    Parameters
    ----------
    int_strength : list
        This is a list of lists. Each list is the result for one species.
        Within each species is a list of normalized integrated strength for each
        reaction in order. In the future, it may be beneficial to also return
        IS, the non-normalized integrated strength.
    out_path : str
        A string of the save path to the simulation folder.
    rxn_num : list
        Reaction of interest to plot the senstivities of.
    spec : list
        A list of specific species of interest in which the sensitivity to
        each reaction per time step is calculated

    Returns
    -------
    None.

    """
    output = []
    for item in int_strength:
        mix = item[0]
        IS = item[1][spec][rxn_num]
        if not np.isnan(IS):
            output.append([mix['temperature'], mix['pressure'], mix['equivalence'],
                           mix['fuel'], IS])
    output = np.array(output)
    out_sorted = output[output[:, 4].argsort()]  # Sort by IS in decreasing order
    header = 'Temperature [K], Pressure [kPa], Phi, Fuel mole fraction, IS_norm'
    np.savetxt(os.path.join(out_path, 'Integrated Strength.csv'),
               out_sorted, delimiter=', ', header=header)


def main(plot_inputs):
    """ Main plotting function. """
    Rxn_interest = plot_inputs['Rxn_interest']
    Four_Plot = plot_inputs['Four_Plot']
    Threshold = plot_inputs['Threshold']
    Folder_name = plot_inputs['Folder_name']

    if Folder_name == '' or Folder_name == 1:
        try:
            with open('last run 0d.pkl', 'rb') as f:
                Folder_name = pickle.load(f)
        except FileNotFoundError:
            print('last run 0d.pkl is missing - I do not know which folder was analyzed most recently.')
            Folder_name = 2
    if Folder_name == 2:
        print('Plotting the most recent simulations.')
        folders = os.listdir('0D_Sensitivity_Results')
        # Assuming this is always sorted in ascending order...
        print(folders)
        Folder_name = [x for x in folders if x[:2] == '20'][-1]
    print('Loading ' + Folder_name)

    # Save the loaded folder name
    with open('last run 0d.pkl', 'wb') as f:
        pickle.dump(Folder_name, f)

    Load_path   = '0D_Sensitivity_Results\\'+Folder_name
    Plot_path   = Load_path+'\\ZeroD_Sensitivity_Plots'

    file = open(Load_path+'\\Case Description.txt','r')
    Description = file.read()
    file.close()
    print(Description)

    Species_Rxn  = pickle.load(open(Load_path+'\\Species_Rxn.pkl', 'rb'))
    Max_Sens_Rxn = pickle.load(open(Load_path+'\\Max_Sens_Rxn.pkl', 'rb'))
    Case_Params  = pickle.load(open(Load_path+'\\Case_Parameters.pkl', 'rb'))
    Array_type = Species_Rxn[4][0]

    with open(os.path.join(Load_path, 'Integrated Strength.pkl'), 'rb') as f:
        Int_Strength = pickle.load(f)

    rxn_plots(Max_Sens_Rxn, Species_Rxn, Array_type, Plot_path, Four_Plot)
    for roi in Rxn_interest:
        reaction_of_interest_plot(Max_Sens_Rxn, Species_Rxn, roi, Array_type, Plot_path)
        integrated_strength_plots(Int_Strength, Species_Rxn, Plot_path, roi,
                                  Array_type, Four_Plot)
        integrated_strength_rxn_plots(Int_Strength, Species_Rxn, Plot_path,
                                      roi, Threshold, Array_type)


if __name__ == "__main__":
    # Open and parse the input file
    inputs = list(yaml.safe_load_all(open('input.yaml', 'r')))
    main(inputs[1])
