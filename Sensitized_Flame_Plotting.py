# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:20:01 2020

@author: Kodo Bear
"""

import os
import sys
import csv
import numpy
import pickle
import yaml
import datetime
import shutil
from operator import itemgetter
import matplotlib.pyplot as plt
import Sensitized_Flame_Experiment as experiment
dirname = os.path.normpath(os.path.dirname(__file__))
plt.style.use(os.path.join(dirname, 'CSULA_Combustion.mplstyle'))
#from tkinter import filedialog

def rxn_plots(f_info, save_path, log, conditions=['P', 'F', 'Phi', 'O2']):
    """
    Plots and saves figures of the reactions that were most sensitive per
    simulation case per independant variable

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    save_path : str
        A string of the save path to the plots folder.
    log : boolean
        If true plots will use a logarithmic scale.
        If false plots will use a linear scale
    conditions : list, optional
        A list of variables to plot the sensitivities against
        The default is ['P', 'F', 'Phi', 'O2'].

    Returns
    -------
    None.

    """
    Tint     = f_info[0]['Conditions'][0]
    Fue      = f_info[0]['Conditions'][6]
    Oxi      = f_info[0]['Conditions'][7]
    Pressure = []
    Phi      = []
    Fuel     = []
    Oxidizer = []
    Max_rxn  = []
    for f in f_info:
        max_rxn  = max_sens(f)
        max_num  = max_rxn[0]
        Max_rxn.append(max_num+1)
        Pressure.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Fuel.append(f['Conditions'][9])
        Oxidizer.append(f['Conditions'][10])

    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)
    ax        = axes.flatten()
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel Mole Fraction '+str(Fue)],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'O2': [Oxidizer, 'Oxygen Mole Fraction '+str(Oxi)],
                 'T': [Tint, 'Temperature [K]'],
                 'Rxn': [Max_rxn, 'Reaction #']}

    for var in conditions:
        if var not in cond_dict.keys():
            fmt = ('WARNING: {} is an invalid variable name for rxn_plots.'
                   ' It must be in {}. Using [P, F, Phi, O2]')
            print(fmt.format(var, list(cond_dict.keys())))
            conditions = ['P', 'F', 'Phi', 'O2']
    for a, condition, string in zip(ax, conditions,
                                    ['(a)', '(b)', '(c)', '(d)']):
        x_key = condition
        y_key = 'Rxn'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
        a.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        a.set_xlabel(cond_dict[x_key][1])
        a.set_ylabel(cond_dict[y_key][1])
        if log:
            a.set_xscale('log')
        fig.suptitle('Initial Temperature: {:.0f} K'.format(Tint))
        plt.savefig(save_path +
                    '\\Max Reactions with Initial Temperature {:.0f} K.png'.format(Tint))
    plt.close(fig)


def rxn_strength_plots(f_info, rxn_int, nrxns, threshold, save_path, log,
                       conditions_a=['P', 'F', 'Phi', 'O2']):
    """
    Plots and saves figures of a specific reactions with its normalized
    sensitivity versus four specified conditions.

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    rxn_int : list
        A list of reactions of interest to plot.
    nrxns : int
        An integer that defines the top n-reactions that are used to define the
        normalized sensitivity
    threshold : float
        An integer representing the minimum normalize sensitiity.
    save_path : str
        A string of the save path to the plots folder.
    log : boolean
        If true plots will use a logarithmic scale.
        If false plots will use a linear scale
    conditions_a : list, optional
        A list of variables to plot the sensitivities agains
        The default is ['P', 'F', 'Phi', 'O2'].

    Returns
    -------
    None.

    """
    Rxn_Eq     = format(f_info[0]['Flame'][0][rxn_int][2])
    Rxn_Number = format(f_info[0]['Flame'][0][rxn_int][0]+1)
    #List for strengths
    Pressure      = []
    Phi           = []
    Fuel          = []
    Oxygen        = []
    Su            = []
    Sens_strength = []
    Tint          = []
    #Lists for above threshold
    P_threshold        = []
    Phi_threshold      = []
    F_threshold        = []
    Oxygen_threshold   = []
    Su_threshold       = []
    Sens_str_threshold = []
    T_threshold        = []

    species_of_interest = {}
    for cond in conditions_a:
        if cond not in ['T', 'F', 'Phi', 'O2', 'Su', 'P']:
            species_of_interest[cond] = []

    for f in f_info:
        average_nrxns = rxn_average(f, nrxns)
        strength = f['Flame'][0][rxn_int][1]/average_nrxns
        Sens_strength.append(strength)
        Tint.append(f['Conditions'][0])
        Pressure.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Fuel.append(f['Conditions'][9])
        Oxygen.append(f['Conditions'][10])
        Su.append(f['Flame'][1])
        for cond in conditions_a:
            if cond not in ['T', 'F', 'Phi', 'O2', 'Su', 'P']:
                # This must be a species name. search for that species
                try:
                    species_of_interest[cond].append(f['Conditions'][5][cond])
                except KeyError:
                    species_of_interest[cond].append(0)

    s_of_interest = {k: [v, k + ' Mole Fraction', None] for k, v in
                     species_of_interest.items()}

    for n in range(len(Sens_strength)):
        if abs(Sens_strength[n]) >= threshold:
            P_threshold.append(f_info[n]['Conditions'][1])
            Phi_threshold.append(f_info[n]['Conditions'][2])
            F_threshold.append(f_info[n]['Conditions'][9])
            Oxygen_threshold.append(f_info[n]['Conditions'][10])
            Su_threshold.append(f['Flame'][1])
            Sens_str_threshold.append(Sens_strength[n])
            T_threshold.append(f_info[n]['Conditions'][0])

    #Dictionary organized as follow 'Key': [Data, axis-label, Data_threshold]
    cond_dict = {'P': [Pressure, 'Pressure [atm]', P_threshold],
                 'F': [Fuel, 'Fuel Mole Fraction', F_threshold],
                 'Phi': [Phi, 'Equivalence Ratio [$\phi$]', Phi_threshold],
                 'O2': [Oxygen, 'Oxygen Mole Fraction', Oxygen_threshold],
                 'T': [Tint, 'Temperature [K]', T_threshold],
                 'Su': [Su, 'Flame Speed [m/s]', Su_threshold],
                 'Strength': [Sens_strength, r'$\hat S_{'+Rxn_Eq+'}$',
                              Sens_str_threshold]}
    cond_dict = {**cond_dict, **s_of_interest}  # Add species mole fractions

    figa, axesa = plt.subplots(nrows=2, ncols=2, sharey=True)
    axa         = axesa.flatten()

    T_initials = set(Tint)
    T_fmt = 'Initial Temperature = ' + ', '.join(['{:.0f} K' for x in T_initials])
    for a, condition, string in zip(axa, conditions_a,
                                    ['(a)', '(b)', '(c)', '(d)']):
        x_key = condition
        y_key = 'Strength'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
        a.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        a.set_xlabel(cond_dict[x_key][1])
        a.set_ylabel(cond_dict[y_key][1])
        if log:
            a.set_xscale('log')

        figa.suptitle(('Reaction {}\n' + T_fmt).format(Rxn_Eq, *T_initials))
        plt.savefig(save_path+'\\Reaction '+Rxn_Number+
                    ' Average Strength against top '+format(nrxns)+' rxns.png')
    plt.close(figa)

    if not len(Sens_str_threshold) == 0:
        figb, axesb  = plt.subplots(nrows=2, ncols=3)
        axb          = axesb.flatten()
        conditions_b = [('P', 'F'), ('P', 'Phi'), ('P','O2'),
                        ('F', 'Phi'), ('F', 'O2'), ('O2', 'Phi')]
        for a, condition in zip(axb, conditions_b):
            x_key = condition[0]
            y_key = condition[1]
            a.plot(cond_dict[x_key][2], cond_dict[y_key][2])
            a.set_xlabel(cond_dict[x_key][1])
            a.set_ylabel(cond_dict[y_key][1])
            if log:
                a.set_xscale('log')
            if log:
                a.set_yscale('log')

            figb.suptitle(('Reaction {}\n Threshold >= {:.1f}, ' + T_fmt).format(Rxn_Eq, threshold, *T_initials))
            fname = ('Reaction {} Average Strength with ' + T_fmt + '.png').format(Rxn_Number, *T_initials)
            plt.savefig(os.path.join(save_path, fname))
        plt.close(figb)
    else:
        print('Reaction: '+str(f_info[0]['Flame'][0][rxn_int][2])+
              ' shows no cases where the sensitiviy is above threshold '+
              str(threshold))

    figc, axc = plt.subplots(nrows=1, ncols=1)
    x_key = 'Su'
    y_key = 'Strength'
    axc.plot(cond_dict[x_key][0], cond_dict[y_key][0])
    axc.set(xlabel=cond_dict[x_key][1], ylabel=cond_dict[y_key][1],
            title='Normalized Sensitivity of Reaction '+Rxn_Eq+
            ' vs. Flame Speed')
    fname = ('Reaction {} Normalized Sensitivity vs Flame speed with ' + T_fmt + '.png').format(Rxn_Number, *T_initials)
    figc.savefig(os.path.join(save_path, fname))
    plt.close(figc)


def rxn_interest_plots(f_info, rxn_int, save_path, log):
    """
    Plots and save figures where reactions of interest where most senstivity
    per independant varaible

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    rxn_int : list
        A list of reactions of interest to plot.
    save_path : str
        A string of the save path to the plots folder.
    log : boolean
        If true plots will use a logarithmic scale.
        If false plots will use a linear scale

    Returns
    -------
    None.

    """
    Pressure = []
    Phi      = []
    Fuel     = []
    Oxygen   = []
    Tint     = f_info[0]['Conditions'][0]
    Rxn_name = f_info[0]['Flame'][0][rxn_int][2]
    Rxn_num  = f_info[0]['Flame'][0][rxn_int][0]+1
    for f in f_info:
        Rxn_sens = f['Flame'][0][rxn_int][1]
        Max_rxn_sens = max_sens(f)
        if Rxn_sens == Max_rxn_sens[1]:
            Pressure.append(f['Conditions'][1])
            Phi.append(f['Conditions'][2])
            Fuel.append(f['Conditions'][9])
            Oxygen.append(f['Conditions'][10])
    if not len(Phi) == 0:
        fig, axes = plt.subplots(nrows=2, ncols=3)
        ax        = axes.flatten()
        #Dictionary organized as follow 'Key': [Data, axis-label]
        cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                     'F': [Fuel, 'Fuel Mole Fraction'],
                     'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                     'O2':[Oxygen, 'Oxygen Fuel Mole Fraction'],
                     'T': [Tint, 'Temperature [K]']}
        conditions = [('P', 'F'), ('P', 'Phi'), ('P','O2'),
                      ('F', 'Phi'), ('F', 'O2'), ('O2', 'Phi')]

        #Three subplots of unique paired independent variables
        for a, condition in zip(ax, conditions):
            x_key = condition[0]
            y_key = condition[1]
            a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
            a.set_xlabel(cond_dict[x_key][1])
            a.set_ylabel(cond_dict[y_key][1])
            if log:
                a.set_xscale('log')
            if log:
                a.set_yscale('log')
            fig.suptitle('Reaction Number: '+format(Rxn_num)+
                          ' Reaction Name: '+Rxn_name+
                          '\nInitial Temperature: '+format(Tint, '.0f')+' K')
            plt.savefig(save_path+'\\Reaction Number '+format(Rxn_num)+
                        ' Max Sensitivity for Parameters with'
                        ' Initial Temperature '+format(Tint, '.0f')+'K.png')
        plt.close(fig)
    else:
        print('Reaction Number '+str(Rxn_num)+', Reaction: '+str(Rxn_name)+
              ' shows no cases where it is most sensitive reaction')


def flame_speed_plots(f_info, save_path, log):
    """
    Plots and saves figures of simulations where the flame speed is above
    the minumum flame speed defined by the user. The plots show two independant
    variables where the conditions have been met.

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    save_path : str
        A string of the save path to the plots folder.
    log : boolean
        If true plots will use a logarithmic scale.
        If false plots will use a linear scale

    Returns
    -------
    None.

    """
    Tint     = f_info[0]['Conditions'][0]
    Pressure = []
    Phi      = []
    Fuel     = []
    Oxidizer = []
    Su       = []
    Flame_T  = []
    for s in f_info:
        Pressure.append(s['Conditions'][1])
        Phi.append(s['Conditions'][2])
        Fuel.append(s['Conditions'][9])
        Oxidizer.append(s['Conditions'][10])
        Su.append(s['Flame'][1])
        Flame_T.append(s['Flame'][3])

    fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True)
    ax        = axes.flatten()
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel Mole Fraction'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'Oxi':[Oxidizer, 'Oxygen Mole Fraction'],
                 'T': [Tint, 'Temperature [K]'],
                 'Su': [Su, 'Su [m/s]'],
                 'Flame_T': [Flame_T, 'Flame Temperature [K]']}
    conditions = ['P', 'F', 'Phi', 'Oxi']
    for a, condition in zip(ax, conditions):
        x_key = condition
        y_key = 'Su'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
                marker='o', mfc='none', mec='k')
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
        a.set_xlabel(cond_dict[x_key][1])
        a.set_ylabel(cond_dict[y_key][1])
        if log:
            a.set_xscale('log')
        fig.suptitle('Initial Temperature: '+format(Tint, '.0f')+' K')
        plt.savefig(save_path+'\\Flame Speed vs. Independant Variables with'
                    ' Initial Temperature '+format(Tint, '.0f')+'K.png')
    plt.close(fig)


def max_rxn_csv(f_info, save_path):
    """
    Creates a two csv files. The Maximum Rxns.csv stores information about
    reactions that were most sensitivy and the number of cases this was true.
    The Maximum Rxns Parameters.csv stores information of the reactions that
    was most sensitivy per case with case information.

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    save_path : str
        A string of the save path to the plots folder.

    Returns
    -------
    Max_rxn_params : dict
        Information of the reaction that was most sensitivy and the parameters
        are stored per case in a list.
    Max_rxns : list
        A list of the most sensitive reactions numbers and the number of cases
        where they were the most sensitive.
    Max_rxn : list
        A list of the reaction that was most sensitive providing its reaction
        number, reaction equation, and sensitivity per case.

    """
    Tint     = f_info[0]['Conditions'][0]
    Pressure = []
    Phi      = []
    Fuel     = []
    Oxygen   = []
    Max_rxn  = []
    for f in f_info:
        max_rxn  = max_sens(f)
        Max_rxn.append([max_rxn[0], max_rxn[2], max_rxn[1]])
        Pressure.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Fuel.append(f['Conditions'][9])
        Oxygen.append(f['Conditions'][10])
    Max_rxn_params = {}
    Max_rxns   = {}
    for n in range(len(Max_rxn)):
        rxn_num = str(Max_rxn[n][0]+1)
        if rxn_num not in Max_rxn_params.keys():
            Max_rxn_params[rxn_num] = [[Max_rxn[n][1], Tint],
                                       [Pressure[n], Phi[n],
                                        Fuel[n], Oxygen[n]]]
            Max_rxns[rxn_num]       = [Max_rxn[n][1], 1, Tint]
        elif rxn_num in Max_rxn_params.keys():
            Max_rxn_params[rxn_num].append([Pressure[n], Phi[n],
                                            Fuel[n], Oxygen[n]])
            Max_rxns[rxn_num][1] += 1
        else:
            print('Something went wrong!')
    csv_file_max = os.path.join(save_path, 'Maximum Rxns.csv')
    csv_file_params = os.path.join(save_path, 'Maximum Rxns Parameters.csv')
    with open(csv_file_max, mode='w') as max_rxns:
        mr_writer = csv.writer(max_rxns, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        mr_writer.writerow(['Rxn Number', 'Rxn Equation',
                            'Max Count', 'Initial Temperature [K]'])
        for keys in Max_rxns:
            mr_writer.writerow([keys, Max_rxns[keys][0],
                                Max_rxns[keys][1], Max_rxns[keys][2]])
    with open(csv_file_params, mode='w') as rxn_params:
        rp_writer = csv.writer(rxn_params, delimiter=',',
                               quotechar='"', quoting=csv.QUOTE_MINIMAL)
        rp_writer.writerow(['Rxn Number', 'Rxn Equation',
                            'Initial Temperature [K]'])
        for keys in Max_rxn_params:
            rp_writer.writerow([keys, Max_rxn_params[keys][0][0],
                                Max_rxn_params[keys][0][1]])
            rp_writer.writerow(['', 'Parameters:','Pressure [atm]',
                                'Equivalence Ratio','Fuel Mole Fraction',
                                'Oxygen Mole Fraction'])
            for p in range(1,len(Max_rxn_params[keys])):
                rp_writer.writerow(['', '', Max_rxn_params[keys][p][0],
                                    Max_rxn_params[keys][p][1],
                                    Max_rxn_params[keys][p][2],
                                    Max_rxn_params[keys][p][3]])

    return Max_rxn_params, Max_rxns, Max_rxn


def top_nrxns_csv(f_info, nrxns, save_path):
    """
    Creates a csv file the provides information about the most sensitivy
    reaction and the parameters per simulation case.

    Parameters
    ----------
    f_info : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    nrxns : int
        An integer that defines the top n-reactions that are used to define the
        normalized sensitivity
    save_path : str
        A string of the save path to the plots folder.

    Returns
    -------
    Top_Rxns_List : list
        An appened list of per simulation case information with the structure:
            [[f['Conditions'][1], f['Conditions'][2],
             f['Conditions'][9], f['Conditions'][10],
             top_rxns(f, nrxns)], [...], ...]

    """
    Tint          = f_info[0]['Conditions'][0]
    Top_Rxns_List = []
    for f in f_info:
        Top_Rxns_List.append([f['Conditions'][1], f['Conditions'][2],
                              f['Conditions'][9], f['Conditions'][10],
                              top_rxns(f, nrxns)])
    csv_file_top_nrxns = os.path.join(save_path, 'Top Rxns per Case.csv')
    with open(csv_file_top_nrxns, mode='w') as top_nrxns:
        tr_writer = csv.writer(top_nrxns, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for i in Top_Rxns_List:
            tr_writer.writerow(['Initial Temperature [K]','Pressure [atm]',
                                'Equivalence Ratio','Fuel Mole Fraction',
                                'Oxygen Mole Fraction'])
            tr_writer.writerow([Tint, i[0], i[1], i[2], i[3]])
            tr_writer.writerow(['','Rxns:','Rxn #',
                                'Sensitivity Strength','Rxn Equation'])
            for j in i[4]:
                tr_writer.writerow(['', '', j[0]+1, j[1], j[2]])
    return Top_Rxns_List


def average_sens_csv(flame, nrxns, save_path):
    """
    Creates a csv providing the average sensitivity of reactions across the
    entire simulation from the first to last case. Two csvs are created one
    of just the top n-reactions and the other of all reactions in the flame
    information sorted in a descending order.

    Parameters
    ----------
    flame : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    nrxns : int
        An integer that defines the top n-reactions that are used to define the
        normalized sensitivity
    save_path : str
        A string of the save path to the plots folder.

    Returns
    -------
    sens_average : list
        A list of reaction information in a descending order of average
        sensitivity.
    topnrxnnumbers : list
        The same list as sens_average but on the top n-reactions are listed.

    """
    values = []
    topnrxnnumbers = []
    for f in flame[0]['Flame'][0]:
        values.append((f[0], 0.0, f[2]))
    dtype = [('Rxn Number', int), ('ANSens', float), ('Rxn Equation', object)]
    sens_average = numpy.array(values, dtype=dtype)
    for i in flame:
        for j in range(len(sens_average)):
            sens_average[j][1] += abs(i['Flame'][0][j][1])/rxn_average(i, nrxns)
    for m in sens_average:
        m[1] = m[1]/len(flame)
    sens_average = numpy.sort(sens_average, order='ANSens')
    sens_average = sens_average[::-1]
    for n in range(nrxns):
        topnrxnnumbers.append(sens_average[n][0])
    for k in range(len(sens_average)):
        sens_average[k][0] += 1
    topnrxnssens = sens_average[0:nrxns]
    all_sp = os.path.join(save_path,
                                 "All Rxns Total Average Sensitivities.csv")
    topn_sp = os.path.join(save_path,
                           "Top "+str(nrxns)+" Rxn Total Average Sensitivities.csv")
    csvformat = ['%d', '%f', '%s']
    csvheader = 'Rxn Number, Average Normalized Sensitivity, Rxn Equation'
    numpy.savetxt(all_sp, sens_average, delimiter=',', fmt=csvformat,
                  header=csvheader, comments='')
    numpy.savetxt(topn_sp, topnrxnssens, delimiter=',', fmt=csvformat,
                  header=csvheader, comments='')
    return sens_average, topnrxnnumbers


def rxn_average(flame, nrxns):
    """
    Calculate the average sensitivity of the top n-reactions for the given
    simulation case.

    Parameters
    ----------
    flame : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    nrxns : int
        An integer that defines the top n-reactions that are used to define the
        normalized sensitivity

    Returns
    -------
    average_sens : float
        Calculate average senstivity of the top n-reactions.

    """
    rxn_sens = []
    for s in flame['Flame'][0]:
        rxn_sens.append(abs(s[1]))
    rxn_sens.sort()
    top_nrxns    = rxn_sens[-nrxns:]
    average_sens = (sum(top_nrxns)/len(top_nrxns))
    return average_sens


def max_sens(flame):
    """
    Determines the most senstivy reaction in a simulation case and provides
    the reaction number, reaction sensitivity, and the reaction equation.

    Parameters
    ----------
    flame : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}

    Returns
    -------
    max_rxn : list
        A list of information of the reaction that was most senstive given in
        the structure:
            [Reaction Number, Reaction Sensitivity, Reaction Equation]

    """
    max_rxn_num  = 0
    max_rxn_sens = 0
    for m in flame['Flame'][0]:
        if abs(m[1]) > abs(max_rxn_sens):
            max_rxn_num  = m[0]
            max_rxn_sens = m[1]
            max_rxn_eq   = m[2]
    max_rxn = [max_rxn_num, max_rxn_sens, max_rxn_eq]
    return max_rxn


def top_rxns(flame, nrxns):
    """
    Creates a list of the most sensitivy n-reactions in descending order.

    Parameters
    ----------
    flame : dict
        Information of the results of the simulation as well as simulation
        conditions with the following structure:
         {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
          'Conditions': [T, p, phi, Fuel, Oxidizer, mix,
                         Fuel_name, Oxidizer_name, Diluent_name,
                         Fue_Percent, Oxi_Percent, Dil_Percent,
                         at]}
    nrxns : int
        An integer that defines the top n-reactions that are used to define the
        normalized sensitivity

    Returns
    -------
    top_rxns_list : list
        A sorted list of reactions and reaction information in order of
        sensitivity in a descending order.

    """
    f = flame['Flame'][0]
    f_sort = []
    for i in f:
        f_sort.append(i)
    for j in range(len(f_sort)):
        f_sort[j][1] = abs(f_sort[j][1])
    f_sort = sorted(f_sort, key=itemgetter(1))
    top_rxns_list = f_sort[-nrxns:]
    sens_sum = 0
    for k in top_rxns_list:
        sens_sum += k[1]
    rxn_str = sens_sum/nrxns
    for n in top_rxns_list:
        for m in f:
            if m[0] == n[0]:
                n[1] = m[1]/rxn_str
    return top_rxns_list


def check_compatible(folder, file1, file2):
    f1 = os.path.join(folder, file1)
    f2 = os.path.join(folder, file2)

    input1 = list(yaml.safe_load_all(open(f1, 'r')))[0]
    input2 = list(yaml.safe_load_all(open(f2, 'r')))[0]

    assert input1['Simulation_Type'] == input2['Simulation_Type']
    assert input1['Mechanism'] == input2['Mechanism']
    assert input1['Mixture_options']['Fuel'] == input2['Mixture_options']['Fuel']
    assert input1['Mixture_options']['Oxidizer'] == input2['Mixture_options']['Oxidizer']
    assert input1['Mixture_options']['Diluent'] == input2['Mixture_options']['Diluent']


def main(Folder_name, Rxn_interest, Four_Plot, Min_speed, Nrxns, Threshold):
    """Main plotting script."""
    multiple = False
    if Folder_name == '' or Folder_name == 1:
        try:
            with open('last run 1d.pkl', 'rb') as f:
                Folder_name = pickle.load(f)
        except FileNotFoundError:
            print('last run 1d.pkl is missing - I do not know which folder was analyzed most recently.')
            Folder_name = 2
    elif Folder_name == 2:
        print('Plotting the most recent simulations.')
        folders = os.listdir('Flame_Sensitivity_Results')
        # Assuming this is always sorted in ascending order...
        Folder_name = [x for x in folders if x[:2] == '20'][-1]
    elif type(Folder_name) is list:
        # Combine and load multiple simulations.
        multiple = True
        # Working directory
        now = datetime.datetime.now()
        directory = now.strftime("%Y_%m_%d %H.%M.%S Combined_sens_plots")
        Load_path = os.path.join('Flame_Sensitivity_Results', directory)
        figure_path = os.path.join(Load_path, 'Flame_Sensitivity_Plots')
        os.makedirs(figure_path)

        shutil.copyfile('input.yaml', os.path.join(Load_path, 'plot_input.yaml'))

        count = 0
        Flame_info = []
        for load_dir in Folder_name:
            full_dir = os.path.join('Flame_Sensitivity_Results', load_dir)
            count += 1
            shutil.copyfile(os.path.join(full_dir, 'input.yaml'),
                            os.path.join(Load_path, 'input ' + str(count) + '.yaml'))

            if count > 1:
                # Check that most recent folder is compatible
                check_compatible(Load_path, 'input 1.yaml', 'input ' + str(count) + '.yaml')

            one_set = experiment.collect_flame_info(full_dir)
            Flame_info.extend(one_set)

    print('Loading ' + str(Folder_name))

    # Save the loaded folder name
    with open('last run 1d.pkl', 'wb') as f:
        pickle.dump(Folder_name, f)

    if not multiple:
        Load_path = os.path.join('Flame_Sensitivity_Results', Folder_name)
        # Open up text file with description of simulation
        with open(os.path.join(Load_path, 'Case Description.txt'), 'r') as f:
            print(f.read())

        # Import flame file found in corresponding folder
        with open(os.path.join(Load_path, 'Flame Information.pkl'), 'rb') as f:
            Flame_info = pickle.load(f)

    Plot_path = os.path.join(Load_path, 'Flame_Sensitivity_Plots')

    # Create two lists of flame and no_flame created from flame_info
    Flame = []
    No_flame = []
    for x in Flame_info:
        if x['Flame'][0] is None:
            No_flame.append(x)
        else:
            Flame.append(x)

    # If Flame list is empty plotting script will not occur
    if len(Flame) == 0:
        print('\nNo Conditions Produced a Flame!')
        sys.exit()

    # Creates list that only appends information where flame speed is above min
    Flame_speed_filter = []
    for x in Flame:
        if x['Flame'][1] >= Min_speed:
            Flame_speed_filter.append(x)

    # Note Flame and No_flame are dictionaries
    # {'Flame': [flame_sens, Su, flame_rho, flame_T, mg, ms],
    #  'Conditions': [Tin, P, Phi, Fuel, Oxidizer, Mix,
    #                 Fuel_name, Oxidizer_name, Diluent_name,
    #                 Fue_Percent, Oxi_Percent, Dil_Percent, at]}

    # Plot Functions
    array_type   = Flame[0]['Conditions'][12]
    if array_type == 'log':
        Logspace = True # If True all plots will use logspace
    elif array_type == 'lin':
        Logspace = False # If False all plots will use linspace
    else:
        print('Error! invalid string for array_type.')

    Max_rxn_cond, Max_rxns_dict, M = max_rxn_csv(Flame_speed_filter, Load_path)
    T_Rxn_List = top_nrxns_csv(Flame_speed_filter, Nrxns, Load_path)
    rxn_plots(Flame_speed_filter, Plot_path, Logspace)
    flame_speed_plots(Flame_speed_filter, Plot_path, Logspace)
    Average_Sensitivities, Topnrxns = average_sens_csv(Flame_speed_filter,
                                                       Nrxns, Load_path)
    if not len(Rxn_interest) == 0:
        for Rxns in Rxn_interest:
            rxn_strength_plots(Flame_speed_filter, Rxns, Nrxns,
                               Threshold, Plot_path, Logspace, Four_Plot)
            rxn_interest_plots(Flame_speed_filter, Rxns, Plot_path, Logspace)
    for Trxns in Topnrxns:
        rxn_strength_plots(Flame_speed_filter, Trxns, Nrxns,
                           Threshold, Plot_path, Logspace, Four_Plot)
        rxn_interest_plots(Flame_speed_filter, Trxns, Plot_path, Logspace)


if __name__ == "__main__":
    # Open and parse the input file
    inputs = list(yaml.safe_load_all(open('input.yaml', 'r')))
    main(**inputs[1])