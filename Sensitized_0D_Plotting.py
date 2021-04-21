# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 15:46:05 2020

@author: Kodo Bear
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpl_patches
dirname = os.path.normpath(os.path.dirname(__file__))
plt.style.use(os.path.join(dirname, 'CSULA_Combustion_test.mplstyle'))
#from tkinter import filedialog


def rxn_plots(max_rxns, species):
    for specie in species[0]:
        Temperature  = []
        Pressure     = []
        Equivalence  = []
        Fuel         = []
        Rxns         = []
        num_ticks = 10
        for ents in max_rxns:
            if ents[1] is None:
                continue
            else:
                Rxns.append(ents[1][0])
            # entry = [k for k in ents[0].values()]
            # a     = species[2][0]
            # phi   = a*(entry[0]/entry[1])
            # Temperature.append(ents[1])
            # Pressure.append(ents[2])
            # Equivalence.append(ents[3])
            # Fuel.append(entry[0])
            # Rxns.append(rxn)
            Temperature.append(ents[0]['temperature'])
            Pressure.append(ents[0]['pressure'])
            Equivalence.append(ents[0]['phi'])
            Fuel.append(ents[0]['fuel'])


        # fs = 15 #Controls font size
        varnames = ['P', 'phi', 'X', 'T']
        # info is a dictionary with an entry for each figure. The entries are
        # [x label, x points, name]
        info = {'P': ['Pressure [kPa]', Pressure, 'Pressure'],
                'phi': ['Equivalence Ratio [$\phi$]', Equivalence,
                        'Equivalence'],
                'X': ['Fuel [mole fraction]', Fuel, 'Fuel'],
                'T': ['Temperature [K]', Temperature, 'Temperature']}
        # fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))
        fig, axes = plt.subplots(nrows=2, ncols=2)
        #Four subplots of rxn versus independant variables
        for var, ax, string in zip(varnames, axes.flatten(),
                                   ['(a)', '(b)', '(c)', '(d)']):
            if Array_Type == 'log':
                ax.set_xscale('log')
            # ax.plot(info[var][1], Rxns, ls='none',
            #         marker='o', mfc='none', mec='k')
            ax.plot(info[var][1], Rxns)
            ax.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
            # ax.set_xlabel(info[var][0], fontsize=fs)
            ax.set_xlabel(info[var][0])
            # ax.set_ylabel('Reaction #', fontsize=fs)
            ax.set_ylabel('Reaction #')
            ax.set_ylim([0, int(max(Rxns))+1])
            if max(Rxns) <= num_ticks:
                ax.set_yticks(np.arange(0,max(Rxns)+1, step=2))
            else:
                ax.set_yticks(np.arange(0,max(Rxns),
                                        step=int(max(Rxns)/num_ticks)))
            ax.grid(True)

        # fig.suptitle('Time Instance: '+format(species[3][0])+
        #              ' s    Species: '+ specie, fontsize=12)
        fig.suptitle('Time Instance: '+format(species[3][0])+
                     ' s    Species: '+ specie)
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.savefig(Plot_path+'\\['+format(specie)+
                                  '] Most Sensitive Reactions.png')
        plt.close(fig)

        #Individual plots of Rxn versus independant variable
        for var in varnames:
            # plt.figure("Sensitivity of "+specie+" to "+info[var][2],
            #            figsize=(15,10))
            plt.figure("Sensitivity of "+specie+" to "+info[var][2])
            if Array_Type == 'log':
                plt.xscale('log')
            # plt.plot(info[var][1], Rxns, ls='none',
            #          marker='o', mfc='none', mec='k')
            plt.plot(info[var][1], Rxns)
            # plt.xlabel(info[var][0], fontsize=fs)
            plt.xlabel(info[var][0])
            # plt.ylabel('Most Sensitive Reaction Number for '+specie,
            #            fontsize=10)
            plt.ylabel('Most Sensitive Reaction Number for '+specie)
            plt.ylim([0,int(max(Rxns))+1])
            if max(Rxns) <= 25:
                plt.yticks(np.arange(0,max(Rxns)+1, step=1))
            else:
                plt.yticks(np.arange(0,max(Rxns), step=int(max(Rxns)/20)))
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white",
                                             ec="white", lw=0, alpha=0)] * 2
            labels = ['Time: '+format(species[3][0])+' s','Species: '+specie]
            plt.legend(handles, labels, loc='best', fontsize='small',
                       fancybox=True, framealpha=0.7,
                       handlelength=0, handletextpad=0)
            plt.grid(True)
            # plt.tight_layout()
            plt.savefig(Plot_path+'\\['+format(specie)+
                        '] '+info[var][2]+' Most Sensitive Reaction.png')
            plt.close()


def integrated_strength_plots(int_strength, species, Plot_path, rxn_num):
    varnames = ['P', 'phi', 'X', 'T']
    # info is a dictionary with an entry for each figure. The entries are
    # [x label, key for int_strength]
    info = {'P': ['Pressure [kPa]', 'pressure'],
            'phi': ['Equivalence Ratio', 'phi'],
            'X': ['Fuel Mole Fraction', 'fuel'],
            'T': ['Temperature [K]', 'temperature']}
    file_format = '[{:}] Normalized Integrated Strength.png'.format
    for i in range(len(species[0])):
        spec = species[0][i]
        # fig, axes = plt.subplots(
        #     nrows=2, ncols=2, sharey=True, figsize=(15, 10),
        #     subplot_kw={'ylabel': r'$\hat{IS}_{' + spec + ', j}$'})
        fig, axes = plt.subplots(nrows=2, ncols=2, sharey=True,
            subplot_kw={'ylabel': r'$\hat{IS}_{' + spec + ', j}$'})
        for var, ax, string in zip(varnames, axes.flatten(),
                                   ['(a)', '(b)', '(c)', '(d)']):
            # ax.set_xlabel(info[var][0], fontsize=15)
            ax.set_xlabel(info[var][0])
            if Array_Type == 'log':
                ax.set_xscale('log')
            ax.grid(True)
            x = []
            y = []
            key = info[var][1]
            for condition in int_strength:
                if condition[1][0] is None:
                    continue
                else:
                    x.append(condition[0][key])
                    norm_IS = condition[1][i]
                    y.append(norm_IS[rxn_num])
            # ax.plot(x, y, ls='none', marker='o', mfc='none', mec='k')
            ax.plot(x, y)
            # if var != 'T':
            #     by, ty = plt.ylim()
            #     y_perc = (ty-by)*0.01
            #     plt.ylim([by-y_perc, ty+y_perc])
            #     bx, tx = plt.xlim()
            #     x_perc = (tx-bx)*0.01
            #     plt.xlim([bx-x_perc, tx+x_perc])
            ax.axhline(0, c='gray', ls='--', marker="None")
            ax.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        fig.suptitle(species[2][0][rxn_num])
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fname = file_format(spec)
        fig.savefig(os.path.join(Plot_path, fname))
        plt.close(fig)


def integrated_strength_rxn_plots(int_strength, species, Plot_path,
                                  rxn_num, threshold):
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
    # fs = 15
    for i in range(len(species[0])):
        spec = species[0][i]
        # fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
        fig, axes = plt.subplots(nrows=2, ncols=3)
        for condition, ax, string in zip(conditions, axes.flatten(), anlabel):
            # ax.set_xlabel(info[condition[0]][0], fontsize=fs)
            ax.set_xlabel(info[condition[0]][0])
            # ax.set_ylabel(info[condition[1]][0], fontsize=fs)
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
            # ax.plot(x, y, ls='none', marker='o', mfc='none', mec='k')
            ax.plot(x, y)
            ax.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
        fig.suptitle(species[2][0][rxn_num]+
                     '    Threshold >='+format(threshold))
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fname = file_format(spec)
        fig.savefig(os.path.join(Plot_path, fname))
        plt.close(fig)


def reaction_of_interest_plot(max_rxns, species, rxn_number):
    for specie in species[0]:
        Temperature  = []
        Pressure     = []
        Equivalence  = []
        Fuel         = []
        for ents in max_rxns:
            if ents[1] is None:
                continue
            elif ents[1][0] == rxn_number:
                # entry = [k for k in ents[0].values()]
                # a     = species[2][0]
                # phi   = a*(entry[0]/entry[1])
                # Temperature.append(ents[1])
                # Pressure.append(ents[2])
                # Equivalence.append(ents[3])
                # Fuel.append(entry[0])
                Temperature.append(ents[0]['temperature'])
                Pressure.append(ents[0]['pressure'])
                Equivalence.append(ents[0]['phi'])
                Fuel.append(ents[0]['fuel'])

        # fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
        fig, axes = plt.subplots(nrows=2, ncols=3)
        ax = axes.flatten()
        # fs = 15 #Controls font size

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
            # a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
            #        marker='o', mfc='none', mec='k')
            a.plot(cond_dict[x_key][0], cond_dict[y_key][0])
            a.annotate(xy=(0.90,0.90), text=string, xycoords='axes fraction')
            # a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
            a.set_xlabel(cond_dict[x_key][1])
            # a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
            a.set_ylabel(cond_dict[y_key][1])
            if Array_Type == 'log':
                a.set_xscale('log')
                a.set_yscale('log')
            a.grid(True)

        fig.suptitle('Reaction Number '+format(rxn_number)+
                      ' '+species[2][0][rxn_number]+
                      '\nTime Instance: '+format(species[3][0])+' s')
        # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(Plot_path+'\\['+specie+'] Reaction Number '
                    +format(rxn_number)+'.png')
        plt.close(fig)

        #Individual plots for unique paired independant variables
        for cond in conditions:
            x_key = cond[0]
            y_key = cond[1]
            # plt.figure(figsize=(15,10))
            plt.figure()
            # plt.plot(cond_dict[x_key][0], cond_dict[y_key][0],
            #          ls='none', marker='o', mfc='none', mec='k')
            plt.plot(cond_dict[x_key][0], cond_dict[y_key][0])
            # plt.xlabel(cond_dict[x_key][1], fontsize=fs)
            plt.xlabel(cond_dict[x_key][1])
            # plt.ylabel(cond_dict[y_key][1], fontsize=fs)
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
            # plt.tight_layout()
            plt.savefig(Plot_path+'\\['+specie+'] Reaction Number '
                        +format(rxn_number)+' '+cond_dict[y_key][2]+
                        ' vs '+cond_dict[x_key][2]+'.png')
            plt.close()


def create_csv(int_strength, out_path, rxn_num, spec=0):
    """ Save the integrated strength in a sorted csv file"""
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


if __name__ == "__main__":
    #TODO: use same file opening method as 1D plotting, where it can use the
    # last folder *simulated*, too.
    Folder_name = input('Please type name of folder. if blank, use last folder:\n')
    if Folder_name == '':
        with open('last run 0d.pkl', 'rb') as f:
            Folder_name = pickle.load(f)
        print('Loading ' + Folder_name)
    Load_folder = '\\'+Folder_name
    with open('last run 0d.pkl', 'wb') as f:
        pickle.dump(Folder_name, f)

    Load_path   = '0D_Sensitivity_Results'+Load_folder
    Plot_path   = Load_path+'\\ZeroD_Sensitivity_Plots'

    file = open(Load_path+'\\Case Description.txt','r')
    Description = file.read()
    file.close()
    print(Description)

    Species_Rxn  = pickle.load(open(Load_path+'\\Species_Rxn.pkl', 'rb'))
    Max_Sens_Rxn = pickle.load(open(Load_path+'\\Max_Sens_Rxn.pkl', 'rb'))
    Case_Params  = pickle.load(open(Load_path+'\\Case_Parameters.pkl', 'rb'))

    Threshold = 2
    Reaction_of_Interest = [15]
    Array_Type = Species_Rxn[4][0]

    with open(os.path.join(Load_path, 'Integrated Strength.pkl'), 'rb') as f:
        Int_Strength = pickle.load(f)

    rxn_plots(Max_Sens_Rxn, Species_Rxn)
    for roi in Reaction_of_Interest:
        reaction_of_interest_plot(Max_Sens_Rxn, Species_Rxn, roi)
        integrated_strength_plots(Int_Strength, Species_Rxn, Plot_path, roi)
        integrated_strength_rxn_plots(Int_Strength, Species_Rxn,
                                      Plot_path, roi, Threshold)
