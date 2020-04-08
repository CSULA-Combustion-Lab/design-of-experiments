# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:20:01 2020

@author: Kodo Bear
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpl_patches
#from tkinter import filedialog

def rxn_plots(f_info, save_path):
    """[Insert Information]"""
    Tint     = f_info[0]['Conditions'][3]
    Pressure = []
    Fuel     = []
    Phi      = []
    Max_rxn  = []
    for f in f_info:
        max_sens = 0
        max_num  = 0
        for m in f['Flame'][0]:
            ms = abs(m[1])
            if ms > max_sens:
                max_num = m[0]
        Max_rxn.append(f['Flame'][0][max_num][1])
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 10), sharey=True)
    ax        = axes.flatten()
    fs        = 15
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel [Mole Fraction]'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'T': [Tint, 'Temperature [K]'],
                 'Rxn': [Max_rxn, 'Rxn Number']}
    conditions = ['P', 'F', 'Phi']
    #Three subplots of max rxn number against independent variables
    for a, condition in zip(ax, conditions):
        x_key = condition
        y_key = 'Rxn'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        a.set_xscale('log')
        a.grid(True)

        fig.suptitle('Initial Temperature: '+format(Tint)+' [K]')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Max Reactions with Initial Temperature '
                    +format(Tint)+'K.png')
    plt.close(fig)
    
    
def rxn_strength_plots(f_info, rxn_int, nrxns, threshold, save_path):
    """[Insert Information]"""
    Tint     = f_info[0]['Conditions'][3]
    Pressure = []
    Fuel     = []
    Phi      = []
    sens_strength = []
    for f in f_info:
        average_nrxns = rxn_average(f, nrxns)
        strength = f['Flame'][0][rxn_int][1]/average_nrxns
        sens_strength.append(strength)
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
    figa, axesa = plt.subplots(nrows=1, ncols=3, figsize=(15, 10), sharey=True)
    axa         = axesa.flatten()
    fs          = 15
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel [Mole Fraction]'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'T': [Tint, 'Temperature [K]'],
                 'Strength': [sens_strength, 'Rxn Strength']}
    conditions_a = ['P', 'F', 'Phi']
    for a, condition in zip(axa, conditions_a):
        x_key = condition
        y_key = 'Strength'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        a.set_xscale('log')
        a.grid(True)

        figa.suptitle('Reaction '+format(f_info[0]['Flame'][0][rxn_int][2])+
                      '\nInitial Temperature: '+format(Tint)+' [K]')
        figa.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Reaction '
                    +format(f_info[0]['Flame'][0][rxn_int][0])+
                    ' Average Strength against top '+format(nrxns)+' rxns'
                    ' with Initial Temperature '+format(Tint)+'K.png')
    plt.close(figa)
    
    figb, axesb = plt.subplots(nrows=1, ncols=3, figsize=(15, 10))
    axb         = axesb.flatten()
    conditions_b = [('P', 'F'), ('P', 'Phi'), ('F', 'Phi')]
    for a, condition in zip(axb, conditions_b):
        x_key = condition[0]
        y_key = condition[1]
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        a.set_xscale('log')
        a.grid(True)

        figa.suptitle('Reaction '+format(f_info[0]['Flame'][0][rxn_int][2])+
                      '\nInitial Temperature: '+format(Tint)+' [K]')
        figa.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Reaction '
                    +format(f_info[0]['Flame'][0][rxn_int][0])+
                    ' Average Strength against top '+format(nrxns)+' rxns'
                    ' with Initial Temperature '+format(Tint)+'K.png')
    plt.close(figa)


def rxn_average(flame, nrxns):
    """[Insert Information]"""
    rxn_sens     = []
    for s in flame['Flame'][0]:
        rxn_sens.append(abs(s[1]))
    rxn_sens.sort()
    top_nrxns = rxn_sens[-nrxns:]
    average_sens = (sum(top_nrxns)/len(top_nrxns))
    return average_sens
 
    
# def rxn_interest_plots(f_info, rxn_int, save_path):
#     """[Insert Information]"""
#     Pressure = []
#     Fuel     = []
#     Phi      = []
#     Tint     = f_info[0]['Conditions'][3]
#     Rxn_name = f_info[0]['Flame'][0][rxn_int][2]
#     Rxn_num  = f_info[0]['Flame'][0][rxn_int][0]
#     for f in f_info:
#         Pressure.append(f['Conditions'][0])
#         Fuel.append(f['Conditions'][1])
#         Phi.append(f['Conditions'][2])
#     fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 10))
#     ax        = axes.flatten()
#     fs        = 15
#     #Dictionary organized as follow 'Key': [Data, axis-label]
#     cond_dict = {'P': [Pressure, 'Pressure [atm]'],
#                  'F': [Fuel, 'Fuel [Mole Fraction]'],
#                  'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
#                  'T': [Tint, 'Temperature [K]']}
#     conditions = [('P', 'F'), ('P', 'Phi'), ('F', 'Phi')]
#     #Three subplots of unique paired independent variables
#     for a, condition in zip(ax, conditions):
#         x_key = condition[0]
#         y_key = condition[1]
#         a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
#                marker='o', mfc='none', mec='k')
#         a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
#         a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
#         a.set_xscale('log')
#         a.set_yscale('log')
#         a.grid(True)

#         fig.suptitle('Reaction Number: '+format(Rxn_num)+
#                       'Reaction Name: '+Rxn_name+
#                       '\nInitial Temperature: '+format(Tint)+' [K]')
#         fig.tight_layout(rect=[0, 0.03, 1, 0.95])
#         plt.savefig(save_path+'\\Reaction Number '+format(Rxn_num)+
#                     ' Initial Temperature '+format(Tint)+'K.png')
#     plt.close(fig)

    
if __name__ == "__main__":
    Folder_name = input('Please type name of folder.'
                        '\n If blank, use last folder:\n')
    if Folder_name == '':
        with open('last run.pkl', 'rb') as f:
            Folder_name = pickle.load(f)
        print('Loading ' + Folder_name)
    Load_folder = '\\'+Folder_name
    with open('last run.pkl', 'wb') as f:
        pickle.dump(Folder_name, f)

    #Paths for loading and saving files
    Load_path = 'Flame_Sensitivity_Results'+Load_folder
    Plot_path = Load_path+'\\Flame_Sensitivity_Plots'

    #Open up text file with description of simulation
    File = open(Load_path+'\\Case Description.txt','r')
    Description = File.read()
    File.close()
    print(Description)

    #Import flame file found in corresponding folder
    with open(os.path.join(Load_path, 'Flame Information.pkl'), 'rb') as f:
        Flame_info = pickle.load(f)
        
    #Create two lists of flame and no_flame created from flame_info    
    Flame = []
    No_flame = []
    for x in Flame_info:
        if x['Flame'] == None:
            No_flame.append(x)
        else:
            Flame.append(x)
            
    #Note flame is a dictionary of {'Flame': [Flame_sens Su], 
    #                               'Conditions': [P, Fuel, Phi, Tin, Mix]}
    
    #Plot Functions
    Rxn_interest = 5 #Reaction number of the reaction of interest
    Nrxns        = 5 #Top n-reactions
    Threshold    = 2 #Threshold for rxn_interst to be above in average strength
    rxn_plots(Flame, Plot_path)
    rxn_strength_plots(Flame, Rxn_interest, Nrxns, Threshold, Plot_path)