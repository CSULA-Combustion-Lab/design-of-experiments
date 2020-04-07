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

def rxn_plots(F_info, save_path):
    """[Insert Information]"""
    Tint     = F_info[0]['Conditions'][3]
    Pressure = []
    Fuel     = []
    Phi      = []
    Max_rxn  = []
    for f in F_info:
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
                    +format(Tint)+' [K].png')
    plt.close(fig)
    
    
def rxn_strength_plots(F_info, nrxns):
    """[Insert Information]"""
    flame_strength = []
    for flame in F_info:
        flame_sens = []
        for fs in nrxns:
            flame_sens.append()
        flame_info = {'Sensitivity': flame_sens,
                      'Conditions': flame['Conditions']}
        flame_strength.append(flame_info)
    
def rxn_interest_plots(F_info, rxn_int, save_path):
    """[Insert Information]"""
    Pressure = []
    Fuel     = []
    Phi      = []
    Tint     = F_info[0]['Conditions'][3]
    Rxn_name = F_info[0]['Flame'][0][rxn_int][2]
    Rxn_num  = F_info[0]['Flame'][0][rxn_int][0]
    for f in F_info:
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 10))
    ax        = axes.flatten()
    fs        = 15
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel [Mole Fraction]'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'T': [Tint, 'Temperature [K]']}
    conditions = [('P', 'F'), ('P', 'Phi'), ('F', 'Phi')]
    #Three subplots of unique paired independent variables
    for a, condition in zip(ax, conditions):
        x_key = condition[0]
        y_key = condition[1]
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        a.set_xscale('log')
        a.set_yscale('log')
        a.grid(True)

        fig.suptitle('Reaction Number: '+format(Rxn_num)+
                      'Reaction Name: '+Rxn_name+
                      '\nInitial Temperature: '+format(Tint)+' [K]')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Reaction Number '+format(Rxn_num)+
                    ' Initial Temperature '+format(Tint)+' [K].png')
    plt.close(fig)
    
    
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

    Load_path   = 'Flame_Sensitivity_Results'+Load_folder
    Plot_path   = Load_path+'\\Flame_Sensitivity_Plots'

    file = open(Load_path+'\\Case Description.txt','r')
    Description = file.read()
    file.close()
    print(Description)

    with open(os.path.join(Load_path, 'Flame Information.pkl'), 'rb') as f:
        flame_info = pickle.load(f)
        
    flame = []
    no_flame = []
    for x in flame_info:
        if x['Flame'] == None:
            no_flame.append(x)
        else:
            flame.append(x)
            
    #Note flame is a dictionary of {'Flame': [Flame_sens Su], 
    #                               'Conditions': [P, Fuel, Phi, Tin, Mix]}
    
    #Plot Functions
    Rxn_Interest = 5 #Reaction number of the reaction of interest
    rxn_plots(flame, Plot_path)