# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:20:01 2020

@author: Kodo Bear
"""

import os
import numpy
import pickle
import matplotlib.pyplot as plt
#from tkinter import filedialog

def rxn_plots(f_info, save_path, log):
    """[Insert Information]"""
    Tint     = f_info[0]['Conditions'][3]
    Pressure = []
    Fuel     = []
    Phi      = []
    Oxygen   = []
    Max_rxn  = []
    for f in f_info:
        max_rxn  = max_sens(f)
        max_num  = max_rxn[0]
        Max_rxn.append(max_num)
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Oxygen.append(f['Conditions'][4][1][1])
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10), sharey=True)
    ax        = axes.flatten()
    fs        = 15
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel Mole Fraction'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'O2': [Oxygen, 'Oxygen Mole Fraction'],
                 'T': [Tint, 'Temperature [K]'],
                 'Rxn': [Max_rxn, 'Rxn Number']}
    conditions = ['P', 'F', 'Phi', 'O2']
    #Three subplots of max rxn number against independent variables
    for a, condition in zip(ax, conditions):
        x_key = condition
        y_key = 'Rxn'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        if log:
            a.set_xscale('log')
        elif x_key == 'P':
            a.set_xscale('log')
        a.grid(True)

        fig.suptitle('Initial Temperature: '+format(Tint)+' [K]')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Max Reactions with Initial Temperature '
                    +format(Tint)+'K.png')
    plt.close(fig)
    
    
def rxn_strength_plots(f_info, rxn_int, nrxns, threshold, save_path, log):
    """[Insert Information]"""
    Tint = f_info[0]['Conditions'][3]
    #List for strengths
    Pressure      = []
    Fuel          = []
    Phi           = []
    Oxygen        = []
    Sens_strength = []
    #Lists for above threshold
    P_threshold        = []
    F_threshold        = []
    Phi_threshold      = []
    Oxygen_threshold   = []
    Sens_str_threshold = []
    for f in f_info:
        average_nrxns = rxn_average(f, nrxns)
        strength = f['Flame'][0][rxn_int][1]/average_nrxns
        Sens_strength.append(strength)
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Oxygen.append(f['Conditions'][4][1][1])
        
    for n in range(len(Sens_strength)):
        if abs(Sens_strength[n]) >= threshold:
            P_threshold.append(f_info[n]['Conditions'][0])
            F_threshold.append(f_info[n]['Conditions'][1])
            Phi_threshold.append(f_info[n]['Conditions'][2])
            Oxygen_threshold.append(f_info[n]['Conditions'][4][1][1])
            Sens_str_threshold.append(Sens_strength[n])
            
    #Dictionary organized as follow 'Key': [Data, axis-label, Data_threshold]
    cond_dict = {'P': [Pressure, 'Pressure [atm]', P_threshold],
                 'F': [Fuel, 'Fuel Mole Fraction', F_threshold],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]', Phi_threshold],
                 'O2': [Oxygen, 'Oxygen Mole Fraction', Oxygen_threshold],
                 'T': [Tint, 'Temperature [K]'],
                 'Strength': [Sens_strength, 'Rxn Strength',
                              Sens_str_threshold]}
    
    figa, axesa = plt.subplots(nrows=2, ncols=2, figsize=(15, 10), sharey=True)
    axa         = axesa.flatten()
    fs          = 15
    conditions_a = ['P', 'F', 'Phi', 'O2']
    for a, condition in zip(axa, conditions_a):
        x_key = condition
        y_key = 'Strength'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        if log:
            a.set_xscale('log')
        elif x_key == 'P':
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
                
    figb, axesb  = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
    axb          = axesb.flatten()
    conditions_b = [('P', 'F'), ('P', 'Phi'), ('P','O2'),
                    ('F', 'Phi'), ('F', 'O2'), ('O2', 'Phi')]
    for a, condition in zip(axb, conditions_b):
        x_key = condition[0]
        y_key = condition[1]
        a.plot(cond_dict[x_key][2], cond_dict[y_key][2], ls='none',
               marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        if log:
            a.set_xscale('log')
        elif x_key == 'P':
            a.set_xscale('log')
        if log:
            a.set_yscale('log')
        elif y_key == 'P':
            a.set_yscale('log')
        a.grid(True)

        figb.suptitle('Reaction '+format(f_info[0]['Flame'][0][rxn_int][2])+
                      '\n Threshold >= '+format(threshold)+
                      ' Initial Temperature: '+format(Tint)+' [K]')
        figb.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Reaction '
                    +format(f_info[0]['Flame'][0][rxn_int][0])+
                    ' Average Strength above_equal to '+format(threshold)+
                    ' with Initial Temperature '+format(Tint)+'K.png')
    plt.close(figb)

      
def rxn_interest_plots(f_info, rxn_int, save_path, log):
    """[Insert Information]"""
    Pressure = []
    Fuel     = []
    Phi      = []
    Oxygen   = []
    Tint     = f_info[0]['Conditions'][3]
    Rxn_name = f_info[0]['Flame'][0][rxn_int][2]
    Rxn_num  = f_info[0]['Flame'][0][rxn_int][0]
    for f in f_info:
        Rxn_sens = f['Flame'][0][rxn_int][1]
        Max_rxn_sens = max_sens(f)
        if Rxn_sens == Max_rxn_sens[1]:
            Pressure.append(f['Conditions'][0])
            Fuel.append(f['Conditions'][1])
            Phi.append(f['Conditions'][2])
            Oxygen.append(f['Conditions'][4][1][1])
            
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))
    ax        = axes.flatten()
    fs        = 15
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
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
                marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        if log:
            a.set_xscale('log')
        elif x_key == 'P':
            a.set_xscale('log')
        if log:
            a.set_yscale('log')
        elif y_key == 'P':
            a.set_yscale('log')
        a.grid(True)
        fig.suptitle('Reaction Number: '+format(Rxn_num)+
                      ' Reaction Name: '+Rxn_name+
                      '\nInitial Temperature: '+format(Tint)+' [K]')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Reaction Number '+format(Rxn_num)+
                    ' Max Sensitivity for Parameters with'
                    ' Initial Temperature '+format(Tint)+'K.png')
    plt.close(fig)
    

def flame_speed_plots(f_info, save_path, log):
    """[Insert Information]"""
    Pressure = []
    Fuel     = []
    Phi      = []
    Dil_frac = []
    Oxygen   = []
    Su       = []
    Diluent  = f_info[0]['Conditions'][4][0][0]
    Tint     = f_info[0]['Conditions'][3]
    for s in f_info:
        Pressure.append(s['Conditions'][0])
        Fuel.append(s['Conditions'][1])
        Phi.append(s['Conditions'][2])
        Dil_frac.append(s['Conditions'][4][0][1])
        Oxygen.append(s['Conditions'][4][1][1])
        Su.append(s['Flame'][1])
        
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10), sharey=True)
    ax        = axes.flatten()
    fs        = 15
    #Dictionary organized as follow 'Key': [Data, axis-label]
    cond_dict = {'P': [Pressure, 'Pressure [atm]'],
                 'F': [Fuel, 'Fuel Mole Fraction'],
                 'Phi':[Phi, 'Equivalence Ratio [$\phi$]'],
                 'Dil':[Dil_frac, Diluent+' Mole Fraction'],
                 'O2':[Oxygen, 'Oxygen Mole Fraction'],
                 'T': [Tint, 'Temperature [K]'],
                 'Su': [Su, 'Su [m/s]']}
    conditions = ['P', 'F', 'Phi', 'O2']
    for a, condition in zip(ax, conditions):
        x_key = condition
        y_key = 'Su'
        a.plot(cond_dict[x_key][0], cond_dict[y_key][0], ls='none',
                marker='o', mfc='none', mec='k')
        a.set_xlabel(cond_dict[x_key][1], fontsize=fs)
        a.set_ylabel(cond_dict[y_key][1], fontsize=fs)
        if log:
            a.set_xscale('log')
        elif x_key == 'P':
            a.set_xscale('log')
        a.grid(True)
        fig.suptitle('Initial Temperature: '+format(Tint)+' [K]')
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(save_path+'\\Flame Speed vs. Independant Variables with'
                    ' Initial Temperature '+format(Tint)+'K.png')
    plt.close(fig)


def max_rxn_text(f_info, save_path):
    """[Insert Information]"""
    Tint     = f_info[0]['Conditions'][3]
    Pressure = []
    Fuel     = []
    Phi      = []
    Oxygen   = []
    Max_rxn  = []
    for f in f_info:
        max_rxn  = max_sens(f)
        Max_rxn.append([max_rxn[0], max_rxn[2]])
        Pressure.append(f['Conditions'][0])
        Fuel.append(f['Conditions'][1])
        Phi.append(f['Conditions'][2])
        Oxygen.append(f['Conditions'][4][1][1])
    Max_rxn_dict = {}
    Max_rxn_eq   = {}
    for n in range(len(Max_rxn)):
        rxn_num = str(Max_rxn[n][0]+1)
        if rxn_num not in Max_rxn_dict.keys():
            Max_rxn_dict[rxn_num] = [[Max_rxn[n][1], Tint],
                                     [Pressure[n], Fuel[n],
                                      Phi[n], Oxygen[n]]]
            Max_rxn_eq[rxn_num]   = [Max_rxn[n][1], 1, Tint]
        elif rxn_num in Max_rxn_dict.keys():
            Max_rxn_dict[rxn_num].append([Pressure[n], Fuel[n],
                                          Phi[n], Oxygen[n]])
            Max_rxn_eq[rxn_num][1] += 1
        else:
            print('Something went wrong!')
    return Max_rxn_dict, Max_rxn_eq


# def rxn_sens_bar_plots(flame, nrxns, spec_cond, save_path):
#     """[Insert Information]"""
#     fig, ax = plt.subplots()
#     ax.grid(axis='x', which='major', ls='--')
#     ax.grid(axis='y', which='minor', c='k')
#     sens     = []
#     pressure = []
#     oxygen   = []
#     for i in spec_cond['P']:
#         for j in spec_cond['O2']:
#             for f in flame:
#                 if f['Condtions'][0] == i and f['Conditions'][4][1][1] == j:
#                     sens.append(f['Flame'][0])
#                     pressure.append()
#                     oxygen.append()
#     sens.sort(key=lambda x: abs(x[1]), reverse=True)
#     sens_plot = sens[:nrxns]
#     ylocs = numpy.arange(nrxns)
#     ax.barh(ylocs, [x[1] for x in sens_plot], align='center')
#     ax.set_yticks(ylocs)
#     ax.set_yticklabels([x[2] for x in sens_plot])
#     ax.set_yticks(ylocs - 0.5, minor=True)
# #    ax.tick_params(axis='y', which='minor', bottom='off')
#     ax.invert_yaxis()
#     ax.axvline(c='k')
#     ax.set_xlabel('Normalized Sensitivity')
#     ax.set_ylim([max(ylocs)+0.5, min(ylocs)-0.5])
#     fig.tight_layout()
#     fig.savefig(os.path.join(save_path, 'Flame Sensitivity Bar Plot.png'))
#     plt.close(fig)

    
def rxn_average(flame, nrxns):
    """[Insert Information]"""
    rxn_sens = []
    for s in flame['Flame'][0]:
        rxn_sens.append(abs(s[1]))
    rxn_sens.sort()
    top_nrxns    = rxn_sens[-nrxns:]
    average_sens = (sum(top_nrxns)/len(top_nrxns))
    return average_sens


def max_sens(flame):
    max_rxn_num  = 0
    max_rxn_sens = 0
    for m in flame['Flame'][0]:
        if abs(m[1]) > max_rxn_sens:
            max_rxn_num  = m[0]
            max_rxn_sens = m[1]
            max_rxn_eq   = m[2]
    max_rxn = [max_rxn_num, max_rxn_sens, max_rxn_eq]
    return max_rxn

    
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
    #Creates list that only appends information where flame speed is above min     
    Min_speed          = 0 #Minimum flame speed, lower limit
    Flame_speed_filter = []
    for x in Flame:
        if x['Flame'][1] >= Min_speed:
            Flame_speed_filter.append(x)
            
    #Note Flame and No_flame are dictionaries 
    # {'Flame': [Flame_sens Su Flame_rho, Flame_temp, Mingrid, Mult_sort], 
    #  'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
    
    #Plot Functions
    Rxn_interest = numpy.arange(69,81) #Reaction number of the reaction of interest
    Nrxns        = 5 #Top n-reactions
    Threshold    = 0.5 #Threshold for rxn_interst to be above in average strength
    Logspace     = True #If true all plots will use logspace
    Spec_Conditions = {'Key': ['P', 'O2'],
                       'O2': [0.25, 0.5],
                       'P': [0.5, 1]}
    Max_rxn_cond, Max_rxn_count = max_rxn_text(Flame_speed_filter, Plot_path)  
    rxn_plots(Flame_speed_filter, Plot_path, Logspace)
    flame_speed_plots(Flame_speed_filter, Plot_path, Logspace)
    for Rxns in Rxn_interest:
        rxn_strength_plots(Flame_speed_filter, Rxns, Nrxns,
                            Threshold, Plot_path, Logspace)
        rxn_interest_plots(Flame_speed_filter, Rxns, Plot_path, Logspace)
    # rxn_sens_bar_plots(Flame_speed_filter, Nrxns, Spec_Conditions, Plot_path)