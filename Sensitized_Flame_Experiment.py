# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:52:04 2020

@author: Kodo Bear
"""

import os
import copy
import errno
import pickle
import time, sys
import numpy as np
import cantera as ct
import itertools as it
path = (r'D:\School\Cal State La\Extracurricular\NOx Combustion Research')
sys.path.insert(0, path)
from utilities import flame
from datetime import datetime
ct.suppress_thermo_warnings() #Suppress cantera warnings!

####Set experiment parameters
mechanism = 'Li_Methylnitrite.cti' #Mechanism file
gas = ct.Solution(mechanism)

#Working directory
flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

#Parameters for main loop
P    = np.logspace(np.log10(0.5), np.log10(2), 4) #Pressure [atm]
Phi  = np.logspace(np.log10(0.25), np.log10(4), 12) #Equivalence ratio
Fuel = np.logspace(0.06, 0.85, 1) #Fuel mole fraction
OtO  = np.logspace(np.log10(0.1), np.log10(0.95), 12) #Oxygen to Oxidizer ratio [Air = .21]


#Initial Temperature
Tint = 323 #Temperature [K]

#Parameters for mixture
fuel_name = 'CH3ONO' #chemical formula of fuel

#Fuel C(x)H(y)O(z)
fuel_index = gas.species(gas.species_index(fuel_name)).composition
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
    
multifuel = False # if true, additional fuels will be added from fuels list
if multifuel:
    # fuel_list should have the fuel name followed by the percentage of the
    # named fuel relative to the total fuel. Fuel percentages should be less
    # than 100%
    fuel_list = ['CH3OH', .50]
    percent_fuels = 0
    for n in range(1,len(fuel_list),2):
        percent_fuels += fuel_list[n]
    if percent_fuels >= 1:
        sys.exit("Error! Multifuel percentage above/equal to 100%")
    mf_percent = 1 - percent_fuels #Main fuel percentage
    x *= mf_percent
    y *= mf_percent
    z *= mf_percent
    for n in range(0,len(fuel_list),2):
        fuel_index = gas.species(gas.species_index(fuel_list[n])).composition
        if 'C' in fuel_index:
            x += fuel_index['C']*fuel_list[n+1]
        if 'H' in fuel_index:
            y += fuel_index['H']*fuel_list[n+1]
        if 'O' in fuel_index:
            z += fuel_index['O']*fuel_list[n+1]

a = x+y/4-z/2       #molar oxygen-fuel ratio
diluent_name = 'N2' #chemical formula of diluent

custom     = False # If true, custom styles used for range and save files
oxidizer   = True # If true, OtO will be used instead of Fuel Mole Fraction
debug      = False # If true, print lots of information for debugging.
save_files = True # If true, save files for plotting script

#Debug Files
DEBUG_FMT = 'Removing condition: T={:.0f}, P={:.0f}, phi={:.3g}, fuel={:.3g}'
#Debug parameters [Pressure, Equivalence Ratio, Fuel or Oxygen, Temperature]
Debug_params = [1, 1, 0.21, 300]
LogLevel     = 1

#Flame Conditions
mingrid = 200
mul_soret = False

#Custom loop
if custom:
    P      = [0.5, 1]
    Fuel   = 0.05
    Oxygen = np.linspace(0.05,0.95,50)
    conditions = {'Parameters': [P, Fuel, Oxygen, Tint],
                  'Mixture': [fuel_name, x, y, z, a, diluent_name],
                  'Flame': [mingrid, mul_soret],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxidizer, custom, multifuel]}
#Normal loop
else:
    conditions = {'Parameters': [P, Phi, Fuel, Tint, OtO],
                  'Mixture': [fuel_name, x, y, z, a, diluent_name],
                  'Flame': [mingrid, mul_soret],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxidizer, custom, multifuel]}

if multifuel:
    conditions.update([('MultiFuels', fuel_list)])

def parallelize(param, cond, fun):
    """[Fill in information]"""
    from multiprocessing import cpu_count, Pool
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
        results.append(pool.apply_async(fun, args=(*x, cond['Parameters'][3],
                                                   cond)))
    pool.close()
    pool.join()

    # Get the results
    datadict = dict()
    casenum  = 0
    for p in results:
        try:
            # Assign it to datadict. This is ordered by the time when each
            # simulation starts, not when they end
            datadict[casenum] = p.get()
        except RuntimeError: # I'm not sure what this is
            print('\nUnknown RunTimeError.')
            datadict[casenum] = None
        casenum += 1
    outlist = [datadict[k] for k in datadict.keys()] # Convert back to list
    return outlist


def flame_sens(P, Phi, F_O, Tin, Cond):
    """[Fill in information]"""
    chem         = Cond['Files'][0]
    tempfile     = Cond['Files'][1]
    a            = Cond['Mixture'][4]
    Fuel_name    = Cond['Mixture'][0]
    Diluent_name = Cond['Mixture'][5]
    debug        = Cond['T/F'][0]
    oxidizer     = Cond['T/F'][1]
    custom       = Cond['T/F'][2]
    multifuel    = Cond['T/F'][3]
    mg           = Cond['Flame'][0]
    ms           = Cond['Flame'][1]
    
    if custom:
        Fuel    = Phi
        Oxygen  = F_O
        Diluent = 1 - Oxygen - Fuel
        OtO     = Oxygen/(Oxygen + Diluent)
        Phi     = a/Oxygen*Fuel
        Mix     = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    elif oxidizer:
        OtO     = F_O
        noxy    = (Phi*OtO)/a
        Oxygen  = OtO/(noxy + 1)
        Fuel    = noxy/(noxy + 1)
        Diluent = (1 - OtO)/(noxy + 1)
        if multifuel:
            mf  = Cond['MultiFuels']
            fuel_percent = 0
            for n in range(1,len(mf),2):
                fuel_percent += mf[n]
            Fuel_main = (1 - fuel_percent)*Fuel
            Mix       = [[Diluent_name, Diluent], ['O2', Oxygen],
                         [Fuel_name, Fuel_main]]
            Mix_check = (Diluent_name+':'+format(Diluent)+', O2:'+
                         format(Oxygen)+', '+Fuel_name+':'+format(Fuel_main))
            for n in range(0,len(mf),2):
                Mix.append([mf[n], mf[n+1]*Fuel])
                Mix_check += ', '+mf[n]+':'+format(mf[n+1]*Fuel)
            gas.X = Mix_check
            Phi_check = gas.get_equivalence_ratio()
        else:
            Mix = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    else:
        Fuel    = F_O
        Oxygen  = a/Phi*Fuel
        Diluent = 1 - Oxygen - Fuel
        OtO     = Oxygen/(Oxygen + Diluent)
        if multifuel:
            mf  = Cond['MultiFuels']
            fuel_percent = 0
            for n in range(1,len(mf),2):
                fuel_percent += mf[n]
            Fuel_main = (1 - fuel_percent)*Fuel
            Mix       = [[Diluent_name, Diluent], ['O2', Oxygen],
                         [Fuel_name, Fuel_main]]
            Mix_check = (Diluent_name+':'+format(Diluent)+', O2:'+
                         format(Oxygen)+', '+Fuel_name+':'+format(Fuel_main))
            for n in range(0,len(mf),2):
                Mix.append([mf[n], mf[n+1]])
                Mix_check += ', '+mf[n]+':'+format(mf[n+1]*Fuel)
            gas.X = Mix_check
            Phi_check = gas.get_equivalence_ratio()
        else:
            Mix = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    
    if multifuel:
        if not np.isclose(Phi, Phi_check):
            flame_info = {'Flame': None,
                          'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
            return flame_info
    
    if Diluent < 0:
        flame_info = {'Flame': None,
                      'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
        return flame_info
   
    if debug:
        logl = Cond['Debug'][1]
        print('\nMixture Composition:')
        print(Mix)     
    else:
        logl = 0
        
    f = flame.Flame(Mix, P, Tin, tempfile, chemfile=chem)
    f.run(mingrid=mg, loglevel=logl, mult_soret=ms)
    if f.flame_result is None:
        flame_info = {'Flame': None,
                      'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
        return flame_info
    f.sensitivity()
    flame_sens = f.sens #Rxn sensitivities 
    Su         = f.flame_result.u[0]  #Flame speed at the front
    flame_T    = f.flame_result.T[-1] #Flame temperature at end
    flame_rho  = f.flame_result.density_mass[0] #Flame density at the front
    flame_info = {'Flame': [flame_sens, Su, flame_rho, flame_T, mg, ms],
                  'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
    return flame_info


def duplicate_reactions(gas):
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


def flame_info_filter(flame_information, duplicate_reactions):
    for f in flame_information:
        if f['Flame'] == None:
            continue
        else:
            for d in duplicate_rxns:
                sum = 0
                for n in duplicate_rxns[d]:
                     sum += f['Flame'][0][n][1]
                     f['Flame'][0][n][1] = 0
                f['Flame'][0][duplicate_rxns[d][0]][1] = sum


if __name__ == "__main__":
    #Start time
    tic = time.time()
    #Initializing
    if custom: 
        totaliterations = len(P)*len(Oxygen)
        paramlist       = []
        for i in P:
            for k in Oxygen:
                paramlist.append((i, Fuel, k))
    elif oxidizer:
        totaliterations = len(P)*len(Phi)*len(OtO)
        paramlist       = list(it.product(P,Phi,OtO))
    else:
        totaliterations = len(P)*len(Phi)*len(Fuel)
        paramlist       = list(it.product(P,Phi,Fuel))

##############################Debug loop#######################################
    if debug:
        print('Debugging in process..')
        print('\nConditions Used:'
              '\nPressure: '+format(Debug_params[0])+' [atm]'
              '\nPhi: '+format(Debug_params[1])+
              '\nFuel: '+format(Debug_params[2])+
              '\nTemperature: '+format(Debug_params[3])+' [K]')
        flame_info_debug = flame_sens(*Debug_params, conditions)
        print('Debuggin complete!')
        toc      = time.time()
        duration = toc - tic
        print('Dubugging time: '+format(duration, '0.5f')+' seconds\n')
   
############################Simulation loop####################################
    else:
        print('Initial number of cases: '+format(len(paramlist)))
        print('\nStart of simulations...')
        sim_start  = time.time()
        flame_info = parallelize(paramlist, conditions, flame_sens)
        sim_end    = time.time()
        sim_time   = sim_end - sim_start
        print('End of simulations')
        print('Simulations took '+format(sim_time, '0.5f')+' seconds.')
        print('\nStart flame information filtering...')
        filter_start = time.time()
        converged    = []
        for x in flame_info:
            if x['Flame'] == None:
                continue
            else:
                converged.append(1)
        flame_info_unfiltered = copy.deepcopy(flame_info)
        duplicate_rxns = duplicate_reactions(gas)
        flame_info_filter(flame_info, duplicate_rxns)
        filter_end  = time.time()
        filter_time = filter_end - filter_start
        print('End of filtering')
        print('Filtering took '+format(filter_time, '0.5f')+ ' seconds.')
        print('Number of cases converged:', len(converged))
        toc      = time.time()
        duration = toc-tic
        print('Total time '+format(duration, '0.5f')+' seconds.\n')
        
#################################Save Files####################################   
    if save_files:
        #Save Path/Parent Directory
        parent_dir = 'Flame_Sensitivity_Results'
        try:
            os.makedirs(parent_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        #######################Debug directory#################################
        if debug:
            print('Start debug file save...')
            debug_dir = 'debug'
            debug_path = os.path.join(parent_dir, debug_dir)
            if not os.path.exists(debug_path):
                try:
                    os.makedirs(debug_path)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
                        
            figures_dir = 'Flame_Sensitivity_Plots'
            save_path   = os.path.join(parent_dir, debug_dir)
            figure_path = os.path.join(save_path, figures_dir)
            if not os.path.exists(figure_path):
                try:
                    os.makedirs(figure_path)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise
            
            print('\nCreating text file...')
            #Text Description
            filename  = 'Case Description.txt'
            filename  = os.path.join(save_path, filename)
            f         = open(filename, "w")
            if oxidizer:
                F_O_text  = "\nOxygen to Oxidizer Fraction: "
            else:
                F_O_text  = "\nFuel Mole Fraction: "
            if multifuel:
                MF_text = ("Multifuel = "+format(multifuel)+"\n"
                           "Fuels\Percentages = "+format((fuel_list))+"\n")
            else:
                MF_text = "Multifuel = "+format(multifuel)+"\n"    
            text_description = ("This file provides debug information.\n The "
                                "following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
                                +MF_text+
                                "==============================\n"
                                "\n==========Parameters==========\n"
                                "Initial Temperature: "+format(Debug_params[3])
                                +" [Kelvin]\nPressure: "
                                +format(Debug_params[0])+" [atm]\n"
                                "Equivalence Ratio: "+format(Debug_params[1])+
                                F_O_text+format(Debug_params[2])+"\n"
                                "==============================\n"
                                "\n==========Flame-Information==========\n"
                                "FLame Speed: "+format(flame_info['Flame'][1])
                                +" [m/s]"
                                "\nMingrid = "+format(mingrid)+
                                "\nMult_Soret = "+format(mul_soret)+
                                "\n=====================================\n"
                                "\n==============Time==============\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "================================")
            f.write(text_description)
            f.close()
            print('\nText file created')
            
            debug_file = os.path.join(debug_path, 'flame_info_debug.pkl')
            with open(debug_file, 'wb') as f:
                pickle.dump(flame_info_debug, f)
            print('End debug file save')
         
        ########################Normal directory###############################
        else:               
            #Create Directory Name
            print('Creating Directory...')
            now = datetime.now()
            dt_string = now.strftime("%Y_%m_%d %H.%M.%S Flame_Speed_Sens")
            directory = dt_string
            save_path = os.path.join(parent_dir, directory)
            os.makedirs(save_path)
    
            figures_dir = 'Flame_Sensitivity_Plots'
            figure_path = os.path.join(save_path, figures_dir)
            os.makedirs(figure_path)
            print('Directory Created')
            
            print('\nCreating text file...')
            #Text Description
            filename  = 'Case Description.txt'
            filename  = os.path.join(save_path, filename)
            f         = open(filename, "w")
            if oxidizer:
                F_O_text = ("\nOxygen to Oxidizer Mole Fraction Range: "
                            +format(OtO)+"\n")
            else:
                F_O_text = "\nFuel Mole Fraction Range: "+format(Fuel)+"\n"
            if multifuel:
                MF_text = ("Multifuel = "+format(multifuel)+"\n"
                           "Fuels\Percentages = "+format((fuel_list))+"\n")
            else:
                MF_text = "Multifuel = "+format(multifuel)+"\n"   
            text_description = ("This file provides simulation information.\n"
                                "The following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
                                +MF_text+
                                "==============================\n"
                                "\n================Parameters================"
                                "\nInitial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure Range: "
                                +format(P)+" [atm]\nEquivalence Ratio Range: "
                                +format(Phi)+F_O_text+
                                "==========================================\n"
                                "\n======Flame Simulation Information======"
                                "\nMingrid = "+format(mingrid)+
                                "\nMult_Soret = "+format(mul_soret)+
                                "\n========================================\n"
                                "\n=============Time/Converged=============\n"
                                "Total Cases: "+format(len(paramlist))+"\n"
                                "Sim time: "+format(sim_time, '0.5f')+" [s]\n"
                                "Cases Converged: "+format(len(converged))+"\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "========================================")
            if custom:
                text_description = ("The following simulation was a custom"
                                " run using the following properties and "
                                " parameters.\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
                                "==============================\n"
                                "\n================Parameters================"
                                "\nInitial Temperature: "+format(Tint)+
                                " [Kelvin]\nPressure Range: "+format(P)+
                                " [atm]\nFuel Range: "+format(Fuel)+
                                " [mole fraction]\nOxygen Range: "
                                +format(Oxygen)+" [mole fraction]\n"
                                "==========================================\n"
                                "\n======Flame Simulation Information======"
                                "\nMingrid = "+format(mingrid)+
                                "\nMult_Soret = "+format(mul_soret)+
                                "\n========================================\n"
                                "\n=============Time/Converged=============\n"
                                "Total Cases: "+format(len(paramlist))+"\n"
                                "Sim time: "+format(sim_time, '0.5f')+" [s]\n"
                                "Cases Converged: "+format(len(converged))+"\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "========================================")
            f.write(text_description)
            f.close()
            print('\nText file created')
            
            print('\nStart file saving...')
            save_time_start = time.time()
    
            file = os.path.join(save_path, 'Flame Information.pkl')
            with open(file, 'wb') as f:
                pickle.dump(flame_info, f)
            save_time_end = time.time()
            print('End file saving')