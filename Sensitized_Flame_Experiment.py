# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:52:04 2020

@author: Kodo Bear
"""

import os
import errno
import time, sys
import numpy as np
import itertools as it
path = (r'D:\School\Cal State La\Extracurricular\NOx Combustion Research')
sys.path.insert(0, path)
from utilities import flame
from datetime import datetime
import pickle

####Set experiment parameters
mechanism = 'Li_model.cti' #Mechanism file

#Working directory
flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

#Parameters for main loop
P    = np.logspace(np.log10(1), np.log10(10), 2) #Pressure [atm]
Phi  = np.linspace(0.1, 1.5, 2) #Equivalence ratio
Fuel = np.linspace(0.1, 0.85, 8) #Fuel mole fraction
OtO  = np.linspace(.077, .21, 2) #Oxygen to Oxidizer ratio [Air = .21]


#Initial Temperature
Tint = 323 #Temperature [K]

#Parameters for mixture
fuel_name = 'CH3OH' #chemical formula of fuel

#Fuel C(x)H(y)O(z)
if fuel_name.find('C',0) == -1:
    x = 0
else:
    x = int(fuel_name[fuel_name.find('C',0)+1]) #moles of carbon in fuel
if fuel_name.find('H',0) == -1:
    y = 0
else:
    y = int(fuel_name[fuel_name.find('H',0)+1]) #moles of hydrogen in fuel
if fuel_name.find('O',0) == -1:
    z = 0
else:
    z = int(fuel_name[fuel_name.find('O',0)+1]) #moles of oxygen in fuel
    
a = x+y/4-z/2       #molar oxygen-fuel ratio
diluent_name = 'N2' #chemical formula of diluent


custom     = False # If true, custom styles used for range and save files
oxidizer   = True # If true, OtO will be used instead of Fuel Mole Fraction
debug      = False # If true, print lots of information for debugging.
save_files = False # If true, save files for plotting script

#Debug Files
DEBUG_FMT = 'Removing condition: T={:.0f}, P={:.0f}, phi={:.3g}, fuel={:.3g}'
#Debug parameters [Pressure, Equivalence Ratio, Fuel or Oxygen, Temperature]
Debug_params = [1, 1, 0.21, 300]
LogLevel     = 1


#Custom loop
if custom:
    P      = [0.5, 1]
    Fuel   = 0.05
    Oxygen = np.linspace(0.05,0.95,50)
    conditions = {'Parameters': [P, Fuel, Oxygen, Tint],
                  'Mixture': [fuel_name, x, y, z, a, diluent_name],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxidizer, custom]}
#Normal loop
else:
    conditions = {'Parameters': [P, Phi, Fuel, Tint, OtO],
                  'Mixture': [fuel_name, x, y, z, a, diluent_name],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxidizer, custom]}


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
        Mix     = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    else:
        Fuel    = F_O
        Oxygen  = a/Phi*Fuel
        Diluent = 1 - Oxygen - Fuel
        OtO     = Oxygen/(Oxygen + Diluent)
        Mix     = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
        
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
    f.run(mingrid=200, loglevel=logl)
    if f.flame_result is None:
        flame_info = {'Flame': None,
                      'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
        return flame_info
    f.sensitivity()
    flame_sens = f.sens #Rxn sensitivities 
    Su         = f.flame_result.u[0]  #Flame speed at the front
    flame_T    = f.flame_result.T[-1] #Flame temperature at end
    flame_rho  = f.flame_result.density_mass[0] #Flame density at the front
    flame_info = {'Flame': [flame_sens, Su, flame_rho, flame_T],
                  'Conditions': [P, Fuel, Phi, Tin, Mix, OtO]}
    return flame_info



if __name__ == "__main__":
    #Start time
    tic = time.time()
    #Initializing
    iteration = 0
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

    #Debug loop
    if debug:
        print('Debugging in process..')
        print('\nConditions Used:'
              '\nPressure: '+format(Debug_params[0])+' [atm]'
              '\nPhi: '+format(Debug_params[1])+
              '\nFuel: '+format(Debug_params[2])+
              '\nTemperature: '+format(Debug_params[3])+' [K]')
        flame_info_debug = flame_sens(*Debug_params, conditions)
        print('Debuggin complete!')
        toc = time.time()
        duration = toc - tic
        print('Dubugging time: '+format(duration, '0.5f')+' seconds\n')
   
    #Simulation loop
    else:
        print('Initial number of cases: '+format(len(paramlist)))
        print('\nStart of simulations...')
        sim_start = time.time()
        flame_info = parallelize(paramlist, conditions, flame_sens)
        sim_end    = time.time()
        sim_time   = sim_end - sim_start
        print('End of simulations')
        print('Simulations took '+format(sim_time, '0.5f')+' seconds.')
        converged = []
        for x in flame_info:
            if x['Flame'] == None:
                continue
            else:
                converged.append(1)
        toc = time.time()
        print('Number of cases converged:', len(converged))
        duration = toc-tic
        print('Total time '+format(duration, '0.5f')+' seconds.\n')
    
    # Save files
    if save_files:
        #Save Path/Parent Directory
        parent_dir = 'Flame_Sensitivity_Results'
        try:
            os.makedirs(parent_dir)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        
        #Debug directory
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
            text_description = ("This file provides debug information.\n The "
                                "following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
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
         
        #Normal directory
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
            text_description = ("This file provides simulation information.\n"
                                "The following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
                                "==============================\n"
                                "================Parameters================"
                                "\nInitial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure Range: "
                                +format(P)+" [atm]\nEquivalence Ratio Range: "
                                +format(Phi)+F_O_text+
                                "==========================================\n"
                                "\n=============Time/Converged=============\n"
                                "Sim time: "+format(sim_time, '0.5f')+" [s]\n"
                                "Cases Converged: "+format(len(converged))+"\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "========================================")
            if custom:
                text_description = ("The following simulation was a custom"
                                " run using the follow properties and "
                                " parameters.\n\n"
                                "==========Properties==========\n"
                                "Mechanism "+mechanism+"\n"
                                "Fuel: "+fuel_name+"\n"
                                "Diluent: "+diluent_name+"\n"
                                "==============================\n"
                                "================Parameters================"
                                "\nInitial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure Range: "+format(P)+
                                " [atm]\nFuel Range: "+format(Fuel)+
                                " [mole fraction]\nOxygen Range: "
                                +format(Oxygen)+" [mole fraction]\n"
                                "==========================================\n"
                                "\n=============Time/Converged=============\n"
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