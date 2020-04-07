# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 09:52:04 2020

@author: Kodo Bear
"""

import os
import errno
import time, sys
import numpy as np
# import cantera as ct
import itertools as it
path = (r'D:\School\Cal State La\Extracurricular\NOx Combustion Research')
sys.path.insert(0, path)
from utilities import flame
from datetime import datetime
# import math
import pickle

####Set experiment parameters
mechanism = 'h2_burke2012.cti' #Mechanism file

#Working directory
flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

#Parameters for main loop
P    = np.logspace(np.log10(1), np.log10(100), 2) #Pressure [atm]
Phi  = np.logspace(np.log10(0.1), np.log10(1), 2) #Equivalence ratio
Fuel = np.logspace(np.log10(0.1), np.log10(0.85), 2) #Fuel mole fraction

#Diluent/Oxygen Ratio
DtO = np.linspace(1, 10, 2) #Maybe use this instead of fuel parameter 

#Initial Temperature
Tint = 273.15 #Temperature [K]

#Parameters for mixture
fuel_name = 'H2' #chemical formula of fuel

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

save_files = True  # If true, save files for plotting script
debug      = False  # If true, print lots of information for debugging.

DEBUG_FMT = 'Removing condition: T={:.0f}, P={:.0f}, phi={:.3g}, fuel={:.3g}'
#Debug parameters [Pressure, Equivalence Ratio, Fuel, Temperature]
Debug_params = [1, 1, 0.42, 300]
LogLevel     = 1

conditions = {'Parameters': [P, Phi, Fuel, Tint],
              'Mixture': [fuel_name, x, y, z, a, diluent_name],
              'Files': [mechanism, flame_temp],
              'Debug': [Debug_params, LogLevel]}

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

def flame_sens(P, Phi, Fuel, Tin, Cond, debug=False):
    """[Fill in information]"""
    chem         = Cond['Files'][0]
    tempfile     = Cond['Files'][1]
    a            = Cond['Mixture'][4]
    Fuel_name    = Cond['Mixture'][0]
    Diluent_name = Cond['Mixture'][5]
    Oxygen       = a/Phi*Fuel
    Diluent      = 1 - Oxygen - Fuel
    Mix = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
   
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
                      'Conditions': [P, Fuel, Phi, Tin, Mix]}
        return flame_info
    f.sensitivity()
    flame_sens = f.sens
    Su         = f.flame_result.u[0]
    flame_info = {'Flame': [flame_sens, Su],
                  'Conditions': [P, Fuel, Phi, Tin, Mix]}
    return flame_info



if __name__ == "__main__":
    #Start time
    tic = time.time()
    #Initializing
    iteration       = 0
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
        flame_info_debug = flame_sens(*Debug_params,
                                      conditions, debug)
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
            text_description = ("This file provides debug information.\n The "
                                "following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Parameters==========\n"
                                "Initial Temperature: "+format(Debug_params[3])
                                +" [Kelvin]\nPressure: "
                                +format(Debug_params[0])+" [atm]\n"
                                "Equivalence Ratio: "+format(Debug_params[1])+
                                "\nFuel Mole Fraction: "
                                +format(Debug_params[2])+"\n"
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
            text_description = ("This file provides simulation information.\n"
                                "The following information are the parameters "
                                "and cases simulated\n\n"
                                "================Parameters================"
                                "\nInitial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure Range: "
                                +format(P)+" [atm]\nEquivalence Ratio Range: "
                                +format(Phi)+"\nFuel Mole Fraction Range: "
                                +format(Fuel)+"\n"
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