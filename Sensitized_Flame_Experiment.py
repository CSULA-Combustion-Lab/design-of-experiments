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
mechanism = 'gri30.cti' #Mechanism file
gas = ct.Solution(mechanism)

#Working directory
flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

#Parameters for main loop
P    = np.logspace(np.log10(1), np.log10(2), 1) #Pressure [atm]
Phi  = np.logspace(np.log10(0.25), np.log10(4), 5) #Equivalence ratio
FtD  = np.logspace(0.06, 0.85, 1) #Fuel to Diluent ratio
OtO  = np.logspace(np.log10(0.1), np.log10(0.95), 5) #Oxygen to Diluent ratio

#Initial Temperature
Tint = 323 #Temperature [K]

#Parameters for mixture (Fuel, Oxidizer, Diluent)
fuel     = 'CH4' #chemical formula of fuel
oxidizer = 'O2' #chemical formula of oxidizer
diluent  = 'N2' #chemical formula of diluent
air      = 'O2:1, N2:3.76' #chemical components for air as an oxidizer

multifuel = False # if true, additional fuels will be added from fuels list
if multifuel:
    fuel = ['CH4', .50 , 'CH3', .50]
    check     = 0
    for c in range(1, len(fuel), 2):
        check += fuel[c]
    if check != 1:
        print('Error in Multifuel.'+
              'Sum of individual fuel percentage must add up to 1!')
        sys.exit()

multioxidizer = True # if true, additional fuels will be added from fuels list
if multioxidizer:
    oxidizer = ['O2', .35 , 'NO2', .65]
    check     = 0
    for c in range(1, len(oxidizer), 2):
        check += oxidizer[c]
    if check != 1:
        print('Error in Multioxidizer.'+
              'Sum of individual fuel percentage must add up to 1!')
        sys.exit()

custom     = False # If true, custom styles used for range and save files
oxi_fuel   = True # If true, OtO will be used instead of Fuel Mole Fraction
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
    P   = 1 #Pressure [atm]
    FtD = 0.05 #Fuel to Diluent Mole Fraction
    OtO = 0.25 #Oxygen to Diluent Mole Fraction
    conditions = {'Parameters': [P, FtD, OtO, Tint],
                  'Mixture': [fuel, diluent, oxidizer, air],
                  'Flame': [mingrid, mul_soret],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxi_fuel, custom, multifuel, multioxidizer]}

#Multifuel loop
else:
    conditions = {'Parameters': [P, Phi, FtD, Tint, OtO],
                  'Mixture': [fuel, diluent, oxidizer, air],
                  'Flame': [mingrid, mul_soret],
                  'Files': [mechanism, flame_temp],
                  'Debug': [Debug_params, LogLevel],
                  'T/F': [debug, oxi_fuel, custom, multifuel, multioxidizer]}


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
        results.append(pool.apply_async(fun, args=(*x, cond)))
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
            datadict[casenum] = {'Flame': [None, 'RunTimeError']}
        casenum += 1
    outlist = [datadict[k] for k in datadict.keys()] # Convert back to list
    return outlist


def flame_sens(p, phi, f_o, cond):
    """[Fill in information]"""
    Tin           = cond ['Parameters'][3]
    chem          = cond['Files'][0]
    tempfile      = cond['Files'][1]
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    debug         = cond['T/F'][0]
    o_f           = cond['T/F'][1]
    custom        = cond['T/F'][2]
    multif        = cond['T/F'][3]
    multio        = cond['T/F'][4]
    mg            = cond['Flame'][0]
    ms            = cond['Flame'][1]
    
# =============================================================================
# Test Section
    Diluent = 1 - f_o
    if custom:
       Fuel = f_o[0]
       Oxidizer = f_o[1]
       
    elif o_f:
       if multif:
           Fuel = ''
           for fl in range(0,len(Fuel_name),2):
                Fuel += Fuel_name[fl]+':'+str(Fuel_name[fl+1])+' '
       else:
           Fuel = Fuel_name
       if multio:
           Oxidizer = ''
           for ol in range(0,len(Oxidizer_name),2):
               Oxidizer += Oxidizer_name[ol]+':'+str(Oxidizer_name[ol+1]*f_o)+' '
           Oxidizer += Diluent_name+':'+str(Diluent)
       else: 
           Oxidizer = Oxidizer_name+':'+str(f_o)+' '+Diluent_name+':'+str(Diluent)

       
    elif not o_f:
       if multif:
           Fuel = ''
           for fl in range(0,len(Fuel_name),2):
                Fuel += Fuel_name[fl]+':'+str(Fuel_name[fl+1]*f_o)+' '
           Fuel += Diluent_name+':'+str(Diluent) 
       else:
           Fuel = Fuel_name+':'+str(f_o)+' '+Diluent_name+':'+str(Diluent)
       if multio:
           Oxidizer = ''
           for ol in range(0,len(Oxidizer_name),2):
               Oxidizer += Oxidizer_name[ol]+':'+str(Oxidizer_name[ol+1])+' '
       else:
           Oxidizer = Oxidizer_name
           
    else:
        print('Error. Check T/F Statments for Incorrect Values.')
        sys.exit()
        
    Mix = mixture_maker(phi, Fuel, Oxidizer)
    Fue_Percent = mixture_percentage(Fuel_name, Mix, multif)
    Oxi_Percent = mixture_percentage(Oxidizer_name, Mix, multio)
# =============================================================================

    # if custom:
    #     Fuel    = Phi
    #     Oxygen  = F_O
    #     Diluent = 1 - Oxygen - Fuel
    #     OtO     = Oxygen/(Oxygen + Diluent)
    #     Phi     = a/Oxygen*Fuel
    #     Mix     = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    # elif oxidizer:
    #     OtO     = F_O
    #     noxy    = (Phi*OtO)/a
    #     Oxygen  = OtO/(noxy + 1)
    #     Fuel    = noxy/(noxy + 1)
    #     Diluent = (1 - OtO)/(noxy + 1)
    #     if multifuel:
    #         mf  = Cond['MultiFuels']
    #         fuel_percent = 0
    #         for n in range(1,len(mf),2):
    #             fuel_percent += mf[n]
    #         Fuel_main = (1 - fuel_percent)*Fuel
    #         Mix       = [[Diluent_name, Diluent], ['O2', Oxygen],
    #                      [Fuel_name, Fuel_main]]
    #         Mix_check = (Diluent_name+':'+format(Diluent)+', O2:'+
    #                      format(Oxygen)+', '+Fuel_name+':'+format(Fuel_main))
    #         for n in range(0,len(mf),2):
    #             Mix.append([mf[n], mf[n+1]*Fuel])
    #             Mix_check += ', '+mf[n]+':'+format(mf[n+1]*Fuel)
    #         gas.X = Mix_check
    #         Phi_check = gas.get_equivalence_ratio()
    #     else:
    #         Mix = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    # elif not oxidizer:
    #     Fuel    = F_O
    #     Oxygen  = a/Phi*Fuel
    #     Diluent = 1 - Oxygen - Fuel
    #     OtO     = Oxygen/(Oxygen + Diluent)
    #     if multifuel:
    #         mf  = Cond['MultiFuels']
    #         fuel_percent = 0
    #         for n in range(1,len(mf),2):
    #             fuel_percent += mf[n]
    #         Fuel_main = (1 - fuel_percent)*Fuel
    #         Mix       = [[Diluent_name, Diluent], ['O2', Oxygen],
    #                      [Fuel_name, Fuel_main]]
    #         Mix_check = (Diluent_name+':'+format(Diluent)+', O2:'+
    #                      format(Oxygen)+', '+Fuel_name+':'+format(Fuel_main))
    #         for n in range(0,len(mf),2):
    #             Mix.append([mf[n], mf[n+1]])
    #             Mix_check += ', '+mf[n]+':'+format(mf[n+1]*Fuel)
    #         gas.X = Mix_check
    #         Phi_check = gas.get_equivalence_ratio()
    #     else:
    #         Mix = [[Diluent_name, Diluent], ['O2', Oxygen], [Fuel_name, Fuel]]
    # else:
    #     print('Error! Check T/F Statments for Incorrect Values.')
    #     sys.exit()
    

    if Diluent < 0:
        Dil_Percent = 0
        flame_info = {'Flame': [None, 'Diluent < 0'],
                      'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent]}
        return flame_info
    else:
        Dil_Percent = 0
        for n in Mix:
            if n[0] == Diluent_name:
                Dil_Percent = n[1]
            else:
                continue
   
    if debug:
        logl = cond['Debug'][1]
        print('\nMixture Composition:')
        print(Mix)     
    else:
        logl = 0
        
    f = flame.Flame(Mix, p, Tin, tempfile, chemfile=chem)
    f.run(mingrid=mg, loglevel=logl, mult_soret=ms)
    if f.flame_result is None:
        flame_info = {'Flame': [None, 'Flame did not converge'],
                      'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent]}
        return flame_info
    f.sensitivity()
    f_sens = f.sens #Rxn sensitivities 
    Su         = f.flame_result.u[0]  #Flame speed at the front
    flame_T    = f.flame_result.T[-1] #Flame temperature at end
    flame_rho  = f.flame_result.density_mass[0] #Flame density at the front
    flame_info = {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
                  'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                 Fuel_name, Oxidizer_name, Diluent_name,
                                 Fue_Percent, Oxi_Percent, Dil_Percent]}
    return flame_info


def duplicate_reactions(gas):
    """[Fill in information]"""
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
    """[Fill in information]"""
    for f in flame_information:
        if f['Flame'][0] is None:
            continue
        else:
            for d in duplicate_rxns:
                sum = 0
                for n in duplicate_rxns[d]:
                     sum += f['Flame'][0][n][1]
                     f['Flame'][0][n][1] = 0
                f['Flame'][0][duplicate_rxns[d][0]][1] = sum


def mixture_maker(phi, fuel, oxidizer):
    """[Fill in information]"""
    gas.set_equivalence_ratio(phi, (fuel), (oxidizer))
    Mix_dict = gas.mole_fraction_dict()
    Mixture = []
    for key, value in Mix_dict.items():
        temp = [key,value]
        Mixture.append(temp)
    return Mixture


def mixture_percentage(fo_list, mix, tf):
    Percentage = 0
    if not tf:
        for m in mix:
            if m[0] == fo_list:
                Percentage += m[1]
            else:
                continue
    else:          
        for n in range(0, len(fo_list), 2):
            for m in mix:
                if m[0] == fo_list[n]:
                    Percentage += m[1]
                else:
                    continue
    return Percentage


if __name__ == "__main__":
#############################Initializing######################################
    """[Fill in information]"""
    tic = time.time() #Main Code Time Start
    if custom: 
        totaliterations = len(P)*len(OtO)
        paramlist       = []
        for i in P:
            for k in OtO:
                paramlist.append((i, FtD, k))
    elif oxidizer:
        totaliterations = len(P)*len(Phi)*len(OtO)
        paramlist       = list(it.product(P,Phi,OtO))
    elif not oxidizer:
        totaliterations = len(P)*len(Phi)*len(FtD)
        paramlist       = list(it.product(P,Phi,FtD))
    else:
        print('Error in Initializing. Check T/F Statments.')
        sys.exit()

##############################Debug loop#######################################
    if debug:
        """[Fill in information]"""
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
        """[Fill in information]"""
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
        converged    = 0
        for x in flame_info:
            if x['Flame'][0] is None:
                continue
            else:
                converged += 1
        flame_info_unfiltered = copy.deepcopy(flame_info)
        duplicate_rxns = duplicate_reactions(gas)
        flame_info_filter(flame_info, duplicate_rxns)
        filter_end  = time.time()
        filter_time = filter_end - filter_start
        print('End of filtering')
        print('Filtering took '+format(filter_time, '0.5f')+ ' seconds.')
        print('Number of cases converged:' +str(converged))
        toc      = time.time() #Main Code Time End
        duration = toc-tic
        print('Total time '+format(duration, '0.5f')+' seconds.\n')
        
#################################Save Files####################################   
    """[Fill in information]"""
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
            """[Fill in information]"""
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
            if oxi_fuel:
                F_O_text  = "\nOxygen to Oxidizer Fraction: "
            else:
                F_O_text  = "\nFuel Mole Fraction: "
            if multifuel:
                MF_text = ("Multifuel = "+format(multifuel)+"\n"
                           "Fuels\Percentages = "+format((fuel))+"\n")
            else:
                MF_text = "Multifuel = "+format(multifuel)+"\n"    
            if multioxidizer:
                MO_text = ("Multioxidizer = "+format(multioxidizer)+"\n"
                           "Oxidizer\Percentages = "+format((oxidizer))+"\n")
            else:
                MO_text = "Multioxidizer = "+format(multioxidizer)+"\n"
            text_description = ("This file provides debug information.\n The "
                                "following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+str(fuel)+"\n"
                                "Oxidizer: "+str(oxidizer)+"\n"
                                "Diluent: "+str(diluent)+"\n"
                                +MF_text+MO_text+
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
            """[Fill in information]"""               
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
            if oxi_fuel:
                F_O_text = ("\nOxygen to Oxidizer Mole Fraction Range: "
                            +format(OtO)+"\n")
            else:
                F_O_text = ("\nFuel Mole Fraction Range: "
                            +format(FtD)+"\n")
            if multifuel:
                MF_text = ("Multifuel = "+format(multifuel)+"\n"
                           "Fuels\Percentages = "+format((fuel))+"\n")
            else:
                MF_text = "Multifuel = "+format(multifuel)+"\n"
            if multioxidizer:
                MO_text = ("Multioxidizer = "+format(multioxidizer)+"\n"
                           "Oxidizer\Percentages = "+format((oxidizer))+"\n")
            else:
                MO_text = "Multioxidizer = "+format(multioxidizer)+"\n"
            text_description = ("This file provides simulation information.\n"
                                "The following information are the parameters "
                                "and cases simulated\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+str(fuel)+"\n"
                                "Oxidizer: "+str(oxidizer)+"\n"
                                "Diluent: "+str(diluent)+"\n"
                                +MF_text+MO_text+
                                "==============================\n"
                                "\n================Parameters================"
                                "\nInitial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure Range: "
                                +format(P)+" [atm]\nEquivalence "
                                +"Ratio Range: "+format(Phi)+F_O_text+
                                "==========================================\n"
                                "\n======Flame Simulation Information======"
                                "\nMingrid = "+format(mingrid)+
                                "\nMult_Soret = "+format(mul_soret)+
                                "\n========================================\n"
                                "\n=============Time/Converged=============\n"
                                "Total Cases: "+format(len(paramlist))+"\n"
                                "Sim time: "+format(sim_time, '0.5f')+" [s]\n"
                                "Cases Converged: "+str(converged)+"\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "========================================")
            if custom:
                """[Fill in information]"""
                text_description = ("The following simulation was a custom"
                                " run using the following properties and "
                                " parameters.\n\n"
                                "==========Properties==========\n"
                                "Mechanism: "+mechanism+"\n"
                                "Fuel: "+str(fuel)+"\n"
                                "Oxidizer: "+str(oxidizer)+"\n"
                                "Diluent: "+str(diluent)+"\n"
                                "==============================\n"
                                "\n================Parameters================"
                                "\nInitial Temperature: "+format(Tint)+
                                " [Kelvin]\nPressure Range: "+format(P)+
                                " [atm]\nFuel Range: "+format(FtD)+
                                " [mole fraction]\nOxygen Range: "
                                +format(OtO)+" [mole fraction]\n"
                                "==========================================\n"
                                "\n======Flame Simulation Information======"
                                "\nMingrid = "+format(mingrid)+
                                "\nMult_Soret = "+format(mul_soret)+
                                "\n========================================\n"
                                "\n=============Time/Converged=============\n"
                                "Total Cases: "+format(len(paramlist))+"\n"
                                "Sim time: "+format(sim_time, '0.5f')+" [s]\n"
                                "Cases Converged: "+str(converged)+"\n"
                                "Run time: "+format(duration, '0.5f')+" [s]\n"
                                "========================================")
            f.write(text_description)
            f.close()
            print('Text file created')
            
            print('\nStart file saving...')
            save_time_start = time.time()
    
            file = os.path.join(save_path, 'Flame Information.pkl')
            with open(file, 'wb') as f:
                pickle.dump(flame_info, f)
            save_time_end = time.time()
            print('End file saving')