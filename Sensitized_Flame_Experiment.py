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

#Set experiment parameters
mechanism = 'Li_model_modified_trioxane.cti' #Mechanism file
gas = ct.Solution(mechanism)

#Working directory
flame_temp = os.path.join(r'Flame_Files', 'temp_flame_files')

#Parameters for main loop
# P and Phi are used in all cases.
# Whether FTD or OtD is used depends on mixture type.
# Each list should follow patter [First Point, Last Point, Number of Points]
# Array_type's (log, lin)
#  log creates a logspace array of parameters
#  lin creates a linspace array of parameters
array_type = 'lin'
Press      = [0.25, 2, 8]  #Pressure [atm]
E_Ratio    = [0.1, 1, 12] #Equivalence ratio
F_to_D     = [0.75, 0.95, 2] #Fuel/(Fuel + Diluent)
O_to_D     = [0.05, 0.95, 12] #Oxidizer/(Oxidizer + Diluent)
if array_type == 'log':
    P    = np.logspace(np.log10(Press[0]), np.log10(Press[1]), Press[2])
    Phi  = np.logspace(np.log10(E_Ratio[0]), np.log10(E_Ratio[1]), E_Ratio[2])
    FtD  = np.logspace(np.log10(F_to_D[0]), np.log10(F_to_D[1]), F_to_D[2])
    OtD  = np.logspace(np.log10(O_to_D[0]), np.log10(O_to_D[1]), O_to_D[2])
elif array_type == 'lin':
    P   = np.linspace(Press[0], Press[1], Press[2])
    Phi = np.linspace(E_Ratio[0], E_Ratio[1], E_Ratio[2])
    FtD = np.linspace(F_to_D[0], F_to_D[1], F_to_D[2])
    OtD = np.linspace(O_to_D[0], O_to_D[1], O_to_D[2])
else:
    print('Error! Check array_type variable for invalid string input')
    sys.exit()

#Initial temperature of unburned mixture
Tint = 373 #Temperature [K]

#Parameters for mixture (Fuel, Oxidizer, Diluent)
fuel     = 'C3H6O3' #chemical formula of fuel
# fuel  = ['CH4', .65 , 'CO2', .35]
oxidizer = 'O2' #chemical formula of oxidizer
# oxidizer = ['O2', .35 , 'NO2', .65]
diluent  = 'N2' #chemical formula of diluent
air      = 'O2:1, N2:3.76' #chemical components for air as an oxidizer

#Flame Conditions
mingrid   = 200
mul_soret = False
loglevel  = 0

#True/False statements
save_files    = True # If true, save files for plotting script

if type(oxidizer) is list:
    multioxidizer = True
elif type(oxidizer) is str:
    multioxidizer = False

if type(fuel) is list:
    multifuel = True
elif type(fuel) is str:
    multifuel = False

#Mixture Type (Debug, Custom, Oxi_Dil, Fue_Dil)
#Provide one of the four types of mixtures into
# variable mixture_type as a string
#  Debug is used for simulating a single flame to check code
#  Custom is under construction
#  Oxi_Dil creates a mixture where the Diluent is a ratio of the Oxidizer used
#  Fue_Dil creates a mixture where the Diluent is a ratio of the Fuel used
mixture_type = 'Oxi_Dil'

#Multifuel mixture percentage of total fuel
# fuel is a list of fuels and their percentages in the total fuel
#  odds are fuel name as a string
#  evens are percetage of previous fuel name in total fuel
# percentages should sum to 1 or script will not run
if multifuel:
    check = 0
    for c in range(1, len(fuel), 2):
        check += fuel[c]
    if not np.isclose(check, 1):
        print('Error in Multifuel.'+
              'Sum of individual fuel percentage must add up to 1!')
        sys.exit()

#Multioxidizer mixture percentage of total oxidizer
# fuel is a list of fuels and their percentages in the total fuel
#  odds are fuel name as a string
#  evens are percetage of previous fuel name in total fuel
# percentages should sum to 1 or script will not run
if multioxidizer:
    check    = 0
    for c in range(1, len(oxidizer), 2):
        check += oxidizer[c]
    if not np.isclose(check, 1):
        print('Error in Multioxidizer.'+
              'Sum of individual fuel percentage must add up to 1!')
        sys.exit()

#Debug loop parameters
if mixture_type == 'Debug':
    DEBUG_FMT = ('Removing condition: T={:.0f}, P={:.0f}, Phi={:.3g}, '+
                 ' FtD={:.3g}, OtD={:.3g}')
    #Debug Mixture [Pressure, Equivalence Ratio, Fuel String,
    #               Oxygen String, Initial Temperature of Unburned Mixture]
    P            = 1 #Pressure [atm]
    Phi          = 1 #Equivalence Ratio
    fuel         = 'CH4:1' #Fuel String
    oxidizer     = air #Oxidizer String
    Tint         = 300 #Initial Temperatur of Unburned Mixture
    loglevel     = 1
    Debug_params = [P, Phi, [fuel, oxidizer]]
    conditions = {'Parameters': [P, Phi, FtD, Tint, OtD, array_type],
                  'Mixture': [fuel, diluent, oxidizer, mixture_type],
                  'Flame': [mingrid, mul_soret, loglevel],
                  'Files': [mechanism, flame_temp],
                  'T/F': [multifuel, multioxidizer]}

#Custom loop
elif mixture_type == 'Custom':
    P   = [0.5, 1] #Pressure [atm]
    Phi = [0.25] #Equivalence Ratio
    FtD = [0.05] #Fuel to Diluent Mole Fraction
    OtD = [0.25] #Oxygen to Diluent Mole Fraction
    conditions = {'Parameters': [P, Phi, FtD, Tint, OtD, array_type],
                  'Mixture': [fuel, diluent, oxidizer, mixture_type],
                  'Flame': [mingrid, mul_soret, loglevel],
                  'Files': [mechanism, flame_temp],
                  'T/F': [multifuel, multioxidizer]}

#Multifuel loop
elif mixture_type == 'Oxi_Dil' or mixture_type == 'Fue_Dil':
    conditions = {'Parameters': [P, Phi, FtD, Tint, OtD, array_type],
                  'Mixture': [fuel, diluent, oxidizer, mixture_type],
                  'Flame': [mingrid, mul_soret, loglevel],
                  'Files': [mechanism, flame_temp],
                  'T/F': [multifuel, multioxidizer]}


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
    at            = cond['Parameters'][5]
    chem          = cond['Files'][0]
    tempfile      = cond['Files'][1]
    Fuel_name     = cond['Mixture'][0]
    Diluent_name  = cond['Mixture'][1]
    Oxidizer_name = cond['Mixture'][2]
    mt            = cond['Mixture'][3]
    multif        = cond['T/F'][0]
    multio        = cond['T/F'][1]
    mg            = cond['Flame'][0]
    ms            = cond['Flame'][1]
    logl          = cond['Flame'][2]

    if not mt == 'Debug':
        Diluent = 1 - f_o #Define Diluent Percentage in Fuel/Oxidizer

    #Four Mixture Types (Debug, Custom, Oxidizer to Diluent, Fuel to Diluent)
    #Debug Mixture is a single mixture to test for errors in code
    #Custom Mixture is under construction
    #o_f True Mixture uses a ratio of Oxidizer to Diluent
    #o_f False Mixture uses a ratio of Fuel to Diluent

    if mt == 'Debug':
        Fuel     = f_o[0]
        Oxidizer = f_o[1]
        Diluent  = 0

    elif mt == 'Custom':
       #Under Construction#
       Fuel = f_o[0]
       Oxidizer = f_o[1]
       #Under Construction#

    elif mt == 'Oxi_Dil':
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

    elif mt == 'Fue_Dil':
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
        print('Error. Check mixture_type variable for invalid string input.')
        sys.exit()

    Mix = mixture_maker(phi, Fuel, Oxidizer)
    Fue_Percent = mixture_percentage(Fuel_name, Mix, multif)
    Oxi_Percent = mixture_percentage(Oxidizer_name, Mix, multio)

    if mt == 'Debug':
        print('\nMixture Composition:')
        print(Mix)

    if Diluent < 0:
        Dil_Percent = 0
        flame_info = {'Flame': [None, 'Diluent < 0'],
                      'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent,
                                     at]}
        return flame_info
    else:
        Dil_Percent = 0
        for n in Mix:
            if n[0] == Diluent_name:
                Dil_Percent = n[1]
            else:
                continue

    f = flame.Flame(Mix, p, Tin, tempfile, chemfile=chem)
    f.run(mingrid=mg, loglevel=logl, mult_soret=ms)
    if f.flame_result is None:
        flame_info = {'Flame': [None, 'Flame did not converge'],
                      'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                     Fuel_name, Oxidizer_name, Diluent_name,
                                     Fue_Percent, Oxi_Percent, Dil_Percent,
                                     at]}
        return flame_info
    f.sensitivity()
    f_sens = f.sens #Rxn sensitivities
    Su         = f.flame_result.u[0]  #Flame speed at the front
    flame_T    = f.flame_result.T[-1] #Flame temperature at end
    flame_rho  = f.flame_result.density_mass[0] #Flame density at the front
    flame_info = {'Flame': [f_sens, Su, flame_rho, flame_T, mg, ms],
                  'Conditions': [Tin, p, phi, Fuel, Oxidizer, Mix,
                                 Fuel_name, Oxidizer_name, Diluent_name,
                                 Fue_Percent, Oxi_Percent, Dil_Percent,
                                 at]}
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
    if mixture_type == 'Debug':
        print('Debug Loop Enabled')
    elif mixture_type == 'Custom':
        print('Custom Loop Enabled')
        #Under Construction
        print('Error! Custom is under construction! Choose another mixture!')
        sys.exit()
        # totaliterations = len(P)*len(OtD)
        # paramlist       = []
        # for i in P:
        #     for k in OtD:
        #         paramlist.append((i, FtD, k))
        #Under Construction
    elif mixture_type == 'Oxi_Dil':
        print('Oxidizer to Diluent Loop Enabled')
        totaliterations = len(P)*len(Phi)*len(OtD)
        paramlist       = list(it.product(P,Phi,OtD))
    elif mixture_type == 'Fue_Dil':
        print('Fuel to Diluent Loop Enabled')
        totaliterations = len(P)*len(Phi)*len(FtD)
        paramlist       = list(it.product(P,Phi,FtD))
    else:
        print('Error in Initializing. Check mixture_type variable.')
        sys.exit()

##############################Debug loop#######################################
    if mixture_type == 'Debug':
        """[Fill in information]"""
        print('Debugging in process..')
        print('\nConditions Used:'
              '\nPressure: '+format(Debug_params[0])+' [atm]'
              '\nPhi: '+format(Debug_params[1])+
              '\nFuel: '+format(Debug_params[2][0])+
              '\nOxidizer: '+format(Debug_params[2][1])+
              '\nUnburned Temperature: '+format(Tint)+' [K]')
        flame_info_debug = flame_sens(*Debug_params, conditions)
        print('\nDebuggin complete!')
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
        print('\nNumber of cases converged:' +str(converged))
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
        if mixture_type == 'Debug':
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
            if mixture_type == 'Oxi_Dil':
                F_O_text  = "\nOxidizer to Diluent Fraction: "
            elif mixture_type == 'Fue_Dil':
                F_O_text  = "\nFuel to Diluent Fraction: "
            else:
                print('Error, check string in text description')
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
                                "Initial Temperature: "+format(Tint)
                                +" [Kelvin]\nPressure: "
                                +format(Debug_params[0])+" [atm]\n"
                                "Equivalence Ratio: "+format(Debug_params[1])+
                                F_O_text+format(Debug_params[2][0])+"\n"
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
            if mixture_type == 'Oxi_Dil':
                F_O_text = ("\nOxidizer to Diluent Mole Fraction Range: "
                            +format(OtD)+"\n")
            elif mixture_type == 'Fue_Dil':
                F_O_text = ("\nFuel Mole Fraction Range: "
                            +format(FtD)+"\n")
            else:
                print('Error, strings in text description')
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
            if mixture_type == 'Custom':
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
                                +format(OtD)+" [mole fraction]\n"
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