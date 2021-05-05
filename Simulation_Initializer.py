# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:32:37 2020

@author: Kodo Bear
"""
import Sensitized_0D_Experiment as zeroD
import Sensitized_Flame_Experiment as oneD
import common_functions as cf
# import BurnerSimulation
import sys

#Initializer
# Simulation Types (0D, 1D, Burner)
## Variable Simulation_Type input is a string
## 0D Simulates a no flame condition tracking species concentraion in time
## 1D Simulates an adiabatic flame tracking species before, at, and after flame
## Burner Simulates a planar flame stabalized at some distance from burner
Simulation_Type = '0D'

#Provide one of types of mixtures into variable mixture_type as a string
#  phi_oxi/dil specifies the equivalence ration and oxidizer to diluent ratio
#  phi_fuel/dil specifies the equivalence ration and fuel to diluent ratio
#  phi_fuel specifies the equivalence ratio and fuel mole fraction
#  phi_oxi specifies the equivalence ratio and oxidizer mole fraction
#  phi_dil specifies the equivalence ratio and diluent mole fraction
#  fuel_dil specifies the fuel and diluent mole fractions
#  oxi_dil specifies the oxidizer and diluent mole fractions
#  oxi_fuel specifies the oxidizer and fuel mole fractions
Mixture_type = 'phi_fuel'

#Parameters (Pressure, Temperature, Variable 1, Variable 2)
## All parameters are lists of three numbers
### First Number: Initial point
### Second Number: End Point
### Third Number: Number of points. If set to 1 only first point is evaluated
## Units
### Pressure: Atmosphere [atm]
### Temperature: Kelvin [K]
### Equivalence: Dimensionless (<1:Fuel Lean, 1:Unity, >1:Fuel Rich)
### fraction_in_oxidizer_or_fuel: Percentage of Total Fuel or Oxidizer [%]
### fuel_fraction: Percentage of Total Mixture [%]
### oxidizer_fraction: Percentage of Total Mixture [%]
Pressure    = [.1, 100, 8]
Temperature = [600, 2500, 8]

# Mixture parameters. Not all of these will be used, depending on Mixture_type
Equivalence                  = [0.00025, 2.5, 10] # Used in mixture types with phi
fraction_in_oxidizer_or_fuel = [0.00001, 0.1, 4]  # X / (X+diluent) where X is O2 or fuel. Used for Oxi_dil or Fue_dil
fuel_fraction                = [0.00001, 0.1, 10] # Used in mixture types with fuel
oxidizer_fraction            = [0.01, 0.9, 15] # Used in mixture types with oxi

### Array_type: 'log' or 'lin', specifying if thermodynamic and mixture
###             variables should vary in log- or linear- space
Array_type = 'lin'

#Set experiment parameters
Mechanism = 'mech-FFCM1_modified.cti' #Mechanism file

#Parameters for mixture (Fuel, Oxidizer, Diluent)
# Fuel and Oxidizer can either be a single chemical string or multichemical list
#  For multichemical follow structure ['Chemical1', % of Total, ...] or use
#  a dictionary.
#  The percentage following the chemical name should sum up to 1
# Diluent should be a single chemical added to the mixture
Fuel     = 'H2' # Chemical formula of fuel
# Fuel     = ['CH4', .75, 'CO2', .25]
Oxidizer = 'O2' # Chemical formula of oxidizer
# Oxidizer = ['O2', .35 , 'NO2', .65]
Diluent  = 'N2' # Chemical formula of diluent

#Zero Dimensional Conditions
SpecificSpecies = ['OH'] #Species of interest for rxn ranking data
Starttime       = 0     # In the case a reading is too early
Endtime         = 0.001 # One milisecond
Delta_T         = 100
threshold       = 1/1000000 # Ignore measurements below this mole fraction

#Flame Conditions
Mingrid   = 200
Mul_soret = False
Loglevel  = 0

#Boolean Statements
Save_files = True # If true, save files for plotting script
Save_time  = True # Only used for zeroD. If true, saves time for GUI

if __name__ == "__main__":

    # Package the mixture parameters. Each mixture type requires two variables
    if Mixture_type in ('phi_oxi/dil', 'phi_fuel/dil'):
        mix_params = (Mixture_type, Equivalence, fraction_in_oxidizer_or_fuel)
    elif Mixture_type == 'phi_fuel':
        mix_params = (Mixture_type, Equivalence, fuel_fraction)
    elif Mixture_type == 'oxi_fuel':
        mix_params = (Mixture_type, oxidizer_fraction, fuel_fraction)
    else:
        raise ValueError('Mixture_type = {} is not supported'.format(Mixture_type))

    if Simulation_Type == '0D':
        print(cf.parameters_string(Pressure, Temperature, mix_params,
                                   Mechanism, Fuel, Oxidizer, Diluent))
        zeroD.run_0D_simulation(Mechanism, Array_type, Pressure, Temperature,
                                Fuel, Oxidizer, Diluent, Mixture_type,
                                mix_params, SpecificSpecies, Starttime,
                                Endtime, Delta_T, threshold, Save_files, Save_time)
    elif Simulation_Type =='1D':
        print(cf.parameters_string(Pressure, Temperature, mix_params,
                                   Mechanism, Fuel, Oxidizer, Diluent))
        oneD.run_flame_simulation(Mechanism, Array_type, Pressure, Temperature,
                                  Fuel, Oxidizer, Diluent, mix_params, Mingrid,
                                  Mul_soret, Loglevel, Save_files)
    elif Simulation_Type == 'Burner':
        # Burner = BurnerSimulation(Parameter_List, Mechanism_List)
        print('Under Construction!')
    else:
        print('Error! Simulation_Type String Does Not Match!'+
              '\nMake sure string matches one of three options!')
        sys.exit()