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
Simulation_Type = '1D'

# TODO: Check mix type names. Keep naming argument similar. Who reads this next must do it!
#Mixture Type (Debug, Custom, Oxi_Dil, Fue_Dil)
#Provide one of the four types of mixtures into
# variable mixture_type as a string
#  Debug is used for simulating a single flame to check code
#  Custom is under construction
#  Oxi_Dil creates a mixture where the Diluent is a ratio of the Oxidizer used
#  Fue_Dil creates a mixture where the Diluent is a ratio of the Fuel used
#  phi_fuel specifies the equivalence ratio and fuel mole fraction
#  Ox_fuel specifies the oxidizer and fuel mole fractions
Mixture_type = 'Ox_Fuel'

#Parameters (Pressure, Equivalence, Temperature, Diluent Percentage)
## All parameters are lists of three numbers
### First Number: Initial point
### Second Number: End Point
### Third Number: Number of points. If set to 1 only first point is evaluated
## Units
### Pressure: Atmosphere [atm]
### Equivalence: Dimensionless (<1:Fuel Lean, 1:Unity, >1:Fuel Rich)
### Temperature: Kelvin [K]
### Diluent_Percentage: Percent [%]
Pressure           = [0.5, 1, 2]
Temperature        = [373, 400, 1]

# Mixture parameters. Not all of these will be used, depending on Mixture_type
Equivalence        = [0.8, 1, 2]
fraction_in_oxidizer_or_fuel = [0.5, 0.7, 2]  # X / (X+diluent) where X is O2 or fuel. Used for Oxi_dil or Fue_dil
fuel_fraction      = [0.1, 0.5, 3]
oxidizer_fraction  = [0.2, 0.6, 3]

### Array_type: 'log' or 'lin', specifying if thermodynamic and mixture
###             variables should vary in log- or linear- space
Array_type = 'log'

#Set experiment parameters
Mechanism = 'Li_model_modified_trioxane.cti' #Mechanism file

#Parameters for mixture (Fuel, Oxidizer, Diluent)
# Fuel and Oxidizer can either be a single chemical string or multichemical list
#  For multichemical follow structure ['Chemical1', % of Total, ...]
#  The percentage following the chemical name should sum up to 1
# Diluent should be a single chemical added to the mixture
Fuel     = 'C3H6O3' #chemical formula of fuel
# Fuel  = ['CH4', .50 , 'CH3OH', .50]
Oxidizer = 'O2' #chemical formula of oxidizer
# Oxidizer = ['O2', .35 , 'NO2', .65]
Diluent  = 'N2' #chemical formula of diluent

#################### Missing code here #########################
#   Add code to define F_to_D and O_to_D. Or, are these necessary?
#################################################################

#Flame Conditions
Mingrid   = 200
Mul_soret = False
Loglevel  = 0

#True/False statements
Save_files = True # If true, save files for plotting script
Save_time  = True # Only used for zeroD. If true, saves time for GUI
if __name__ == "__main__":

    # Package the mixture parameters. Each mixture type requires two variables
    if Mixture_type in ('Oxi_Dil', 'Fue_Dil'):
        mix_params = (Mixture_type, Equivalence, fraction_in_oxidizer_or_fuel)
    elif Mixture_type == 'phi_fuel':
        mix_params = (Mixture_type, Equivalence, fuel_fraction)
    elif Mixture_type == 'Ox_Fuel':
        mix_params = (Mixture_type, oxidizer_fraction, fuel_fraction)
    else:
        raise ValueError('Mixture_type = {} is not supported'.format(Mixture_type))
        
    # TODO: 0D needs to take in mix_params tuple as seen above.
    # Example of implementation can be seen in 1D code below.
    if Simulation_Type == '0D':
        print(cf.parameters_string(Pressure, Temperature, mix_params))
        zeroD.run_0d_simulation(Mechanism, Array_type, Pressure, Temperature,
                                Fuel, Oxidizer, Diluent,
                                mix_params, Save_files, Save_time)
        print('Under-Construction!')
    elif Simulation_Type =='1D':
        print(cf.parameters_string(Pressure, Temperature, mix_params))
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