# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:32:37 2020

@author: Kodo Bear
"""
# import code_2_5 as zeroD
import Sensitized_Flame_Experiment as oneD
# import BurnerSimulation
import sys

#Initializer
# Simulation Types (0D, 1D, Burner)
## Variable Simulation_Type input is a string
## 0D Simulates a no flame condition tracking species concentraion in time
## 1D Simulates an adiabatic flame tracking species before, at, and after flame
## Burner Simulates a planar flame stabalized at some distance from burner
Simulation_Type = '1D'

#Mixture Type (Debug, Custom, Oxi_Dil, Fue_Dil)
#Provide one of the four types of mixtures into
# variable mixture_type as a string
#  Debug is used for simulating a single flame to check code
#  Custom is under construction
#  Oxi_Dil creates a mixture where the Diluent is a ratio of the Oxidizer used
#  Fue_Dil creates a mixture where the Diluent is a ratio of the Fuel used
Mixture_type = 'Oxi_Dil'

#Parameters (Pressure, Equivalence, Temperature, Diluent Percentage)
## All three parameters are lists of three numbers
### First Number: Initial point
### Second Number: End Point
### Third Number: Number of point. If set to 1 only first point is evaluated
## Units
### Pressure: Atmosphere [atm]
### Equivalence: Dimensionless (<1:Fuel Lean, 1:Unity, >1:Fuel Rich)
### Temperature: Kelvin [K]
### Diluent_Percentage: Percent [%]
### Array_type: 'log' or 'lin', specifying if thermodynamic and mixture
###             variables should vary in log- or linear- space
Pressure           = [0.5, 1, 2]
Equivalence        = [0.05, 1, 2]
Temperature        = [373, 400, 1]
Diluent_Percentage = [0.05, 0.95, 2]
Array_type = 'log'

#Set experiment parameters
Mechanism = 'Li_model_modified_trioxane.cti' #Mechanism file

#Initial temperature of unburned mixture
Tint = 373 #Temperature [K]

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
Air      = 'O2:1, N2:3.76' #chemical components for air as an oxidizer

#################### Missing code here #########################
#   Add code to define F_to_D and O_to_D. Or, are these necessary?
#################################################################

#Flame Conditions
Mingrid   = 200
Mul_soret = False
Loglevel  = 0

#True/False statements
Save_files = True # If true, save files for plotting script

if __name__ == "__main__":

    if Simulation_Type == '0D':
        # zeroD.run_flow_reactor_simulation(Parameter_List, Mechanism_List)
        print('Under-Construction!')
    elif Simulation_Type =='1D':
        oneD.run_flame_simulation(Mechanism, Array_type, Pressure, Equivalence,
                                Diluent_Percentage, Tint, Fuel, Oxidizer, Diluent,
                                Air, Mingrid, Mul_soret, Loglevel,
                                Mixture_type, Save_files)
    elif Simulation_Type == 'Burner':
        # Burner = BurnerSimulation(Parameter_List, Mechanism_List)
        print('Under Construction!')
    else:
        print('Error! Simulation_Type String Does Not Match!'+
              '\nMake sure string matches one of three options!')
        sys.exit()