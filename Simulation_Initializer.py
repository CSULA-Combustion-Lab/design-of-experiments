# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 13:32:37 2020

@author: Kodo Bear
"""
# import code_2_5
import Sensitized_Flame_Experiment
# import BurnerSimulation

#Initializer
# Simulation Types (0D, 1D, Burner)
## Variable Simulation_Type input is a string
## 0D Simulates a no flame condition tracking species concentraion in time
## 1D Simulates an adiabatic flame tracking species before, at, and after flame
## Burner Simulates a planar flame stabalized at some distance from burner
Simulation_Type = '1D'

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
Pressure           = [0.5, 1, 2]
Equivalence        = [0.05, 1, 12]
Temperature        = [300, 400, 2]
Diluent_Percentage = [0.05, 0.95, 12]

#Mechanism File
Mechanism = 'Li_model.cti'
Fuel      = 'CH3OH'
Oxidizer  = 'O2'
Diluent   = 'N2'

if __name__ == "__main__":
    Parameter_List = (Pressure, Equivalence, Temperature, Diluent_Percentage)
    Mechanism_List = (Mechanism, Fuel, Oxidizer, Diluent)
    
    if Simulation_Taype == '0D':
        0D = 0DSimulation(Parameter_List, Mechanism_List)
        
    elif Simulation_Type =='1D':
        1D = 1DSimulation(Parameter_List, Mechanism_List)
        
    elif Simulation_Type = 'Burner':
        Burner = BurnerSimulation(Parameter_List, Mechanism_List)
    else:
        print('Error! Simulation_Type String Does Not Match!'+
              '\nMake sure string matches one of three options!')
        sys.exit()