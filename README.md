# design-of-experiments
Optimally design flow reactor experimental conditions with Sensiztized_0D_Experiment.py
Optimally design flame experimental conditions with Sensitized_Flame_Experiment.py

## Simulation Initializer Usage
Open Simulation_Initializer.py
-User shall provide simulation type to be run to Simulation_Type
--Current build has two options (0D, 1D)
---0D will run zero-dimensional simulation
---1D will run one-dimensional simulation

-User shall provide mixture type to be used to Mixture_type
--Current build has 8 options (phi_oxi/dil, phi_fuel/dil, phi_fuel, phi_oxi, phi_dil, fuel_dil, oxi_dil, oxi_fuel)
---phi_oxi/dil specifies the equivalence ration and oxidizer to diluent ratio
---phi_fuel/dil specifies the equivalence ration and fuel to diluent ratio
---phi_fuel specifies the equivalence ratio and fuel mole fraction
---phi_oxi specifies the equivalence ratio and oxidizer mole fraction
---phi_dil specifies the equivalence ratio and diluent mole fraction
---fuel_dil specifies the fuel and diluent mole fractions
---oxi_dil specifies the oxidizer and diluent mole fractions
---oxi_fuel specifies the oxidizer and fuel mole fractions

-User shall provide a list of parameters to be used in creating thermodynamic and mixture properties
--All parameters are lists of three numbers
---First Number: Initial point
---Second Number: End Point
---Third Number: Number of points. If set to 1 only first point is evaluated
--Units
---Pressure: Atmosphere [atm]
---Temperature: Kelvin [K]
---Mixture parameters. Not all of these will be used, depending on Mixture_type
----Equivalence: Dimensionless (<1:Fuel Lean, 1:Unity, >1:Fuel Rich)
----fraction_in_oxidizer_or_fuel: Percentage of Total Fuel or Oxidizer [%]
----fuel_fraction: Percentage of Total Mixture [%]
----oxidizer_fraction: Percentage of Total Mixture [%]

-User shall provide the type of scale for the array of parameters created.
--Array_type current options ('log', 'lin')

-User shall provide the Fuel, Oxidizer, and Diluent Chemical Formulas
--Fuel and Oxidizer can either be a single chemical string or multichemical list
---For multichemical follow structure ['Chemical1', % of Total, ...]
---The percentage following the chemical name should sum up to 1
--Diluent should be a single chemical added to the mixture

-User shall provide Zero-Dimensional conditions
--SpecificSpecies is a list of species of interest to collective senstivite against each reaction and time step
--Starttime is an integer used in the case a reading is too early
--Endtime is a float that defines the end time of the simulation
--Delta_T is an integer that checks to make sure the temperature rise is not great enought to enduce a flame.
--threshold is a float that is determined by the equipment minimum ppm reading

-User shall provide One-DImensional conditions (Flame)
--Mingrid is an integer that defines the number of points that will be solved for
--Mul_soret is a boolean. If true increase the accuracy of the simulation by calculating multicomponent diffusion and soret effect.
--Loglevel is a integer from 0 to 10. The greater the number the more information provided from the simulation.

-User shall provide boolean statements.
--Save_files if True will save results from simulation used in plotting
--Save_Time only for Zero-Dimensional simulation, if True saves mole_fraction data used in GUI

## Plotting Scripts
There are two plotting scripts avalaible. A 0D and 1D plotting script.

For Zero-Dimensional plotting open: Sensitized_0D_Plotting.py
-User shall provide reactions of interest in Rxn_interest list.
--If Rxn_interest is left blank no plots of specific reactions will be created.
-User shall provide minimum threshold of sensitivitiy to Threshold.
--If the sensitivity of a reaction is greater than, or equal, to Threshold data will be plotted for that case.

For One-Dimensional plotting open: Sensitized_Flame_Plotting.py
-User shall provide reactions of interest in Rxn_interest list.
--If Rxn_interest is left blank no plots of specific reactions will be created.
-User shall provide four strings to list Four_Plot.
--Options are 'T', 'F', 'Phi', 'O2', 'Su', 'P', or any species name.
--If left blank options 'P', 'F', 'Phi', 'O2' will be used
-User shall provide minimum flame speed in m/s to Min_speed.
--Default Min_speed is set to 0
--If the flame speed is above the specified value data will be plotted for that case
-User shall provide number of top reactions to normalize data to Nrxns.
-- Nrxns will define how many of the top reactions in decending order that are used to normalize sensitivities.
-User shall provide minimum threshold of sensitivitiy to Threshold.
--If the sensitivity of a reaction is greater than, or equal, to Threshold data will be plotted for that case.

Run Plotting Scripts
-After setting up variables for plotting script user is prompted with a message
-'Please type name of folder, 1, or 2.'
-' 1 or blank: Use the last folder that was analyzed.'
-' 2: Use the last folder that was simulated:'
-User shall provide the name of the folder, or the number requested, to pull data from simulation folder.

## License
[CSULA] (https://www.calstatela.edu/faculty/jeff-santner)