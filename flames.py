#!/usr/bin/env python

"""Functions used to calculate and analyze flame simulations using cantera.

Created by Jeffrey Santner
Modified by Cody Ising
"""

import os
import cantera
import gc
import glob
from cantera import ck2cti


class Flame(object):
    """
    Defines an object for a particular flame simulation.

    This object contains the flame simulation results as well
    as the initial conditions.
    Also contains modules for reporting results

    Parameters
    ----------
    mix : list or Mixture
        If list, should be in the form [[name, moles],[...]].
    P : float
        Unburned pressure (atm)
    Tin : float
        Temperature (Kelvin) of the unburned gas
    workingdir : str
        location for temporary files
    chemfile : str
        Path to chemistry input file (cantera format)
        If None, search for chemistry file in workingdir
    name : str
        String to identify a simulation.
    """

    def __init__(self, mix, P, Tin, workingdir, chemfile=None, name=None):
        if not chemfile:
            self.chemfile = find_chemfile(workingdir)
        else:
            self.chemfile = chemfile
        self.mixture = Mixture.check_or_create(mix, self.chemfile)
        self.P = P
        self.Tin = Tin
        self.workingdir = workingdir
        self.name = name
        self.error = None  # Store internal error codes
        self.flame_result = None

    def run(self, mingrid=200, loglevel=0, restart=None, mult_soret=False):
        """Run a flame simulation using cantera.

        Much of this code was originally taken from
        cantera/examples/onedim/adiabatic_flame.py

        mingrid : int
            Minimum number of grid points in solution.
        loglevel : int
            Amount of diagnostic information.
            0 produces no diagnostics, greater integers produce more.
        restart : str
            Name of simulation to restart from.
            This is not a file path - an xml file in the workingdir
            with the name "restart.xml" will be used.
        mult_soret : bool
            Whether or not to use multicomponent diffusion and
            include the Soret effect.
        """
        self.needs_saving = True

        try:
            os.remove(self.chemfile[:-4]+'.xml')
        except OSError:
            pass

        reactants = self.mixture.reactants
        gas = self.mixture.cantera_mixture(self.Tin, self.P)

        # Initial conditions and convergence parameters
        initial_grid = [0, 0.02, 0.04, 0.06, 0.08, 0.5]
        # initial_grid in meters. One point at 40% of the domain is included
        # regardless because that is the default flame location.
        f = cantera.FreeFlame(gas, initial_grid)  # Create flame object

        tol_ss = [1.0e-4, 1.0e-9]  # values from chemkin
        tol_ts = [1.0e-4, 1.0e-9]  # values from chemkin
        jac_age = 10  # Reuses of each jacobian before being recomputed
        time_step = 1e-6  # Time step in seconds
        num_time_steps = [2, 5, 10, 20, 50]  # Number of time steps to take
        refine_criteria = {'slope': 0.85, 'curve': 0.99, 'prune': 0.01,
                           'ratio': 2}
#        f.show_solution() #debug

        if restart:
            log('Restoring from ' + restart, loglevel)
            f.restore(os.path.join(self.workingdir, restart+'.xml'),
                      loglevel=loglevel-1)
            refine_criteria = f.get_refine_criteria()

        retry = 0
        while retry < 3:
            try:
                if retry > 0:
                    # Loose refinement criteria
                    refine_criteria = {'slope': 0.99, 'curve': 0.99,
                                       'prune': 0.01, 'ratio': 2}
                    time_step = time_step * 5.0  # Increase time step
                    log('Flame Calculation Error, adjusting initial ' +
                        'conditions and trying again', loglevel)
                    for ind in range(1, len(initial_grid)):
                        # Refine initial grid by inserting a point halfway
                        # between each existing point
                        initial_grid.insert(
                            2*ind - 1,
                            (initial_grid[2*ind-2] +
                             initial_grid[2*ind-1])*0.5)
                    f = cantera.FreeFlame(gas, initial_grid)

                # Keep these settings within the loop, so that
                # they can be changed on retries.
#                print('Main Loop')
                f.flame.set_steady_tolerances(default=tol_ss)
                f.flame.set_transient_tolerances(default=tol_ts)
                f.set_max_jac_age(jac_age, jac_age)
                f.set_time_step(time_step, num_time_steps)
                f.set_refine_criteria(**refine_criteria)
                f.set_max_time_step(50)
                # Cantera would iterate forever
                # at the previous (very small) timestep

                # Set properties of the upstream fuel-air mixture.
                # Must be in loop in case of restore from different mixture
                f.inlet.T = self.Tin
                f.inlet.X = reactants
                f.P = self.P * cantera.one_atm

                if not restart:
                    # Solve with the energy equation disabled for stability
                    f.energy_enabled = False
                    f.solve(loglevel=loglevel-1, refine_grid=True)
                # Enable the energy equation for the following loop
                f.energy_enabled = True
                f.solve(loglevel=loglevel-1, refine_grid=True)

                if mult_soret:
                    f.transport_model = 'Multi'  # 'Mix' is default
                    f.soret_enabled = True  # False is default
                # Refine the grid and check for grid independence.
                self._grid_independence(f, mingrid, loglevel)
                log('Finished calculation - S_u =  {:.2f} cm/s'.format(f.u[0] * 100),
                    loglevel)
                if self.name:
                    try:
                        os.remove(
                            os.path.join(self.workingdir, self.name+'.xml'))
                    except OSError:
                        pass
                    # Save solution to restart using f.restore()
                    f.save(os.path.join(self.workingdir, self.name+'.xml'),
                           loglevel=loglevel-1)
                self.flame_result = f
                return
            except Exception as e:  # Except all errors
                log(e, loglevel)
                restart = None  # Stop using previous condition.
                retry += 1
                # Retry with different initial conditions.

        # After three tries, give up.
        self.flame_result = None
        self.error = True

    def sensitivity(self):
        """Calculate sensitivities of flame speed to reaction rates."""
        if self.flame_result is None:
            print('Flame cannot be plotted. It did not converge')
            return
        sens = self.flame_result.get_flame_speed_reaction_sensitivities()
        # Organize sensitivity
        self.sens = []
        for m in range(self.flame_result.gas.n_reactions):
            self.sens.append([m, sens[m],
                              self.flame_result.gas.reaction_equation(m)])

    def _grid_independence(self, flame, mingrid, loglevel=0):
        """Refine until grid-independence and >= mingrid points.

        flame:
            cantera.FreeFlame object
        mingrid:
            minimum number of grid points
        loglevel:
            integer specifying amount of information to print.
        """
        grid = flame.flame.n_points
        speed = flame.u[0] / 10  # Make the while loop below run at least once.
        while grid < mingrid or abs(speed-flame.u[0])/speed > 0.05:
            speed = flame.u[0]  # save previous speed
            flame = self._refine(flame, mingrid)  # Adjust refinement params
            msg = ('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n' +
                   '{} grid points, {} needed\n' +
                   'refining to slope = {:.3f}, curve = {:.3f}, prune = {:.4f}\n' +
                   '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
            refine_criteria = flame.get_refine_criteria()
            log(msg.format(flame.flame.n_points, mingrid,
                           refine_criteria['slope'], refine_criteria['curve'],
                           refine_criteria['prune'], ), loglevel - 1)
            flame.solve(loglevel=loglevel - 1, refine_grid=True)

            grid = flame.flame.n_points  # Final number of points

            log('Grid independence? Su = {:.2f} cm/s with {:} points'.format(
                    flame.u[0]*100, grid), loglevel)

    def _refine(self, fl, mingrid):
        """
        Adjust refinement criteria to get to mingrid points.

        Returns cantera.FreeFlame object
        """
        refine_criteria = fl.get_refine_criteria()
        grid = fl.flame.n_points

        # For example, if you want to double the number of grid
        # points, divide grad and curv by two.
        # If factor > 3, or grid is approaching max (1000),
        # loosen criteria to remove points
        factor = grid / mingrid
        if 2 > factor > 0.7 and grid < 900:
            factor = 0.7  # Do some refining
        elif factor < 0.1:
            factor = 0.1  # Don't refine too much all at once

        # Refine the mesh criteria.
        for key in ['slope', 'curve', 'prune']:
            refine_criteria[key] *= factor
        if factor > 1:
            refine_criteria['prune'] *= factor  # Make sure pruning happens??
        max_criteria = max([refine_criteria[key] for key in
                            ['slope', 'curve', 'prune']])
        if max_criteria > 1.0:
            for key in ['slope', 'curve', 'prune']:
                refine_criteria[key] *= 1 / max_criteria

        fl.set_refine_criteria(**refine_criteria)
        return fl

    def open(self, file, loglevel=0):
        """Open flame solution, check that the conditions match."""
        gas = cantera.Solution(self.chemfile)
        f = cantera.FreeFlame(gas)  # Create flame object.
        f.restore(file, loglevel=loglevel)

        # Check that the file opened has matching P, T, and mixture.
        if (
                abs(f.T[0]/self.Tin-1) > 0.0001 or  # temperature mismatch
                abs(f.P/cantera.one_atm/self.P - 1) > 0.0001):  # pressure mismatch
            raise RuntimeError(
                'Error opening flame object. ' +
                'Requested T={}, P={}. Found T={:5.3f}, P={:5.3f}'.format(
                self.Tin, self.P, f.T[0], f.P/cantera.one_atm))
        for species, desired in self.mixture.fuel:
            found = f.X[f.gas.species_index(species), 0]
            if desired > 1e-7 or found > 1e-7:
                if abs(desired/found - 1) > 0.001:
                    raise RuntimeError(
                        'Error opening flame object. ' +
                        'Requested {} mole fraction = {:6.4f}. Found mole fraction of {:6.4f}'.format(
                        species, desired, found))
        log('Loading flame solution from ' + file, loglevel)
        self.flame_result = f

    def save(self, outfile, save_override=False, loglevel=0):
        """Save global results to text file."""
        if not self.flame_result:
            print('Flame result not found')
            return
        if not self.needs_saving and not save_override:
            return
        f = self.flame_result
        categories = ['Name', 'T_0 (K)', 'P (atm)', 'rho_0 (kg/m^3)',
                      'T_f (K)', 'rho_f (g/cm^3)', 'S_u (cm/s)', 'S_b (cm/s)']

        out = [self.name, f.inlet.T, f.P/cantera.one_atm,
               f.density_mass[0], f.T[-1], f.density_mass[-1],
               f.u[0]*100, f.u[-1]*100]
        if os.path.isfile(outfile):
            file = open(outfile, 'a+')

        else:
            file = open(outfile, 'w+')
            file.write('\t'.join(categories)+'\n')
        file.write('\t'.join(['{}'.format(x) for x in out])+'\n')
        file.close()
        log('Saved flame result to ' + outfile, loglevel)

    def __repr__(self):
        return('Flame({}, {}, {}, {})'.format(self.mixture,
               self.P, self.Tin, self.workingdir))

    @classmethod
    def restore_or_run(cls, mix, P, Tin, workingdir, name, rest_name,
                       loglevel=0, mingrid=200, chemfile=None,
                       mult_soret=False, permissive=True):
        """Open a file and use it or refine it.

        Open the file specified by the workingdir and rest_name.
        If it's sufficient, return that flame object.
        If the grid needs to be improved or the P, T, mix is wrong,
        restart from it.

        NOTE: There is some kind of memory leak in this function, don't use it
        with large mechanisms. It seems that if the previous solution raises an
        error in flame.open (incorrect T, P, or mix) it isn't removed from
        memory before attempting to run the flame in the except block. So,
        a memory overload is possible.

        Parameters
        ----------
        mix : list or Mixture
            If list, should be in the form [[name, moles],[...]].
        P : float
            Unburned pressure (atm)
        Tin : float
            Temperature (Kelvin) of the unburned gas
        workingdir : str
            location for temporary files
        name : str
            String to identify a simulation.
        rest_name : str
            This is the name of the simulation to open and attempt to use.
        loglevel : int
            Amount of diagnostic information.
            0 produces no diagnostics, greater integers produce more.
        mingrid : int
            Minimum number of grid points in solution.
        chemfile : str
            Path to chemistry input file (cantera format)
            If None, search for chemistry file in workingdir
        mult_soret : bool
            Whether or not to use multicomponent diffusion and
            include the Soret effect.
        permissive : bool
            Default behavior is True. If True, rest_name will be used as the
            initial guess to simulate the desired conditions, even if rest_name
            does not have the correct initial P, T, mixture, or number of grid
            points. If False, this code will not attempt to run any
            simulations, it will just return the flame object if it is
            appropriate, or raise an error

        Returns
        -------
        Returns a Flame object
        """
        f = os.path.join(workingdir, rest_name+'.xml')
        if not os.path.isfile(f):
            raise IOError('File not found: ' + f)
        flame = cls(mix, P, Tin, workingdir, chemfile, name)
        try:  # Use previous solution without re-calculating
            flame.open(f, loglevel-1)
            if flame.flame_result.flame.n_points < mingrid:
                raise RuntimeError('Not resolved enough. {:} vs. {:} points'.
                                   format(flame.flame_result.flame.n_points,
                                          mingrid))
            if name == rest_name:
                flame.needs_saving = False  # Flame has not been modified
            else:
                flame.needs_saving = True
        except RuntimeError as e:
            # Previous calculation had
            # different P, T, or mixture. Or, not enough grid points.
            # Restart from it.
            del flame  # Fix memory issue? No, doesn't help
            gc.collect()
            log(e, loglevel)
            if permissive:
                flame = cls(mix, P, Tin, workingdir, chemfile, name)
                flame.run(mingrid, loglevel, rest_name, mult_soret)
            else:
                raise

        return(flame)


class Mixture(object):
    """Defines a mixture of species.

    Parameters
    ----------
    fuel : list
        Define the mixture composition. List of the form
        [['species name', mole fraction], ['species name', mole fraction], ...]
    chemfile : str
        Path to a chemistry file.
    """

    def __init__(self, fuel, chemfile):
        self.fuel = fuel
        self._normalize()
        self.chemfile = chemfile
        self.cantera_mix = None  # initialize variables
        self._phi = None  # initialize variables
        self._fuels = None  # Initialize
        self._oxidizer = None  # Initialize

    def _normalize(self):
        """Normalize mole fractions to sum to one."""
        total = sum([x[1] for x in self.fuel])
        self.fuel = [[row[0], row[1]/total] for row in self.fuel]
        reactant_pieces = []
        for r in self.fuel:
            reactant_pieces.append(r[0] + ':' + str(r[1]))
        self.reactants = ', '.join(reactant_pieces)

    def cantera_mixture(self, T, P):
        """Return a cantera solution object.

        T: temperature in Kelvin
        P: Pressure in atm
        """
        self._normalize()
        if self.chemfile[-4:] != '.cti':
            raise RuntimeError('Chemfile must be in cantera format')
        try:
            gas = cantera.Solution(self.chemfile)
        except RuntimeError as err:  # chemfile has been deleted
            print(err)
            return None

        gas.TPX = T, P*cantera.one_atm, self.reactants
        self.cantera_mix = gas
        return gas

    def _equivalence_ratio(self):
        if self.cantera_mix is None:
            gas = self.cantera_mixture(300, 1)
            if gas is None:  # chemfile is missing
                self._phi = None
                return
        else:
            gas = self.cantera_mix  # Don't need to reload the kinetic model
        self._phi = self.cantera_mix.get_equivalence_ratio()

    def _components(self):
        """Calculate the components of the mixture."""
        if self.cantera_mix is None:
            gas = self.cantera_mixture(300, 1)
            if gas is None:  # chemfile is missing
                return
        else:
            gas = self.cantera_mix  # Don't need to reload the kinetic model
        fuels = []
        oxidizers = []
        for row in self.fuel:
            species = row[0]
            C = gas.n_atoms(species, 'C')
            H = gas.n_atoms(species, 'H')
            O = gas.n_atoms(species, 'O')
            if C == 1 and O == 2 and H == 0:  # This is CO2
                oxidizers.append(species)
            elif C == 0 and O == 1 and H == 2:  # This is H2O
                oxidizers.append(species)
            elif C == 0 and O == 2 and H == 0:  # This is O2
                self.oxygen = species  # Unneccesary?
                oxidizers.append(species)
            elif C + H + O == 0:  # This is inert
                oxidizers.append(species)
            else:
                fuels.append(species)  # This is a fuel component

        self._fuels = fuels
        self._oxidizer = oxidizers

    def _change_phi(self, phi_targ=1.0):
        """Change the equivalence ratio.

        Change the equivalence ratio by changing fuel content, holding
        diluent/o2 ratio constant. This could also be done using
        Solution.set_equivalence_ratio in cantera
        """
        fuel_comps = self.fuels
        for row in self.fuel:
            if row[0] in fuel_comps:
                row[1] *= phi_targ / self.phi
        self._normalize()
        self.cantera_mix = None  # Force re-calculation
        self._phi = None
        if abs(self.phi - phi_targ) > 1e-5:
            print('Error: Could not correctly change equivalence ratio.')

    def __str__(self):
        outstr = ''
        for row in self.fuel:
            outstr += row[0] + ':\t'+'{:5.3f}'.format(row[1])+'\n'
        return outstr

    def __repr__(self):
        return('Mixture({}, {})'.format(self.fuel, self.chemfile))

    @property
    def phi(self):
        """Equivalence ratio property that can use a getter and setter."""
        if self._phi is None:
            self._equivalence_ratio()
        return self._phi

    @phi.setter
    def phi(self, phi_targ):
        self._change_phi(phi_targ)

    @property
    def fuels(self):
        """Find the fuel components."""
        if self._fuels is None:
            self._components()
        return self._fuels

    @property
    def oxidizer(self):
        """Find the mixture components."""
        if self._oxidizer is None:
            self._components()
        return self._oxidizer

    @classmethod
    def check_or_create(cls, mix, chemfile):
        """
        If mix is a list, create a mixture object
        If mix is already a mixture object, check that the chemfile is correct
        """
        if type(mix) is list:
            return(cls(mix, chemfile))
        elif type(mix) is cls:  # mix is already a member of this class
            if mix.chemfile == chemfile:
                return(mix)
            else:
                raise NameError('The chemistry file initially used for the mixture is not the same as the desired chemistry file. {} vs. {}'.format(mix.chemfile, chemfile))
        else:
            raise RuntimeError('The mixture is not a list or a mixture object'+
                               ', it is type: ' + str(type(mix)))


def find_chemfile(workingdir, verbose=0):
    """Find the chemistry file in this working directory.

    If the chemistry file is in chemkin format, convert it to cantera
    Preference order is: chem.cti, chem.inp, *.cti, *.inp

    Returns the path to a cantera chemistry file in the workingdir
    """
    # Search for chemistry file.
    if os.path.isfile(os.path.join(workingdir, 'chem.cti')):
        if verbose > 0:
            print('Chemfile is ' + os.path.join(workingdir, 'chem.cti'))
        return(os.path.join(workingdir, 'chem.cti'))
    elif os.path.isfile(os.path.join(workingdir, 'chem.inp')):
        chemfile = os.path.join(workingdir, 'chem.inp')
    elif glob.glob(os.path.join(workingdir, '*.cti')):
        return(glob.glob(os.path.join(workingdir, '*.cti'))[0])
    elif glob.glob(os.path.join(workingdir, '*.inp')):
        chemfile = glob.glob(os.path.join(workingdir, '*.inp'))
    else:
        raise OSError('No chemistry file found in '+workingdir)

    # This section runs if a chemkin format file is used.
    path, file = os.path.split(chemfile)
    modelname, ext = os.path.splitext(file)
    convert_mech(chemfile)
    if verbose:
        print('Converting chemistry to cantera format')
    return(os.path.join(path, modelname+'.cti'))


def convert_mech(chemfile):
    """Convert chemkin chemistry file to cantera chemistry file."""
    path, file = os.path.split(chemfile)
    modelname, ext = os.path.splitext(file)

    thermo = None
    trans = None
    if os.path.isfile(os.path.join(path, 'therm.dat')):
        thermo = os.path.join(path, 'therm.dat')
    if os.path.isfile(os.path.join(path, 'tran.dat')):
        trans = os.path.join(path, 'tran.dat')

    arg = ['--input=' + chemfile]
    if thermo is not None:
        arg.append('--thermo='+thermo)
    if trans is not None:
        arg.append('--transport='+trans)
    arg.append('--permissive')
    ck2cti.main(arg)


def log(msg, level):
    if level > 0:
        print(msg)
