#!/usr/bin/env python

""" Classes and functions used to calculate and analyze flame simulations
using cantera.

Created by Jeffrey Santner
Modified by Cody Ising 
"""

import os
import cantera
import numpy
# import matplotlib.pyplot as plt
import gc
from . import Classes  # File with custom classes and functions
# MOD_PATH = os.path.dirname(Classes.__file__)
# plt.style.use(os.path.join(MOD_PATH, 'custom_styles'))


def time_out(signum, frame):
    raise TimeOutError("Time out error")
# From http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call


class TimeOutError(Exception):
    "Timeout error"
    pass


class Flame(object):
    """
    Defines an object for a particular flame simulation.

    This object contains the flame simulation results as well
    as the initial conditions.
    Also contains modules for reporting results

    Parameters
    ----------
    mix : list or Classes.Mixture
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
            self.chemfile = Classes.find_chemfile(workingdir)
        else:
            self.chemfile = chemfile
        self.mixture = Classes.Mixture.check_or_create(mix, self.chemfile)
        self.P = P
        self.Tin = Tin
        self.workingdir = workingdir
        self.name = name
        self.error = None  # Store internal error codes
        self.flame_result = None

    def run(self, mingrid=200, loglevel=0, restart=None, mult_soret=False,
            highTbypass=False):
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
        highTbypass : bool
            Even if the temperature is > 500 K, use the low-T version.
        """
        if self.Tin > 500 and highTbypass is False:
            self.run_highT(mingrid, loglevel, restart, mult_soret)
            return

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
#                print('Main Loop 1')
                f.energy_enabled = True
                f.solve(loglevel=loglevel-1, refine_grid=True)
#                print('Main Loop 2')

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
                if 'Time out'in str(e):  # If it's a timeout message
                    self.flame_result = None
                    raise
                else:  # Except any other error, usually a cantera error
                    log(e, loglevel)
                    restart = None  # Stop using previous condition.
                    retry += 1
                    # Retry with different initial conditions.

        # After three tries, give up.
        self.flame_result = None
        self.error = True

    def run_highT(self, mingrid=200, loglevel=0, restart=None,
                  mult_soret=False):
        """
        Calculate flame speed with high initial temperature.

        Adjust the flame stabilization position to check that flame speed is
        independent of position. This ensures that the flame is not actually
        an ignition wave.

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

        # If restarting, the loop below will only run once. Otherwise, it will
        # run until the flame can be initialized.
        T_initial = self.Tin
        while T_initial > 300:
            log('Attempting with T = {:.0f}'.format(T_initial), loglevel)
            f = Flame(self.mixture, self.P, T_initial, self.workingdir,
                      self.chemfile, self.name)
            f.run(mingrid, loglevel, restart, False, True)
            if f.flame_result is not None:
                break
            else:
                T_initial -= 50
        if f.flame_result is None:
            raise TypeError('Could not solve for flame at low T')
        old_f = f.flame_result

        moved, refined= (1, 10)
        while not numpy.isclose(moved, refined, rtol=0.02):
            f = self._move_position(int(mingrid/2), old_f, loglevel)
            log('Location-independent flame speed found. Refining', loglevel)
            moved = f.u[0] * 100
            if mult_soret and f.transport_model != 'Multi':
                # Turn on multicomponent diffusion and Soret effect for final grid
                # refinement.
                f.transport_model = 'Multi'  # 'Mix' is default
                f.soret_enabled = True  # False is default
                f.solve(loglevel=loglevel-1)
                log('Adding multicomponent diffusion and Soret effect, Su = {:.2f} cm/s'.format(f.u[0] * 100), loglevel)
            self._grid_independence(f, mingrid, loglevel)
            old_f = f
            refined = f.u[0] * 100
#            log('After refining, Su = {:.2f} cm/s'.format(refined), loglevel)


        log('Finished calculation - S_u =  {:.2f} cm/s'.format(f.u[0] * 100),
            loglevel)
        if self.name is not None:  # Save to restart using f.restore()
            try:
                os.remove(
                        os.path.join(self.workingdir, self.name+'.xml'))
            except OSError:
                pass
            f.save(os.path.join(self.workingdir, self.name+'.xml'),
                   loglevel=loglevel-1)
        self.flame_result = f
        self.error = False
#        except Exception as e:  # Except all errors
#            if 'Time out' in str(e):  # If it's a timeout message
#                self.flame_result = None
#                raise
#            else:  # Except any other error, usually a cantera error
#                log(e, loglevel)
#                # Give up
#                self.flame_result = None
#                self.error = True

    def sensitivity(self):
        """ Calculate sensitivities of flame speed to reaction rates. """
        if self.flame_result is None:
            print('Flame cannot be plotted. It did not converge')
            return
        sens = self.flame_result.get_flame_speed_reaction_sensitivities()
        # Organize sensitivity
        self.sens = []
        for m in range(self.flame_result.gas.n_reactions):
            self.sens.append([m, sens[m],
                              self.flame_result.gas.reaction_equation(m)])

    # def plot(self, variables, save=False):
    #     """ Plot the profile of the inputted variables"""
    #     if self.flame_result is None:
    #         print('Flame cannot be plotted. It did not converge')
    #         return
    #     if type(variables) is not list:
    #         variables = [variables]
    #     fig, ax = plt.subplots()
    #     f = self.flame_result
    #     x = f.grid
    #     for var in variables:
    #         if var == 'T':
    #             ax.plot(x, f.T, '.', label='T [K]')
    #             if len(variables) == 1:
    #                 ax.set_ylabel('T [K]')
    #         elif var == 'u':
    #             ax.plot(x, f.u, '.', label='Velocity [m/s]')
    #             if len(variables) == 1:
    #                 ax.set_ylabel('Velocity [m/s]')
    #         else:
    #             ind = f.gas.species_index(var)
    #             ax.plot(x, f.X[ind], '.', label=var)
    #             if len(variables) == 1:
    #                 ax.set_ylabel(var)
    #     if len(variables) > 1:
    #         ax.legend()
    #     ax.set_xlabel('Position [m]')
    #     fig.tight_layout()
    #     if save:
    #         plt.savefig(os.path.join(self.workingdir, 'flame profile.png'))
    #         plt.close(fig)
    #     else:
    #         plt.show()

    def _check_too_close(self, fl):
        "Check that flame isn't too close to front boundary."
        for species, desired in self.mixture.fuel:
            found = fl.X[fl.gas.species_index(species), 0]
            if desired > 1e-7 or found > 1e-7:
                if abs(desired/found - 1) > 0.001:
                    raise RuntimeError(
                        'Flame is too close to the front of domain')
        if abs(fl.T[0]/self.Tin - 1) > 0.001:
                raise RuntimeError(
                    'Flame is too close to the front of domain')

    def _grid_independence(self, flame, mingrid, loglevel=0):
        """ Refine until flame is independent of the grid, and has at least
        mingrid points.
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
            # Check that flame isn't too close to front boundary
            self._check_too_close(flame)

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

    def _move_position(self, mingrid, old_f, loglevel=0):
        """ Move flame position to check that it is independent.
        Used in run_highT."""
        tol_ss = [1.0e-4, 1.0e-9]  # values from chemkin
        tol_ts = [1.0e-4, 1.0e-9]  # values from chemkin
        jac_age = 10  # Reuses of each jacobian before being recomputed
        time_step = 1e-6  # Time step in seconds
        num_time_steps = [2, 5, 10, 20, 50]  # Number of time steps to take

        transport = old_f.transport_model
        soret = old_f.soret_enabled

        new_speed = old_f.u[0]
        old_speed = new_speed * 10
        refine_criteria = old_f.get_refine_criteria()

        retry = -1
        # Flame locations to attempt on errors.
        locations = [0.06, 0.04, 0.02, 0.01, 0.005, 0.003, 0.001, 0.0005]

        counter = 0
        while abs(old_speed - new_speed)/old_speed > 0.01:
            counter += 1
            try:
                if self.error:  # Retry run
                    # Loose refinement criteria
                    refine_criteria = {'slope': 0.99, 'curve': 0.99,
                                       'prune': 0.01, 'ratio': 2}
                    time_step = time_step * 5.0  # Increase time step
                    log('Flame Calculation Error, adjusting refine criteria ' +
                        'and trying again', loglevel)

#                    for ind in range(1,len(initial_grid)):
#                        # Refine initial grid by inserting a point halfway
#                        # between each existing point
#                        initial_grid.insert(
#                            2*ind - 1,
#                            (initial_grid[2*ind-2] + initial_grid[2*ind-1])*0.5)

                    flame_loc_meters = locations[retry]  # Move flame forwards
                else:
                    flame_loc_ind = numpy.searchsorted(old_f.T,
                                                       (old_f.T[0] +
                                                        old_f.T[-1]) / 2)
                    flame_loc_meters_prev = old_f.grid[flame_loc_ind]
                    if counter == 1:
                        adjust = 1 / 0.85
                    elif self.error:
                        adjust = adjust ** 0.6
                    else:
                        adjust = 0.85
                    flame_loc_meters = flame_loc_meters_prev * adjust

                #if flame_loc_meters <= 0.002:
                #    flame_loc_meters = locations[retry]

                log('Moving flame from {:.3f} cm to {:.3f} cm'.format(
                    flame_loc_meters_prev*100, flame_loc_meters*100), loglevel)

                delta = flame_loc_meters_prev - flame_loc_meters
                ignore_ind = numpy.searchsorted(old_f.grid, delta)
                end_ind = numpy.searchsorted(old_f.grid, old_f.grid[-1] +
                                             delta)
                new_grid = ([0] + [x - delta for x in
                            old_f.grid[ignore_ind:end_ind]] + [old_f.grid[-1]])

                # Create new flame
                gas = self.mixture.cantera_mixture(self.Tin, self.P)

                f = cantera.FreeFlame(gas, new_grid)
                f.set_initial_guess()  # Call it now so that it doesn't get called later?

#                    locs = [x/old_f.grid[-1] for x in old_f.grid[ignore_ind:]] + [1.01]

                # Set the Temperature, velocity, and species profiles
                locs = [x/new_grid[-1] for x in new_grid]
                f.set_profile('u', locs, ([old_f.u[ignore_ind]] +
                                          list(old_f.u[ignore_ind:end_ind]) +
                                          [old_f.u[-1]]))
                f.set_profile('T', locs, ([f.inlet.T] +
                                          list(old_f.T[ignore_ind:end_ind]) +
                                          [old_f.T[-1]]))
                for n in range(f.gas.n_species):
                    f.set_profile(f.gas.species_name(n), locs,
                                  ([f.inlet.Y[n]] +
                                   list(old_f.Y[n, ignore_ind:end_ind]) +
                                   [old_f.Y[n, -1]]))
                f.set_fixed_temperature(0.5 * (old_f.T[0] + old_f.T[-1]))

                # Plot original profiles with new profiles
#                    fig, ax = plt.subplots(3, 1, sharex=True)
#                    ax[0].plot(old_f.grid, old_f.X[1,:], 'r.')
#                    ax[0].plot(f.grid, f.X[1, :], 'bx')
#
#                    ax[1].plot(old_f.grid, old_f.u, 'r.')
#                    ax[1].plot(f.grid, f.u, 'bx')
#
#                    ax[2].plot(old_f.grid, old_f.T, 'r.')
#                    ax[2].plot(f.grid, f.T, 'bx')

                # Keep these settings within the loop, so that
                # they can be changed on retries.
                f.flame.set_steady_tolerances(default=tol_ss)
                f.flame.set_transient_tolerances(default=tol_ts)
                f.set_max_jac_age(jac_age, jac_age)
                f.set_time_step(time_step, num_time_steps)
                f.set_refine_criteria(**refine_criteria)
                f.set_max_time_step(50)
                # Cantera would iterate forever
                # at the previous (very small) timestep

                # Set properties of the upstream fuel-air mixture.
                f.inlet.T = self.Tin
                f.inlet.X = self.mixture.reactants
                f.P = self.P * cantera.one_atm
                f.inlet.mdot = gas.density * old_f.u[0]

                f.energy_enabled = True
                f.transport_model = transport
                f.soret_enabled = soret
                f.solve(loglevel=loglevel-1, refine_grid=True)

                # Profiles after solving - Temporary
#                ax[0].plot(f.grid, f.X[1,:], 'g.')
#                ax[1].plot(f.grid, f.u, 'g.')
#                ax[2].plot(f.grid, f.T, 'g.')
#                plt.show()

                self._grid_independence(f, mingrid, loglevel)
#                log('At this position, S_u = {:.2f} cm/s with {} grid points'.format(f.u[0] * 100, f.flame.n_points),
#                    loglevel)
#                plt.plot(f.grid, f.T, 'bx')
#                plt.show()
                self.error = False
                old_speed = new_speed
                new_speed = f.u[0]
                old_f = f  # Previously, this was only done on restart runs. Why?
            except Exception as e:  # Except all errors
                self.error = True
                log(e, loglevel)
                retry += 1
                if retry == 6:
                    # Give up
                    log('Cannot solve this condition', loglevel)
                    self.flame_result = None
                    self.error = True
                    return
                if 'Time out' in str(e):
                    self.flame_result = None
                    raise
#                elif 'Flame is too close' in str(e):
#                else:
#                    self.error = True
#                    log(e, loglevel)
        return f


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
            refine_criteria['prune'] *= factor # Make sure pruning happens??
        max_criteria = max([refine_criteria[key] for key in
                            ['slope', 'curve', 'prune']])
        if max_criteria > 1.0:
            for key in ['slope', 'curve', 'prune']:
                refine_criteria[key] *= 1 / max_criteria

        fl.set_refine_criteria(**refine_criteria)
        return fl

    def open(self, file, loglevel=0):
        "Open flame solution, check that the conditions match."
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
        '''
        Open a file and use it or refine it.

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
        mix : list or Classes.Mixture
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
        '''
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


# def plot_from_xml(file, variables, show_fig=False):
#     """
#     Plot the profile of the inputted variables.

#     file: str
#         Path to an xml file containing a flame simulation
#     variables: list
#         List of variables to plot. Right now it only works for 'T' and 'u'
#     show_fig: bool
#         Either show the figure or save it to a file in the same directory
#         as file.
#     """
#     plt.style.use(r'C:\Users\jsantne\Documents\Research\Mild Ignition with SSG\From VM\Cantera Simulations\Cantera Codes\custom_styles')
#     var_info = {'T': 'T [K]', 'u': 'Velocity [m/s]'}

#     path, fname = os.path.split(file)
#     name, ext = os.path.splitext(fname)
#     if type(variables) is str:
#         variables = [variables]
#     chemfile = Classes.find_chemfile(path)
#     gas = cantera.Solution(chemfile)  # Dummy mixture
#     f = cantera.FreeFlame(gas)
#     f.restore(file, loglevel=0)  # Open the cantera flame object
#     x = f.grid
#     fig, axes = plt.subplots(nrows=len(variables), sharex=True, squeeze=False)
#     for var, ax in zip(variables, axes.flatten()):
#         y = getattr(f, var)  # Get the y-values to be plotted
#         ax.plot(x, y, '.')
#         if var in var_info:
#             ax.set_ylabel(var_info[var])
#         else:
#             ax.set_ylabel(var)
#     ax.set_xlabel('Position [m]')
#     if show_fig:
#         plt.show()
#     else:
#         plt.savefig(os.path.join(path, name + ' profile.png'))
#         plt.close(fig)


# def elevated_P_T_off_stoich(mix, P_targ, T_targ, workingdir, chemfile=None,
#                             name=None, P_start=1, T_start=300, phi_start=1,
#                             final_grid=500, loglevel=0):
#     """
#     Simulate a flame at 'difficult' conditions.

#     This function is used when a flame must be calculated
#     at elevated temperature and pressure at any equivalence ratio.

#     It first solves for a flame at P_start, T_start, phi_start,
#     and then slowly changes Pressure, temperature, and equivalence ratio
#     to approach P_targ, T_targ, and the original mix composition.

#     Parameters
#     ----------
#     mix : list or Classes.Mixture
#         If list, should be in the form [[name, moles],[...]].
#     P_targ : float
#         Target pressure for simulation in atm
#     T_targ : float
#         Target temperature for simulation in K
#     workingdir : str
#         location for temporary files
#     chemfile : str
#         Path to chemistry input file (cantera format)
#         If None (default), search for chemistry file in workingdir
#     name : str
#         Save xml file with this name. If a file with this name exists,
#         restart from it.
#     P_start : float
#         Pressure to use for first simulation (atm). This should be a pressure
#         where the simulation is "easy." Default is 1 atm.
#     T_start : float
#         Temperature to use for first simulation (K). This should be a
#         temperature where the simulation is "easy." Default is 300 K.
#     phi_start : float
#         Equivalence ratio (phi) to use for first simulation. This should be a
#         phi where the simulation is "easy." Default is 1.
#     final_grid : int
#             Minimum number of grid points in solution. Default is 500.
#     loglevel : int
#             Amount of diagnostic information. Default is 0.
#             0 produces no diagnostics, greater integers produce more.
#     """

# #    NOTE: Try Newton's method in iteration. Use derivative of state
# #    variables with respect to T and P to create initial guess for next T and P

#     if name:  # Try to restore from this file if it exists.
#         try:
#             flame = Flame.restore_or_run(
#                 mix, P_targ, T_targ, workingdir, name,
#                 name, loglevel-1, final_grid, chemfile)
#             log('Restoring or running from ' +
#                 os.path.join(workingdir, name+'.xml'), loglevel)
#             return flame
#         except IOError:
#             pass

#     if chemfile is None:
#         chemfile = Classes.find_chemfile(workingdir)

#     mixture = Classes.Mixture.check_or_create(mix, chemfile)
#     phi_targ = mixture.phi
#     if (abs(P_targ / P_start - 1) < 0.3 and abs(T_targ / T_start - 1) < 0.3
#         and abs(phi_targ - phi_start) < 0.2):# or T_targ > 500:
#         # Not at elevated P, T or off stoich. just run a normal flame.

#         log('P = {:5.2f} atm, T ={:6.1f} K, phi ={:5.2f}'.format(P_targ,
#             T_targ, mixture.phi), loglevel)

#         flame = Flame(mixture, P_targ, T_targ, workingdir, chemfile, name)
#         flame.run(final_grid, loglevel-1, mult_soret=True)
#         return flame

#     # First just try at given conditions
#     log('Attpemting given condition, P = {:5.2f} atm, T ={:6.1f} K, phi ={:5.2f}'.format(
#             P_targ, T_targ, mixture.phi), loglevel)
#     flame = Flame(mixture, P_targ, T_targ, workingdir, chemfile, name)
#     try:
#         flame.run(final_grid, loglevel-1, mult_soret=True)
#     except Exception as e:
#         log(str(e), loglevel)
#     else:
#         if flame.flame_result:
#             return flame

#     mixture.phi = phi_start
#     if loglevel > 0:
#         print('Starting from:')
#         print(mixture)
#         print('P = {:5.2f} atm, T ={:6.1f} K, phi ={:5.2f}'
#               .format(P_start, T_start, phi_start))

#     flame = Flame(mixture, P_start, T_start, workingdir, chemfile, name='temp')
#     flame.run(mingrid=100, loglevel=loglevel-1)
#     if not flame.flame_result:
#         log("Couldn't converge at T = {}, P = {}, phi={}" .format(T_start,
#             P_start, phi_start), loglevel)
#         return flame
#     step = 0.2
#     frac = 0.2
#     minstep = 0.005
#     count = 1
#     last_success = 0
#     speeds = [flame.flame_result.u[0]]
#     fracs = [0]
#     while step > minstep:
#         if frac > 1:
#             frac = 1
#         P = P_start*(P_targ/P_start)**frac  # Pressure changes exponentially.
#         T = (T_targ - T_start)*frac + T_start  # Temperature changes linearly.
#         phi = (phi_targ - phi_start)*frac + phi_start  # Phi changes linearly.
#         mixture.phi = phi  # Set mixture equivalence ratio to new phi.
#         log('P = {:5.2f} atm, T ={:6.1f} K, phi ={:5.2f}'.format(P, T, phi),
#             loglevel)
#         if count > 1:
#             # There are two previous solutions.
#             # Extrapolate to alter the initial velocity for the next run.
#             ratio = 1+((1-speeds[count-2]/speeds[count-1]) *
#                        (frac-fracs[count-1]) /
#                        (fracs[count-1]-fracs[count-2]))
#             f = cantera.FreeFlame(cantera.Solution(flame.chemfile))
#             f.restore(os.path.join(workingdir, 'temp.xml'),
#                       loglevel=loglevel - 1)
#             f.set_profile('u', f.grid/f.grid[-1], f.u*ratio)
# #            points = range(0, len(f.grid), 5) # Every 5th point.
# #            if len(a) - 1 not in points:
# #                points += [len(a)-1]
# #            g = [f.grid[x]/f.grid[-1] for x in points]
# #            new_u = [f.u[x]*ratio for x in points]
# #            f.set_profile('u', g, new_u) #sparser. Not necessary?
#             try:
#                 os.remove(os.path.join(workingdir, 'temp2.xml'))
#             except OSError:
#                 pass
#             # Save flame with adjusted velocity.
#             f.save(os.path.join(workingdir, 'temp2.xml'), loglevel=loglevel-1)
#             rest_name = 'temp2'
#         else:
#             rest_name = 'temp'
# #        if frac == 1:
# #            # Simulate target condition using refined grid, multicomponent
# #            # and Soret diffusion.
# #            flame = Flame(mix, P_targ, T_targ, workingdir, chemfile, name=name)
# #            flame.run(final_grid, loglevel, rest_name, True, highTbypass=True)
# #            os.remove(os.path.join(workingdir, rest_name+'.xml'))
# #        else:
#         flame = Flame(mixture, P, T, workingdir, chemfile, name='temp')
#         flame.run(mingrid=100, restart=rest_name, loglevel=loglevel,
#                   highTbypass=True)
#         if flame.flame_result:
#             if frac == 1:
#                 # Refine grid and add mult/soret
#                 flame = Flame(mix, P_targ, T_targ, workingdir, chemfile,
#                               name=name)
#                 flame.run(final_grid, loglevel, rest_name, True)
#                 os.remove(os.path.join(workingdir, rest_name+'.xml'))

#                 return flame
#             else:
#                 fracs.append(frac)
#                 speeds.append(flame.flame_result.u[0])
#                 last_success = frac
#                 frac += step
#                 count += 1
#         else:
#             step = step / 2
#             frac = last_success + step
#     flame.flame_result = None  # Simulation didn't succeed, remove result.
#     flame.name = name
#     return flame


# def sweep_P_T(mix, workingdir, pressures=None, temperatures=None,
#               P_T_list=None, mingrid=500, loglevel=0, fname='flame output.txt',
#               timeout=900, T_start=300, phi_start=1.0, fname_prefix='',
#               no_restarts=False):
#     """
#     Solve flames over a range of P, T with the same mixture.

#     Parameters
#     ----------
#     mix: list
#         list in the form [[name, moles],[...]] or Classes.mixture object
#     workingdir: str
#         Location to save output file, flame *.xml files,
#         and search for chemistry file
#     pressures: list
#         List of pressures (atm) to calculate flames
#     temperatures: list
#         List of pressures (K) to calculate flames
#     P_T_list: list
#         List of (pressure, temperature) tuples to calculate flames
#     loglevel: int
#         Integer determining verbosity of output
#     fname: str
#         Save results to this text file.
#     timeout: float
#         maximum time for a single simulation (seconds). If timeout=0, only load
#         simulations, don't run them.
#     T_start: float
#         to be used in elevated_P_T_off_stoich
#     phi_start: float
#         to be used in elevated_P_T_off_stoich
#     fname_prefix: str
#         string to prefix all flame filenames
#     no_restarts: Never restart any simulations from previous simulations.

#     Returns
#     -------
#     list of tuples in the form: (T_0 (K), P (atm), rho_0 (kg/m^3), T_f (K),
#                                  rho_f (g/cm^3), S_u (m/s), S_b (m/s),
#                                  max(dT/dx) (K/mm))
#     """
#     import signal
#     import parse

#     err_msg = 'sweep_P_T requires either (pressures, temperatures) OR P_T_list'
#     if numpy.all(P_T_list) and (numpy.all(temperatures) or numpy.all(pressures)):
#         raise SyntaxError(err_msg)
#     elif not (numpy.all(pressures) or numpy.all(temperatures) or numpy.all(P_T_list)):
#         raise SyntaxError(err_msg)
#     elif numpy.all(pressures) != numpy.all(temperatures):
#         raise SyntaxError(err_msg)
#     elif numpy.all(pressures) and numpy.all(temperatures):
#         import itertools
#         P_T_list = list(itertools.product(pressures, temperatures))

#     name_format = fname_prefix+'{:.1f} atm, {:.0f} K'
#     success_num = 0
#     output = []
#     for P, T in P_T_list:
#         name = name_format.format(P, T)

#         log('Beginning flame simulation for P={:.1f} atm, T_0={:.0f} K'.format(P, T),
#             loglevel)
#         error = True
#         if name + '.xml' in os.listdir(workingdir):
#             # Open old simulation
#             log('Opening old simulation', loglevel)
#             if timeout == 0:
#                 # Open old simulations and read their results, but do not
#                 # run new simulations.
#                 try:
#                     flame = Flame.restore_or_run(mix, P, T, workingdir, name,
#                                                  name, loglevel-1, mingrid,
#                                                  mult_soret=True,
#                                                  permissive=False)
#                 except RuntimeError:
#                     flame.flame_result = None
#             else:
#                 flame = Flame.restore_or_run(mix, P, T, workingdir, name, name,
#                                              loglevel-1, mingrid,
#                                              mult_soret=True)
#             error = False
#         elif no_restarts:
#             log('Starting flame simulation without restarting', loglevel)
#             flame = Flame(mix, P, T, workingdir, name=name)
#             if os.name == 'posix':
#                 # With Linux
#                 signal.signal(signal.SIGALRM, time_out)
#                 signal.alarm(timeout)  # Send alarm signal after timeout secs.
#             else:
#                 log('Timeout does not work with Windows', loglevel)
#             try:
#                 flame.run(mingrid, loglevel-1, mult_soret=True)
#                 error = False
#             except TimeOutError:
#                 log('Time Out Error', loglevel)
#                 error = True
#                 flame = Flame(mix, P, T, workingdir, name=name)
#                 flame.error = True
#             if os.name == 'posix':
#                 signal.alarm(0)  # Cancel the alarm.

#             if flame.flame_result is None:
#                 error = True
#             else:
#                 error = False
#         elif timeout == 0:
#             # Open old simulations and read their results, but do not
#             # run new simulations.
#             log('timeout=0, so new simulations are not performed',
#                 loglevel)
#             continue

#         # If there was an error loading a previous solution, or an
#         # error starting without using a previous solution, we will load a
#         # previous, similar solution to restart from.
#         elif error:
#             # Find successful flame calculations in the workingdir.
#             P_T_successes = []
#             for file in os.listdir(workingdir):
#                 name_old, ext = os.path.splitext(file)
#                 if ext == '.xml':
#                     result = parse.parse(fname_prefix+'{:g} atm, {:g} K',
#                                          name_old)
#                     if result:
#                         P_T_successes.append(result.fixed)

#             if not P_T_successes:  # This is the first simulation in the sweep.
#                 log('Starting flame simulation without restarting', loglevel)
#                 flame = elevated_P_T_off_stoich(
#                     mix, P, T, workingdir, name=name,
#                     final_grid=mingrid, loglevel=loglevel-1,
#                     T_start=T_start, phi_start=phi_start)
#                 error = False

#             else:  # Look for a successful simulation to load or restart from.
#                 #goodness=[abs(P_s/P-1)+4*abs(T_s/T-1) for P_s, T_s in P_T_successes] #weight temperature by 4, seems to help

#                 # Factor of 2 in pressure similar to 100 K change in temperature
#                 goodness = [abs(numpy.log10(P_s/P)) + 2*abs(T_s-T)/333
#                             for P_s, T_s in P_T_successes]

#                 # Index of previous successful condition
#                 # that is closest to the present condition
#                 g, ind = min((val, index) for (index, val) in
#                              enumerate(goodness))
#                 rest_P, rest_T = P_T_successes[ind]
#                 rest_name = name_format.format(rest_P, rest_T)

#                 if os.name == 'posix':
#                     signal.signal(signal.SIGALRM, time_out)
#                     signal.alarm(timeout)  # Send signal after timeout secs.
#                 try:
#                     flame = Flame.restore_or_run(mix, P, T, workingdir, name,
#                                                  rest_name, loglevel-1,
#                                                  mingrid, mult_soret=True)
#                     error = False
#                 except TimeOutError:
#                     log('Time Out Error', loglevel)
#                     error = True
#                     flame = Flame(mix, P, T, workingdir, name=name)
#                     flame.error = True
#                 if os.name == 'posix':
#                     signal.alarm(0)  # Cancel the alarm.

#                 #remove temporarily
# #                if flame.error: #there was a cantera error while restoring
# #                    error=False
# #                    if loglevel>0:
# #                        print('Restarting from previous condition failed. Attempting to start from scratch for '+name)
# #                    signal.signal(signal.SIGALRM, time_out)
# #                    signal.alarm(timeout) #send alarm signal after this many seconds.
# #                    try:
# #                        flame=elevated_P_T_off_stoich(mix, P, T, workingdir, name=name, final_grid=mingrid, loglevel=loglevel-1, T_start=T_start, phi_start=phi_start) #try to run this condition without restarting from a previous, good condition
# #                    except TimeOutError:
# #                        if loglevel > 0:
# #                            print('Time Out Error')
# #                        error=True
# #                    signal.alarm(0) #cancel the alarm

#         if flame.flame_result and not error:
#             flame.save(os.path.join(workingdir, fname), save_override=False,
#                        loglevel=loglevel)
# #            P_T_successes.append((P,T))
#             f = flame.flame_result
#             output.append((f.inlet.T, f.P/cantera.one_atm,
#                            f.density_mass[0], f.T[-1],
#                            f.density_mass[-1], f.u[0], f.u[-1],
#                            max(abs(numpy.diff(f.T)/numpy.diff(f.grid)))/1000))
#             success_num += 1
#             log('\tSuccess', loglevel)
#         else:
#             log('\tFailure', loglevel)
#     log('Failure for {:3.1f}% of cases: '.format(
#             100*(1-float(success_num)/float(len(P_T_list)))), loglevel)
# #        print(sorted(list(set(P_T_list)-set(P_T_successes))))
#     return output


# def sens_plot(mix, P, T, workingdir, nrxns=15, chemfile=None, name=None,
#               mingrid=200, loglevel=0, restart=None, mult_soret=False):
#     """ Create a plot of the sensitivity of flame speed to the top nrxn
#     reactions.

#     All variable names (except nrxn) as in Flame and Flame.run
#     Parameters
#     ----------
#     mix : list or Classes.Mixture
#         If list, should be in the form [[name, moles],[...]].
#     P : float
#         Unburned pressure (atm)
#     T : float
#         Temperature (Kelvin) of the unburned gas
#     workingdir : str
#         location for temporary files and output figure
#     nrxns : int
#         Number of reactions to plot
#     chemfile : str
#         Path to chemistry input file (cantera format)
#         If None, search for chemistry file in workingdir
#     name : str
#         String to identify a simulation.
#     mingrid : int
#             Minimum number of grid points in solution.
#     loglevel : int
#         Amount of diagnostic information.
#         0 produces no diagnostics, greater integers produce more.
#     restart : str
#         Name of simulation to restart from.
#         This is not a file path - an xml file in the workingdir
#         with the name "restart.xml" will be used.
#     mult_soret : bool
#         Whether or not to use multicomponent diffusion and
#         include the Soret effect.
#     """
#     fl = Flame(mix, P, T, workingdir, chemfile, name)
#     fl.run(mingrid, loglevel, restart, mult_soret)
#     fl.sensitivity()
#     fig, ax = plt.subplots()
#     ax.grid(axis='x', which='major', ls='--')
#     ax.grid(axis='y', which='minor', c='k')
#     sens = fl.sens
#     sens.sort(key=lambda x: abs(x[1]), reverse=True)
#     sens_plot = sens[:nrxns]
#     ylocs = numpy.arange(nrxns)
#     ax.barh(ylocs, [x[1] for x in sens_plot], align='center')
#     ax.set_yticks(ylocs)
#     ax.set_yticklabels([x[2] for x in sens_plot])
#     ax.set_yticks(ylocs - 0.5, minor=True)
# #    ax.tick_params(axis='y', which='minor', bottom='off')
#     ax.invert_yaxis()
#     ax.axvline(c='k')
#     ax.set_xlabel('Normalized Sensitivity')
#     ax.set_ylim([max(ylocs)+0.5, min(ylocs)-0.5])
#     fig.tight_layout()
#     fig.savefig(os.path.join(workingdir, 'Flame Sensitivity.png'))
#     plt.close(fig)


# def sweep_P_T_plot(T_P_su, outfile=None, cbar_label=None):
#     print('sweep_P_T_plot is no longer available. ' +
#           'Use plot_tools.contour_plot instead')


def log(msg, level):
    if level > 0:
        print(msg)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print('\n\nERROR: This function must be called with a directory as an argument\n\n')
        sys.exit()

    workingdir = sys.argv[1]
    execfile(os.path.join(workingdir, 'input.py'))
    flame=Flame(fuel, P, Tin, workingdir, chemfile=chemfile, name=name)
    flame.run(mingrid=mingrid, loglevel=loglevel, restart=restart)
    flame.save(os.path.join(workingdir, 'output.txt'))
