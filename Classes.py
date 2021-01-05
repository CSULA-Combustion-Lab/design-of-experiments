# Modified for use in python 3

import os
import glob
import cantera
from cantera import ck2cti


class Mixture(object):
    """Defines a mixture of species

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
        self.alpha = None  # initialize variables

    def _normalize(self):
        "Normalize mole fractions to sum to one."
        total = sum([x[1] for x in self.fuel])
        self.fuel = [[row[0], row[1]/total] for row in self.fuel]
        reactant_pieces = []
        for r in self.fuel:
            reactant_pieces.append(r[0] + ':' + str(r[1]))
        self.reactants = ', '.join(reactant_pieces)

    def cantera_mixture(self, T, P):
        '''
        Return a cantera solution object.
        T: temperature in Kelvin
        P: Pressure in atm
        '''
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
                self.alpha = None
                return
        else:
            gas = self.cantera_mix  # Don't need to reload the kinetic model
        self._phi = self.cantera_mix.get_equivalence_ratio()

    def _components(self):
        "Calculate the components of the mixture"
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
        """Change the equivalence ratio by changing fuel content,
        holding diluent/o2 ratio constant.
        This could also be done using Solution.set_equivalence_ratio in cantera
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
        "This is called if you print(mixture)"
        outstr = ''
        for row in self.fuel:
            outstr += row[0] + ':\t'+'{:5.3f}'.format(row[1])+'\n'
        return outstr

    def __repr__(self):
        return('Mixture({}, {})'.format(self.fuel, self.chemfile))

    @property
    def phi(self):
        '''Equivalence ratio property that can use a getter and setter'''
        if self._phi is None:
            self._equivalence_ratio()
        return self._phi

    @phi.setter
    def phi(self, phi_targ):
        self._change_phi(phi_targ)

    @property
    def fuels(self):
        ''' Find the fuel components. '''
        if self._fuels is None:
            self._components()
        return self._fuels

    @property
    def oxidizer(self):
        ''' Find the fuel components. '''
        if self._oxidizer is None:
            self._components()
        return self._oxidizer

    @classmethod
    def check_or_create(cls, mix, chemfile):
        '''
        If mix is a list, create a mixture object
        If mix is already a mixture object, check that the chemfile is correct
        '''
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
    '''
    Find the chemistry file in this working directory.
    If the chemistry file is in chemkin format, convert it to cantera
    Preference order is: chem.cti, chem.inp, *.cti, *.inp

    Returns the path to a cantera chemistry file in the workingdir
    '''
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
    "Convert chemkin chemistry file to cantera chemistry file."
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
