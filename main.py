""" Description here
and here
and here
 """

import sys
import os
import cantera
cantera.suppress_thermo_warnings()

# import ignition tools from utilitites repository
dirname = os.path.dirname(os.path.normpath(os.path.dirname(__file__)))
sys.path.append(dirname)
from utilities import ignition
