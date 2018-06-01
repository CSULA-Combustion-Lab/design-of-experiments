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


def fun1():
    pass

def fun2():
    pass



if __name__ == "__main__":
    # This is the code that will run when you call this file.
    fun1()
    fun2()