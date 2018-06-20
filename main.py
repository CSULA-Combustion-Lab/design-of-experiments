""" This is just a change to the discription to test 
pushing a change in the code.  
 """

import sys
import os
import cantera
cantera.suppress_thermo_warnings()

# import ignition tools from utilitites repository
dirname=os.path.dirname(os.path.abspath(os.path.curdir))
dirname=dirname+'/utilities'
sys.path.append(dirname)
import ignition

#create gas object
gas1 = cantera.Solution('gri30.xml')

#set state
gas1.TP = 800, 8*101325   

#set gas compsotion; X - mole fraction;  Y - mass fraction
gas1.X = 'CH4:0.5, O2:2, N2:7.5'

ignition.run(self)

def fun1():
    pass

def fun2():
    pass



if __name__ == "__main__":
    # This is the code that will run when you call this file.
    fun1()
    fun2()