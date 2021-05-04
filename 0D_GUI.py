#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 14:04:17 2018

@author: boss
"""
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib import style
import cantera
import pickle
import os

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog

LARGE_FONT= ("Verdana", 12)  #type of font and size
EXTRALG_FONT= ("Verdana", 18)
EXTRALG_FONT2= ("Verdana", 14)
style.use("ggplot") #ggplot, dark_background,

# Choose directory of stored data
# For some reason, filedialog.askdirectory makes the rest of the gui
# behave strangely. For now, I'm just hard-coding the directory
# load_dir = filedialog.askdirectory(initialdir='Outputs')
load_dir = r'C:\Users\jsantne\Documents\GitHub\Python-Code\Outputs\11_02_2020 21.25.16 Stanford'
print(load_dir)

# Load data
with open(os.path.join(load_dir, 'All_Ranks.pkl'), 'rb') as f:
    All_Ranks = pickle.load(f)

with open(os.path.join(load_dir, 'Max_Sens_Rxn.pkl'), 'rb') as f:
    score2_Max_sens_rxn = pickle.load(f)

with open(os.path.join(load_dir, 'Species_rxn.pkl'), 'rb') as f:
    species_rxn_file = pickle.load(f)

with open(os.path.join(load_dir, 'Case_Parameters.pkl'), 'rb') as f:
    Parameters = pickle.load(f)

with open(os.path.join(load_dir, 'Mole_fractions.pkl'), 'rb') as f:
    All_tTP_AMo = pickle.load(f)

with open(os.path.join(load_dir, 'All_time_sens.pkl'), 'rb') as f:
    All_time_Sens = pickle.load(f)


# Organize some of the loaded information
T = list(set([x['temperature'] for x in Parameters]))
P = list(set([x['pressure'] for x in Parameters]))
Phi = list(set([x['equivalence'] for x in Parameters]))
Fuel = list(set([x['fuel'] for x in Parameters]))
SS = species_rxn_file[0]

t_final = All_time_Sens[0][-1][0]
if t_final < 1.5e-3:
    t_scale = 1e6
    t_label = r'Time [$\mu$s]'
elif t_final < 1.5:
    t_scale = 1e3
    t_label = 'Time [ms]'
else:
    t_scale = 1
    t_label = 'Time [s]'

# THIS INFORMATION SHOULD COME FROM LOADED DATA SOMEHOW
mechanism = 'mech-FFCM1_modified.cti'
gas1 =cantera.Solution(mechanism)
S = gas1.species_names  # species_names
RR = gas1.reaction_equations()

root = tk.Tk()

#labels,Buttons,frames, optionmenus
# TO DO: use str.format() to improve T, P, phi, fuel drop-down menus.
# May have to use np.isclose() when searching for conditions
topFrame = tk.Frame(root)
topFrame.pack()
botFrame = tk.Frame(root)
botFrame.pack()

label1 = tk.Label(topFrame, text="Temperature [K]", font=LARGE_FONT)
label1.grid(row=1,column=0)
var1=tk.StringVar()
var1.set(T[0])
set1 = tk.OptionMenu(topFrame,var1, *T)
set1.configure(font=LARGE_FONT)
set1.grid(row=1,column=1)

label2 = tk.Label(topFrame, text="Pressure [kPa]", font=LARGE_FONT)
label2.grid(row=1,column=2)
var2=tk.StringVar()
var2.set(P[0])
set2 = tk.OptionMenu(topFrame,var2, *P)
set2.configure(font=LARGE_FONT)
set2.grid(row=1,column=3)

label3 = tk.Label(topFrame, text="Equivalence Ratio", font=LARGE_FONT)
label3.grid(row=2,column=0)
var3=tk.StringVar()
var3.set(Phi[0])
set3 = tk.OptionMenu(topFrame,var3, *Phi)
set3.configure(font=LARGE_FONT)
set3.grid(row=2,column=1)

label4 = tk.Label(topFrame, text="Fuel Mole Fraction", font=LARGE_FONT)
label4.grid(row=2,column=2)
var4=tk.StringVar()
var4.set(Fuel[0])
set4 = tk.OptionMenu(topFrame,var4, *Fuel)
set4.configure(font=LARGE_FONT)
set4.grid(row=2,column=3)

label5 = tk.Label(topFrame, text="Sensitivity of:", font=LARGE_FONT)
label5.grid(row=3,column=0)

var5=tk.StringVar()
var5.set(SS[0])
set5 = tk.OptionMenu(topFrame,var5, *SS)
set5.configure(font=LARGE_FONT)
set5.grid(row=3,column=1)

label6 = tk.Label(topFrame, text="Additional Species:", font=LARGE_FONT)
label6.grid(row=3, column=2)

var6=tk.StringVar()
set6 = tk.OptionMenu(topFrame,var6, *S)
set6.configure(font=LARGE_FONT)
set6.grid(row=3,column=3)

var7=tk.StringVar()
set7 = tk.OptionMenu(topFrame,var7, *S)
set7.configure(font=LARGE_FONT)
set7.grid(row=3,column=4)

var8=tk.StringVar()
set8 = tk.OptionMenu(topFrame,var8, *S)
set8.configure(font=LARGE_FONT)
set8.grid(row=4,column=3)

var9=tk.StringVar()
set9 = tk.OptionMenu(topFrame,var9, *S)
set9.configure(font=LARGE_FONT)
set9.grid(row=4,column=4)

label7 = tk.Label(topFrame, text="Time [s]", font=LARGE_FONT)
label7.grid(row=4,column=0)

entry3=tk.Entry(topFrame,width=5)
entry3.grid(row=4,column=1)

#functions
def do():
    Temp = float(var1.get())
    Press = float(var2.get())
    phi = float(var3.get())
    fuell = float(var4.get())
    specie1 = var5.get()
    specie2 = var6.get()
    specie3 = var7.get()
    specie4 = var8.get()
    specie5 = var9.get()
    species =[specie1,specie2,specie3,specie4,specie5]
    target = [phi, fuell, Temp, Press]
    for i in range(0,len(All_Ranks)):
        xx=All_Ranks[i][0]
        test = [xx['equivalence'], xx['fuel'], xx['temperature'], xx['pressure']]
        if np.all(np.isclose(test, target)):
            mixcount = i
            break
    try:
        mixcount
    except UnboundLocalError:
        # Should clear the figures in the future
        print('Condiiton not found:')
        print(str(target))
        return

    # print(mixcount)
    speciecount=[0,0,0,0,0]
    for i in range(0, len(S)):
        if specie1 == S[i]:
            speciecount[0] = i
        if specie2 == S[i]:
            speciecount[1] = i
        if specie3 == S[i]:
            speciecount[2] = i
        if specie4 == S[i]:
            speciecount[3] = i
        if specie5 == S[i]:
            speciecount[4] = i
    plotb(mixcount, speciecount,species)
    sensitive_rxn(mixcount,species,speciecount)

button1 = tk.Button(topFrame, text="Plot it!", command=do, font=LARGE_FONT)
button1.grid(row=4, column=2)

def sensitive_rxn(mixcount, species,speciecount):
    err = 100
    loc = 0
    try:
        time = float(entry3.get())
    except ValueError:
        return
    y=All_tTP_AMo[mixcount]
    for i in range(0,len(y)-1):
        y1=y[i] #iterates through time steps of chosen simulation
        y11=y1[0] #grabs first item in above list (should be time)
        if np.absolute(y11 - time) <= err:
            loc = i
            err = np.absolute(y11 - time)
    # print(loc) # prints location of least error

    specnum = 0
    specspecies=var5.get()
    for i in range (0, len(SS)):
        if specspecies == SS[i]:
            specnum = i
            # print(SS[i])
    # print(specspecies)
    # print(specnum)
    # print(loc+specnum*len(y))
    x=All_Ranks[mixcount]
    x1=x[loc*len(SS)+specnum+1]
    rxn_1=int(x1.index(1))
    rxn_2=int(x1.index(2))
    rxn_3=int(x1.index(3))
    rxn_4=int(x1.index(4))
    rxn_5=int(x1.index(5))
    rxn_=[rxn_1,rxn_2,rxn_3,rxn_4,rxn_5]
    rxn__ = [x+1 for x in rxn_] #use for labels of plot lines
    R=(rxn_,rxn__)
    # print(rxn_)
    # print(rxn__)
    vert_line.set_xdata(y[loc][0] * t_scale)
    plota(mixcount, species, speciecount, rxn_, rxn__, loc)


def plotb(mixcount, speciecount,species):
    time=[]
    molfracs = []
    y=All_tTP_AMo[mixcount]
#    y=All_tTP_SMo_AMo_SMa_AMa[mixcount]
    for i in range(0,len(y)-1):
        y1=y[i]
        time.append(y1[0] * t_scale)
        y22_=y1[3]  # Originally took abs value
        row = [y22_[int(x)] for x in speciecount]
        molfracs.append(row)
    molfracs = np.array(molfracs)

    maxtime=max(time)
    maxX = molfracs.max()

    if maxX < 1e-3:
        b.set_ylabel('Mole Fraction [ppm]')
        X_mult = 1e6
        maxX *= X_mult
    else:
        b.set_ylabel('Mole Fraction')
        X_mult = 1
    lines = [line1_1, line1_2, line1_3, line1_4, line1_5]
    for line, ind in zip(lines, range(5)):
        if species[ind] == '':
            continue
        line.set_xdata(time)
        line.set_ydata(X_mult * molfracs[:, ind])
        line.set_label(species[ind])
        line

    b.legend()
    b.axis((0-0.1*maxtime, maxtime+0.1*maxtime, 0-0.1*maxX, maxX+0.1*maxX))
    b.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4))
    canvas.draw()


def plota(mixcount, species, speciecount, rxn_, rxn__, loc):
    time=[]
    sens=[]
    sens1=[]
    sens2=[]
    sens3=[]
    sens4=[]
    y=All_time_Sens[mixcount]
    for i in range(0,len(y)-1):
        y1=y[i]
        y21=y1[0] #time
        y22=y1[1] #sensitivity
        y22_= y22[int(speciecount[0]), int(rxn_[0])]
        y22_1= y22[int(speciecount[0]),int(rxn_[1])]
        y22_2= y22[int(speciecount[0]),int(rxn_[2])]
        y22_3= y22[int(speciecount[0]),int(rxn_[3])]
        y22_4= y22[int(speciecount[0]),int(rxn_[4])]
        time.append(y21 * t_scale)
        sens.append(y22_)
        sens1.append(y22_1)
        sens2.append(y22_2)
        sens3.append(y22_3)
        sens4.append(y22_4)

    maxtime = max(time)
    all_sens = np.array([sens, sens1, sens2, sens3, sens4]).T
    minsens = all_sens[loc, :].min()
    maxsens = all_sens[loc, :].max()

    line.set_xdata( time)
    line.set_ydata( sens )
    line.set_label(RR[rxn_[0]])
    line

    line1.set_xdata( time)
    line1.set_ydata( sens1 )
    line1.set_label(RR[rxn_[1]])
    line1

    line2.set_xdata( time)
    line2.set_ydata( sens2 )
    line2.set_label(RR[rxn_[2]])
    line2

    line3.set_xdata( time)
    line3.set_ydata( sens3 )
    line3.set_label(RR[rxn_[3]])
    line3

    line4.set_xdata( time)
    line4.set_ydata( sens4 )
    line4.set_label(RR[rxn_[4]])
    line4

    a.set_ylabel('Sensitivity of ' + species[0])
    a.legend()
    y_delta = 0.2 * max(abs(minsens), abs(maxsens))
    a.axis([-0.1*maxtime, 1.1 * maxtime,
           minsens - y_delta, maxsens + y_delta])
    a.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4))
    canvas.draw()

#figures
f = Figure(figsize=(12,5), dpi=100)
#first plot
a = f.add_subplot(121)
a.set_xlabel(t_label)
a.set_ylabel('Sensitivity')
#second plot
b = f.add_subplot(122)
b.set_xlabel(t_label)

canvas = FigureCanvasTkAgg(f, botFrame)
canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand = True)
#plota
line, = a.plot([], [], "r-")
line1, = a.plot([], [], "b--")
line2, = a.plot([], [], "k-.")
line3, = a.plot([], [], "g-")
line4, = a.plot([], [], "m:")
a.axhline(y=0, color='k', ls='-')
vert_line = a.axvline(x=0, c='gray', ls='--')

#plotb
line1_1, = b.plot([], [], "r-")
line1_2, = b.plot([], [], "b--")
line1_3, = b.plot([], [], "k-.")
line1_4, = b.plot([], [], "g-")
line1_5, = b.plot([], [], "m:")


toolbar = NavigationToolbar2Tk(canvas, botFrame)
toolbar.update()
canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand = True)



root.mainloop()