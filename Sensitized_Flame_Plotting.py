# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:20:01 2020

@author: Kodo Bear
"""

import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpl_patches
#from tkinter import filedialog


if __name__ == "__main__":
    Folder_name = input('Please type name of folder.'
                        '\n If blank, use last folder:\n')
    if Folder_name == '':
        with open('last run.pkl', 'rb') as f:
            Folder_name = pickle.load(f)
        print('Loading ' + Folder_name)
    Load_folder = '\\'+Folder_name
    with open('last run.pkl', 'wb') as f:
        pickle.dump(Folder_name, f)

    Load_path   = 'Flame_Sensitivity_Results'+Load_folder
    Plot_path   = Load_path+'\\Flame_Sensitivity_Plots'

    file = open(Load_path+'\\Case Description.txt','r')
    Description = file.read()
    file.close()
    print(Description)

    with open(os.path.join(Load_path, 'Flame Information.pkl'), 'rb') as f:
        flame_info = pickle.load(f)