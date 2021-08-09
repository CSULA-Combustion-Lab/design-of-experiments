"""GUI to calculate gas properties based on cantera.

Designed with similar functionality to GasEq

Created by Jeffrey Santner
"""

import os
import tkinter
from tkinter import ttk
from tkinter import filedialog
from tkinter import scrolledtext
import numpy as np
import yaml
import pickle
import matplotlib.pyplot as plt
import cantera
import common_functions as cf
import Simulation_Initializer
import Sensitized_Flame_Experiment


fmt = '{:.3g}'.format  # Number format used throughout


class GUI(ttk.Frame):
    """Main window."""


    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent
        self.init_gui()

    def init_gui(self):
        """Build the GUI."""
        self.root.title('Flame Sensitivity Plotter')
        self.grid(column=0, row=0, sticky='nsew')

        # Choose folder
        self.folder_button = ttk.Button(
            self, text='Choose condition',
            command=self.pick_folder)
        self.folder_button.grid(column=0, row=0)
        self.workingdir_var = tkinter.StringVar()
        b = ttk.Entry(self, textvariable=self.workingdir_var)
        b.grid(row=0, column=1)
        # Initialize Variables
        # Set things to None when needed

        # Other buttons
        plot_button = ttk.Button(self, text='Plot', command=self.plot)
        plot_button.grid(row=0, column=2)
        ttk.Label(self, text='Min. # of reactions:').grid(row=0, column=3)
        self.n = tkinter.IntVar()
        n_entry = ttk.Entry(self, textvariable=self.n)
        self.n.set(7)
        n_entry.grid(row=0, column=4)

        # Build empty grid for conditions
        self.max_cases = 7
        height = self.max_cases
        width = 5

        # Initialize grid headings
        self.headings = []
        for i in range(0, width):
            x = tkinter.StringVar()
            b = ttk.Entry(self, textvariable=x)
            b.grid(row=2, column=i)
            self.headings.append(x)
        [x.set(y) for x, y in zip(self.headings, ('', 'T [K]', 'P [atm]', 'var1', 'var2'))]

        cells = {}
        cell_data = {}
        for i in range(height):
            for j in range(width):
                # b = ttk.Entry(self, textvariable=x)
                if j == 0:
                    x = tkinter.StringVar()
                    b = ttk.Entry(self, textvariable=x)
                else:
                    x = tkinter.DoubleVar()
                    b = ttk.OptionMenu(self, x, ())
                    b.configure(state='disabled')
                b.grid(row=i+3, column=j)
                cells[(i, j)] = b
                cell_data[(i, j)] = x
        self.cells = cells
        self.cell_data = cell_data

        self.info = scrolledtext.ScrolledText(self, width= 50, height=12)
        self.info.grid(column=7, row=0, rowspan=10)


    def pick_folder(self):
        """Chose the folder. """
        # TEMP
        # self.workingdir = r'C:\Users\jsantne\Documents\GitHub\design-of-experiments\Flame_Sensitivity_Results\Methanol Tighter range'
        self.workingdir = filedialog.askdirectory(
                title='Select folder',
                initialdir=os.path.join(os.getcwd(), 'Flame_Sensitivity_Results'))
        self.workingdir_var.set(os.path.split(self.workingdir)[1])
        self.set_dropdown_menus_and_info()


    def set_dropdown_menus_and_info(self):
        """Setup the dropdown menus to choose mixture and thermodynamic variables. """
        # Open and parse the input file
        inputs = list(
            yaml.safe_load_all(open(
                self.workingdir + '/input.yaml', 'r')))
        sim_inputs = inputs[0]

        mix_options = sim_inputs['Mixture_options']
        mix_params = Simulation_Initializer.package(mix_options)

        # Unpackage other items
        Pressure = mix_options['Pressure']
        Temperature = mix_options['Temperature']
        fuel = mix_options['Fuel']
        oxi = mix_options['Oxidizer']
        dil = mix_options['Diluent']
        arrtype = mix_options['Array_type']
        Mechanism = sim_inputs['Mechanism']
        Mixture_type = mix_options['Mixture_type']
        problem_type = sim_inputs['Simulation_Type']

        # Column headings
        var1, var2 = Mixture_type.split('_')
        self.headings[3].set(var1)
        self.headings[4].set(var2)

        # Other information
        # TODO: create table with mechanism, array type, fuel, oxidizer, diluent
        text = cf.parameters_string(problem_type, Pressure, Temperature,
                                    mix_params, Mechanism, fuel, oxi, dil)
        self.info.insert(tkinter.INSERT, text)

        # Dropdowns
        def logspace(start, stop, number):
            """Custom, simplified logspace function."""
            return np.logspace(np.log10(start), np.log10(stop), number)

        if arrtype == 'log':
            function = logspace
        elif arrtype == 'lin':
            function = np.linspace

        # Temperatures
        for row in range(self.max_cases):
            for col, param in zip(range(1, 5), (Temperature, Pressure, mix_params[1], mix_params[2])):
                dropdown = self.cells[(row, col)]
                var = self.cell_data[(row, col)]
                dropdown.configure(state='normal')  # Enable drop down
                _reset_option_menu(dropdown, var, function(*param))


        # Save things for use later
        self.mixture_type = Mixture_type
        self.chemfile = cf.model_folder(Mechanism)
        # condi = Sensitized_Flame_Experiment.initialization(
        #     self.chemfile, arrtype, Pressure, Temperature, fuel, oxi, dil,
        #     mix_params, 0, False, 0)
        # self.paralist = cf.case_maker(condi)

    def plot(self):
        with open(self.workingdir + '/Flame Information.pkl', 'rb') as f:
            Flame_info = pickle.load(f)
        labels = []
        sens = []
        for row in range(self.max_cases):
            case, label = find_case(self.cell_data, row, self.headings,
                             self.mixture_type, self.chemfile, Flame_info)
            if case is not None and case['Flame'][0] is not None:
                labels.append(label)
                sens.append(case['Flame'][0])
                self.cell_data[(row, 0)].set(u'\u2713')  # Check mark
            else:
                self.cell_data[(row, 0)].set("Impossible or didn't converge.")
        if len(sens) > 0:
            sens_plot(sens, labels, self.n.get())


def sens_plot(sens, labels, n=7):
    # Sort them
    for condition in sens:
        condition.sort(key=lambda x: abs(x[1]), reverse=True)
    reactions = []
    rank = 0
    while len(reactions) < n:
        # This while loop will find at least the "n" most sensitive reactions
        # among all conditions.
        for condition in sens:
            rxn = condition[rank][2]
            if rxn not in reactions:
                reactions.append(rxn)
        rank += 1
    sens_plot = []
    for rxn in reactions:
        # This loop will populate the sensitivities to plot
        intermediate = []
        for condition in sens:
            intermediate.append([s[1] for s in condition if s[2] == rxn][0])
        sens_plot.append(intermediate)
    sens_plot = np.array(sens_plot)

    n = len(reactions)  # Adjust in case extra reactions were found
    fig, ax = plt.subplots()
    barheight = 1 / (len(labels) + 1)
    ax.grid(axis='x', which='major', ls='--')
    # ax.grid(axis='y', which='minor', c='k', ls='.')

    for i, label in enumerate(labels):
        ylocs = np.arange(n) + i * barheight
        ax.barh(ylocs, sens_plot[:, i], height=barheight, label=label,
                align='center')

    for y in np.arange(n):
        ax.axhline(y - barheight, c='k', ls='--')

    ax.set_yticks(np.arange(n) + barheight * (len(labels) - 1) / 2)
    ax.set_yticklabels(reactions, rotation=30, fontsize=15)
    ax.set_yticks(np.arange(n), minor=True)
#    ax.tick_params(axis='y', which='minor', bottom='off')
    ax.invert_yaxis()
    ax.axvline(c='k', ls='-')  # Why doesn't this show up?
    ax.set_xlabel('Normalized Sensitivity')
    ax.legend()
    # ax.set_ylim([max(ylocs)+0.5, min(ylocs)-0.5])
    # fig.tight_layout()
    fig.show()





def find_case(cell_data, row, headings, mix_type, chemfile, flame_info):
    """ Find the target case in flame_info. """
    T = cell_data[(row, 1)].get()
    P = cell_data[(row, 2)].get()
    var1 = cell_data[(row, 3)].get()
    var2 = cell_data[(row, 4)].get()

    if any([x == 0.0 for x in [T, P, var1, var2]]):
        return None, None

    # use cf.case_maker to create the target mixture
    flame_cond = flame_info[0]['Conditions']
    cond = {'Parameters': [[P, 1e10, 1], [T, 1e10, 1],
                           [mix_type, [var1, 1e10, 1], [var2, 1e10, 1]], 'lin'],
            'Mixture': [cf.normalize_mixture(flame_cond[6]),
                        cf.normalize_mixture(flame_cond[8]),
                        cf.normalize_mixture(flame_cond[7])],
            'Files': [chemfile, None]}
    try:
        target_mix = cf.case_maker(cond)[0]
    except IndexError:  # Impossible Mixture
        print('Mixture is impossible with {}, {}, {}, {}'.format(T, P, var1, var2))
        return None, None
    for case in flame_info:
        cond = case['Conditions']
        mix = cond[5]
        if np.allclose(target_mix[:2], cond[1::-1]):
            species = target_mix[2].keys()
            array_1 = np.array([target_mix[2][key] for key in species])
            array_2 = np.array([mix[key] for key in species])

            Su = case['Flame'][1] * 100
            if type(Su) is str:
                Su = 0.0

            label = series_label(T, P, headings[3].get(), var1,
                                 headings[4].get(), var2, Su)
            if np.allclose(array_1, array_2, rtol=3e-3):
                return case, label

    print('This case was not found')
    print(target_mix)
    return None, None

def series_label(T, P, h1, v1, h2, v2, Su):
    headings_dict = {'phi': r'$\phi$', 'oxi/dil': 'ox/(ox+dil)',
                     'fuel/dil': 'fuel/(fuel+dil)'}
    try:
        heading1 = headings_dict[h1]
    except KeyError:
        heading1 = h1

    try:
        heading2 = headings_dict[h2]
    except KeyError:
        heading2 = h2

    label = r'T={} K, P={} atm, {}={}, {}={}, $S_u$={:.2g} cm/s'.format(
                T, P, heading1, v1, heading2, v2, Su)
    return label

def _reset_option_menu(options_menu, om_variable, options):
        menu = options_menu["menu"]
        menu.delete(0, "end")
        for string in options:
            string = fmt(string)
            menu.add_command(label=string,
                             command=lambda value=string:
                                  om_variable.set(value))
        if len(options) == 1:
            om_variable.set(fmt(options[0]))


if __name__ == "__main__":
    root = tkinter.Tk()
    GUI(root)
    root.mainloop()
