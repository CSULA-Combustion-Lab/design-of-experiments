# -*- coding: utf-8 -*-
"""
To use these tests, go to the module directory (design-of-experiments)
 in the anaconda prompt and type "pytest"

"""

import os
import sys
import numpy as np
import numpy.testing as nptest
import time
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
import common_functions as cf


def function(array, dum=None):
    time.sleep(1e-3)
    return np.sum(array)


def test_parallelize():

    n = 1000
    params = [[np.random.randint(0, 50, 10)] for i in range(n)]

    start = time.time()
    output_s = []
    for param in params:
        output_s.append(function(param))
    series = time.time() - start

    start = time.time()
    output_p = cf.parallelize(params, None, function)
    parallel = time.time() - start

    assert parallel / series < 0.8
    assert output_p == output_s


rand = np.random.rand  # Define random number function for convenience

def test_case_maker_phi_fuel():
    phi = [0.5 + rand(), 2, 1]
    fuel = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_fuel', phi, fuel)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': ['H2', 'N2', 'O2', 'phi_fuel'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['H2'], fuel[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    assert all([x >= 0 for k, x in mixture.items()])

def test_case_maker_oxi_dil():
    phi = [0.3 + rand(), 2, 1]
    ox_in_oxidizer = [0.2 + rand() / 2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('Oxi_Dil', phi, ox_in_oxidizer)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': ['H2', 'N2', 'O2', 'Oxi_Dil'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    assert all([x >= 0 for k, x in mixture.items()])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    nptest.assert_allclose(mixture['O2']/(mixture['N2'] + mixture['O2']),
                           ox_in_oxidizer[0])

def test_case_maker_fue_dil():
    phi = [0.3 + rand(), 2, 1]
    f_in_fuel = [0.2 + rand() / 2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('Fue_Dil', phi, f_in_fuel)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': ['H2', 'N2', 'O2', 'Oxi_Dil'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    assert all([x >= 0 for k, x in mixture.items()])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    nptest.assert_allclose(mixture['H2']/(mixture['N2'] + mixture['H2']),
                           f_in_fuel[0])

def test_multi_f_o():
    phi = [0.3 + rand(), 2, 1]
    fuel = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_fuel', phi, fuel)

    multif = ['H2', .3, 'CO', 0.7]
    multio = ['O2', 0.95, 'AR', 0.05]

    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [multif, 'N2', multio, 'phi_fuel'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['H2'] + mixture['CO'], fuel[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'] + mixture['AR'] + mixture['CO'])
    nptest.assert_allclose(phi[0], (mixture['H2'] + mixture['CO']) / mixture['O2'] * 0.5)
    nptest.assert_allclose(0.3 / 0.7, mixture['H2'] / mixture['CO'])
    nptest.assert_allclose(0.95 / 0.05, mixture['O2'] / mixture['AR'])
    assert all([x >= 0 for k, x in mixture.items()])

