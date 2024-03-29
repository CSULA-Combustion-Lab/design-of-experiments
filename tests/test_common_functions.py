# -*- coding: utf-8 -*-
"""
To use these tests, go to the module directory (design-of-experiments)
 in the anaconda prompt and type "pytest"

Make sure you have downloaded module pytet before attempting to use it.

Instructions:
    cd "[The Direcotry for design-of-experiments on your computer]"
    pytest

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
    time.sleep(0.5)
    return np.sum(array)


def test_parallelize():

    n = 10
    params = [[np.random.randint(0, 50, 10)] for i in range(n)]

    start = time.time()
    output_s = []
    for param in params:
        output_s.append(function(param))
    series = time.time() - start
    print(start, series)

    start = time.time()
    output_p = cf.parallelize(params, None, function)
    parallel = time.time() - start
    print(start, parallel)

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
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'phi_fuel'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['H2'], fuel[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    assert all([x >= 0 for k, x in mixture.items()])

def test_case_maker_phi_oxi():
    phi = [0.5 + rand(), 2, 1]
    oxidizer = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_oxi', phi, oxidizer)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'phi_oxi'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['O2'], oxidizer[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    assert all([x >= 0 for k, x in mixture.items()])

def test_case_maker_phi_dil():
    phi = [0.5 + rand(), 2, 1]
    diluent = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_dil', phi, diluent)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'phi_dil'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['N2'], diluent[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(phi[0], mixture['H2'] / mixture['O2'] * 0.5)
    assert all([x >= 0 for k, x in mixture.items()])

def test_case_maker_phi_oxidil():
    phi = [0.3 + rand(), 2, 1]
    ox_in_oxidizer = [0.2 + rand() / 2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_oxi/dil', phi, ox_in_oxidizer)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'phi_oxi/dil'],
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

def test_case_maker_phi_fuedil():
    phi = [0.3 + rand(), 2, 1]
    f_in_fuel = [0.2 + rand() / 2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_fuel/dil', phi, f_in_fuel)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'phi_fuel/dil'],
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

def test_case_maker_oxi_fuel():
    oxi = [rand()/2, 1, 1]
    fuel = [rand() / 2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('oxi_fuel', oxi, fuel)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'oxi_fuel'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]
    print(mixture)

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    assert all([x >= 0 for k, x in mixture.items()])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(oxi[0], mixture['O2'])
    nptest.assert_allclose(fuel[0], mixture['H2'])

def test_case_maker_fuel_dil():
    fuel = [rand() / 2, 1, 1]
    dil = [rand()/2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('fuel_dil', fuel, dil)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'fuel_dil'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]
    print(mixture)

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    assert all([x >= 0 for k, x in mixture.items()])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(dil[0], mixture['N2'])
    nptest.assert_allclose(fuel[0], mixture['H2'])

def test_case_maker_oxi_dil():
    oxi = [rand() / 2, 1, 1]
    dil = [rand()/2, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('oxi_dil', oxi, dil)
    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2':1}, {'N2':1}, {'O2':1}, 'oxi_dil'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]
    print(mixture)

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    assert all([x >= 0 for k, x in mixture.items()])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'])
    nptest.assert_allclose(dil[0], mixture['N2'])
    nptest.assert_allclose(oxi[0], mixture['O2'])

def test_multi_f_o():
    phi = [0.3 + rand(), 2, 1]
    fuel = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_fuel', phi, fuel)

    conditions = {'Parameters': [P, T, mix_params, 'log'],
                  'Mixture': [{'H2': .3, 'CO': 0.7}, {'N2':0.5, 'HE': 0.5}, {'O2': 0.95, 'AR': 0.05}, 'phi_fuel'],
                  'Files': [cf.model_folder('grimech30.cti'), None]}
    paramlist = cf.case_maker(conditions)
    case = paramlist[0]
    mixture = case[2]

    nptest.assert_allclose(case[:2], [P[0], T[0]])
    nptest.assert_allclose(mixture['H2'] + mixture['CO'], fuel[0])
    nptest.assert_allclose(1, mixture['H2'] + mixture['O2'] + mixture['N2'] + mixture['AR'] + mixture['CO'] + mixture['HE'])
    nptest.assert_allclose(phi[0], (mixture['H2'] + mixture['CO']) / mixture['O2'] * 0.5)
    nptest.assert_allclose(0.3 / 0.7, mixture['H2'] / mixture['CO'])
    nptest.assert_allclose(0.95 / 0.05, mixture['O2'] / mixture['AR'])
    assert all([x >= 0 for k, x in mixture.items()])


def test_mixture_percentage():
    phi = [0.3 + rand(), 2, 1]
    fuel = [0.2 + rand() / 5, 1, 1]
    P = [rand(), 100, 1]
    T = [500 * rand(), 1000, 1]
    mix_params = ('phi_fuel', phi, fuel)

    multif = {'H2': .3, 'CO': 0.7}
    singlf = {'H2': 1}

    list_conditions = {'Parameters': [P, T, mix_params, 'log'],
                       'Mixture': [multif, {'N2':1}, {'O2':1}, 'phi_fuel'],
                       'Files': [cf.model_folder('grimech30.cti'), None]}
    stri_conditions = {'Parameters': [P, T, mix_params, 'log'],
                       'Mixture': [singlf, {'N2':1}, {'O2':1}, 'phi_fuel'],
                       'Files': [cf.model_folder('grimech30.cti'), None]}
    list_paramlist = cf.case_maker(list_conditions)
    list_case = list_paramlist[0]
    list_mixture = list_case[2]
    stri_paramlist = cf.case_maker(stri_conditions)
    stri_case = stri_paramlist[0]
    stri_mixture = stri_case[2]

    lst_comp = cf.mixture_percentage(multif, list_mixture)
    str_comp = cf.mixture_percentage(singlf, stri_mixture)

    nptest.assert_allclose(list_case[:2], [P[0], T[0]])
    nptest.assert_allclose(stri_case[:2], [P[0], T[0]])
    nptest.assert_allclose(lst_comp, fuel[0])
    nptest.assert_allclose(str_comp, fuel[0])
    assert all([x >= 0 for k, x in list_mixture.items()])
    assert all([x >= 0 for k, x in stri_mixture.items()])


def test_mixture_normalize():
    # Using a list
    lst = ['comp1', rand(), 'comp2', rand(), 'comp3', rand()]
    total = sum(lst[1::2])
    mix_dict = cf.normalize_mixture(lst)
    nptest.assert_allclose([x/total for x in lst[1::2]], [mix_dict[k] for k in lst[::2]])
    nptest.assert_allclose(1, sum(mix_dict.values()))

    # Using a dict
    dct = {'comp1': rand(), 'comp2': rand(), 'comp3': rand()}
    total = sum(dct.values())
    mix_dict = cf.normalize_mixture(dct)
    nptest.assert_allclose([x/total for x in dct.values()], list(mix_dict.values()))
    nptest.assert_allclose(1, sum(mix_dict.values()))

    # using a string
    assert cf.normalize_mixture('test') == {'test': 1.0}
