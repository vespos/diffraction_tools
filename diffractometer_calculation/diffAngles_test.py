# -*- coding: utf-8 -*-
"""
Created on Fri May 31 16:35:21 2019

@author: espov
"""

import numpy as np

import diffAngles as diff
import main_4C

""" --------------- PREAMBLE ---------- """
a = 5.387 # Angstrom
b = 5.383
c = 7.609
alpha = 90 # deg
beta = 90
gamma = 90

alp = 5

E = 8.3 # keV

hkl = np.array([1,0,3])
N = np.array([1,1,0])

""" ----------------------------------- """

a = np.array([a,b,c])
aa = np.array([alpha,beta,gamma])
ttheta = main_4C.main_4C(hkl, a, aa, N, E, alp)