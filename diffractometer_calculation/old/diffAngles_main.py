# -*- coding: utf-8 -*-
"""
Diffraction angle calculator

4C: based on Busing(1967), Mochrie(1988)
"""


import numpy as np

import diffAngles as diff

""" --------------- PREAMBLE ---------- """
a = 5.387 # Angstrom
b = 5.383
c = 7.609
alpha = 90 # deg
beta = 90
gam = 90

alp = 5

E = 8.3 # keV

hkl = np.array([1,0,3])
N = np.array([1,1,0])

""" ----------------------------------- """

wavelength = 12.3984/E
K = 2*np.pi/wavelength

a = np.array([a,b,c]) # lattice const
aa = np.array([alpha,beta,gam]) # lattice angles
a0,a1,a2 = diff.vectorFromLengthAndAngles(a,aa)
a,aa,b,ba = diff.lengthAndAngles(a0,a1,a2)
U,B = diff.UBmat(a, aa, N)

ttheta, theta, chi, phi, omega = diff.hklToAngles_4C(hkl, N, E, U, B, 1, alp)

print('ttheta = '+str(np.rad2deg(ttheta)))
print('theta = '+str(np.rad2deg(theta)))
print('chi = '+str(np.rad2deg(chi)))
print('phi = '+str(np.rad2deg(phi)))
print('omega = '+str(np.rad2deg(omega)))



























