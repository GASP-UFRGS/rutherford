#!/usr/bin/env python3

# import physical constants module
from scipy.constants import pi as pi
from scipy.constants import epsilon_0 as epsilon0

kund = 4*pi*epsilon0
k = 1/kund

print('We have kund = {0} and k = {1}'.format(kund,k))
