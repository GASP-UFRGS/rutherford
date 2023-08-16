import numpy as np
import sys
import periodictable
from card_reader import read_card, _raise_missing_card_error 
from scipy.constants import epsilon_0, pi, e, c, alpha, hbar, electron_mass

try:
    card_name = sys.argv[1]
except IndexError:
    _raise_missing_card_error()

# Constants

parameters = read_card(card_name)
kinEn = parameters.get('kinEn')
zTarget = parameters.get('zTarget') 
zProj = parameters.get('zProj') 
angUnit = parameters.get('angUnit') 
angStart = parameters.get('angStart') 
angEnd = parameters.get('angEnd')
mott = parameters.get('mott') 
recoil = parameters.get('recoil')
impactParameter = parameters.get('impactParameter') 
cross_section_variable = parameters.get('difCrossSec')

# Outputs Nuclear mass in atomic mass units (u).
element = periodictable.elements[zTarget]
massTarget = element.mass
massTarget = massTarget*1.6605402E-27 # Converts mass to kg


# If projectile is electron
if zProj != 0:
    massProj = periodictable.elements[zProj].mass * 1.6605402E-27
else:
    zProj = 1
    massProj = electron_mass


kinEn = kinEn*e # Converts energy of incoming particles to Joules.
kconst = 1/(4*pi*epsilon_0)
fm = 1e15 # conversion factor to femtometer.
D = (kconst*zProj*zTarget*e**2/kinEn) * fm # Minimum distance between incident particles and target in fm.


# Functions

def impact_parameter(theta, D):
    """
    Returns impact parameter when given the scattering angle.
    """ 
    return D/(2*np.tan(theta/2))



def scattering_angle(bparam, D):
    """
    Returns scattering angle when given the impact parameter.
    """
    return 2*np.arctan(D/(2*bparam))



def scattering_differential_Ruth(theta, D, cross_section_variable):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if cross_section_variable == 'cos':
        difCrossSec_Ruth = (2*pi*D**2/(1-np.cos(theta))**2)

    if cross_section_variable == 'theta': 
        difCrossSec_Ruth = (D**2*pi*np.cos(theta/2)/(4*np.sin(theta/2)**3))
    
    if cross_section_variable == 'omega':  
        difCrossSec_Ruth = D**2/(16*np.sin(theta/2)**4)     

    return difCrossSec_Ruth



def scattering_differential_Mott(theta, difCrossSec_Ruth, kinEn, massTarget):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    difCrossSec_Mott = difCrossSec_Ruth * np.cos(theta/2)**2 

    return difCrossSec_Mott 

def scattering_differential_Recoil(theta, difCrossSec_Mott, kinEn, massTarget):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    difCrossSec_Recoil = difCrossSec_Mott * (1/(1+(((1-np.cos(theta))*kinEn)/(massTarget*c**2))))

    return difCrossSec_Recoil


# Calculations

theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.

# Converts theta to radians if necessary
if angUnit == 'degrees':
    theta_in = np.radians(theta_in)

# Calculates Impact parameter
if impactParameter == 'true':
    b_out = impact_parameter(theta_in, D) 

# Calculates diferential coss section
if cross_section_variable in ['cos', 'theta', 'omega']:

    difCrossSec_Ruth = scattering_differential_Ruth(theta_in, D, cross_section_variable) #Differential scattering cross section.

    if mott == 'true':
        difCrossSec_Mott = scattering_differential_Mott(theta_in, difCrossSec_Ruth, kinEn, massTarget) # Mott correction cross section.

    if recoil == "true":   
         difCrossSec_Recoil = scattering_differential_Recoil(theta_in, difCrossSec_Mott, kinEn, massTarget) # Target recoil correction cross section.


# Write to file

file_path = "output.dat"

header = 'theta'
data = np.degrees(theta_in)

if impactParameter == 'true':
    header += ',b_out'
    data = np.column_stack((data, b_out))

if cross_section_variable in ['cos', 'theta', 'omega']:
    header += ',difCrossSec_Ruth'
    data = np.column_stack((data, difCrossSec_Ruth))

    if mott == "true":
        header += ',difCrossSec_Mott'
        data = np.column_stack((data, difCrossSec_Mott))

    if recoil == "true":
        header += ',difCrossSec_Recoil'
        data = np.column_stack((data, difCrossSec_Recoil))

np.savetxt(file_path, data, delimiter=",", header = header, comments="")


print("Data has been written to", file_path)