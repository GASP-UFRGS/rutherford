import numpy as np
import sys
from card_reader import read_card, _raise_missing_card_error
from scipy.constants import epsilon_0, pi, e, alpha

try:
    card_name = sys.argv[1]
except IndexError:
    _raise_missing_card_error()

# Constants

kinEn, zTarget, zProj, angUnit, angStart, angEnd = read_card(card_name)

kinEn = kinEn*e # Energy of the particles in J.
kconst = 1/(4*pi*epsilon_0)
r0 = (kconst*zProj*zTarget*e**2/kinEn)*1e15 # Minimum distance between incident particles and target in fm.

# Functions

def impact_parameter(theta, r0, angle_unit):
    """
    Returns impact parameter when given the scattering angle.
    """

    if angle_unit == 'radians':
        bparam = r0/(2*np.tan(theta/2))
    elif angle_unit == 'degrees':
        bparam = r0/(2*np.tan(np.radians(theta)/2))
    return bparam

def scattering_angle(bparam, r0, angle_unit):
    """
    Returns scattering angle when given the impact parameter.
    """

    if angle_unit == 'radians':
        return 2*np.arctan(r0/(2*bparam))
    elif angle_unit == 'degrees':
        return np.degrees(2*np.arctan(r0/(2*bparam)))

def scattering_differential(theta, r0, angle_unit):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if angle_unit == 'radians':
        dsig = np.pi*r0**2*np.cos(theta/2)
        dsigdtheta = dsig/(4*np.sin(theta/2)**3)
    elif angle_unit == 'degrees':
        theta = np.radians(theta)
        dsig = np.pi*r0**2*np.cos(theta/2)
        dsigdtheta = dsig/(4*np.sin(theta/2)**3)
    return dsigdtheta

def scattering_differential_Mott(theta, angle_unit):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if angle_unit == 'radians':
        dsig = e**4*np.cos(theta/2)**3
        dsigdtheta = dsig/(16*kinEn**2*np.pi*np.sin(theta/2)**3)
    elif angle_unit == 'degrees': #CGS
        theta = np.radians(theta)
        dsig = np.pi*alpha**2*np.cos(theta/2)**3
        dsigdtheta = dsig/((kinEn*10e7)**2*np.sin(theta/2)**3)
    return dsigdtheta

# Calculations

theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.
b_out = impact_parameter(theta_in,r0,angUnit) # Impact parameter calculated.

dsig_dtheta = scattering_differential(theta_in,r0,angUnit) #Differential scattering cross section
dsig_dtheta_Mott = scattering_differential_Mott(theta_in,angUnit) #Differential scattering cross section

# Write to file

file_path = "output.dat"

data = np.column_stack((theta_in, b_out, dsig_dtheta, dsig_dtheta_Mott))
np.savetxt(file_path, data, delimiter=",")

# Write angle unit
with open(file_path, "a") as file:
     file.write("# " + angUnit + "\n")

print("Data has been written to", file_path)
