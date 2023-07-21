import numpy as np
import sys
from card_reader import read_card, _raise_missing_card_error 
from scipy.constants import epsilon_0, pi, e, alpha, hbar, c

try:
    card_name = sys.argv[1]
except IndexError:
    _raise_missing_card_error()


# Constants

kinEn, zTarget, zProj, angUnit, angStart, angEnd, mott, mass = read_card(card_name)

kinEn = kinEn*e # Energy of the particles in J.
kconst = 1/(4*pi*epsilon_0)
D = (kconst*zProj*zTarget*e**2/kinEn)*1e15 # Minimum distance between incident particles and target in fm.
D2 = (zProj*zTarget*alpha*hbar*c/(2*kinEn))*1e15

# Either cos or omega. This exists to test conversion between equations
var = 'omega'


# Functions

def impact_parameter(theta, D, angle_unit):
    """
    Returns impact parameter when given the scattering angle.
    """ 
    if angle_unit == 'degrees':
        theta = np.radians(theta)

    return D/(2*np.tan(theta/2))



def scattering_angle(bparam, D, angle_unit):
    """
    Returns scattering angle when given the impact parameter.
    """
    if angle_unit == 'degrees':
        theta = np.radians(theta)

    return 2*np.arctan(D/(2*bparam))



def scattering_differential_Ruth(theta, D, angle_unit):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if angle_unit == 'degrees':
        theta = np.radians(theta)

    if var == 'cos':
        difCrossSec = (2*pi*D**2/(1-np.cos(theta))**2)#(np.sin(theta)) #Rolf
    
    if var == 'omega':        
        difCrossSec = D**2/(16*np.sin(theta/2)**4)                                  

    # Only dOmega and dcos are trustworthy.
    return difCrossSec


def scattering_differential_Mott(theta, D, angle_unit):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if angle_unit == 'degrees':
        theta = np.radians(theta)


    if var == 'cos':
        # Rohlf
        difCrossSec = (2*pi*D**2/(1-np.cos(theta))**2)
        difCrossSec_Mott = difCrossSec*((1+np.cos(theta))/(2*(1+(((1-np.cos(theta))*kinEn)/(mass*c**2)))))

    if var == 'omega':
        difCrossSec = D2**2/(4*np.sin(theta/2)**4) 
        difCrossSec_Mott = difCrossSec*np.cos(theta/2)**2   

    return difCrossSec_Mott 


# Calculations

theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.
b_out = impact_parameter(theta_in, D, angUnit) # Impact parameter calculated.

dsig_dtheta = scattering_differential_Ruth(theta_in, D, angUnit) #Differential scattering cross section

if mott == 'y':
    dsig_dtheta_Mott = scattering_differential_Mott(theta_in, D, angUnit) # Mott correction cross section


if var == 'cos':
    theta_in = np.cos(np.radians(theta_in))


if angUnit == 'radians':
    theta_in = np.degrees(theta_in)

# Write to file

file_path = "output.dat"

data = np.column_stack((theta_in, b_out, dsig_dtheta, dsig_dtheta_Mott))

np.savetxt(file_path, data, delimiter=",")


print("Data has been written to", file_path)
