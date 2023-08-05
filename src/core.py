import numpy as np
import sys
from card_reader import read_card, _raise_missing_card_error 
from scipy.constants import epsilon_0, pi, e, c

try:
    card_name = sys.argv[1]
except IndexError:
    _raise_missing_card_error()


# Constants

kinEn = parameters.get('kinEn')
zTarget = parameters.get('zTarget') 
zProj = parameters.get('zProj') 
angUnit = parameters.get('angUnit') 
angStart = parameters.get('angStart') 
angEnd = parameters.get('angEnd') 
mott = parameters.get('mott') 
recoil = parameters.get('recoil') 
massTarget = parameters.get('massTarget') 

kinEn = kinEn*e # Converts energy of incoming particles to Joules.
kconst = 1/(4*pi*epsilon_0)
fm = 1e15 # conversion factor to femtometer.
D = (kconst*zProj*zTarget*e**2/kinEn) * fm # Minimum distance between incident particles and target in fm.


# Either cos, theta or omega. 
var = parameters.get('var')


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



def scattering_differential_Ruth(theta, D):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if var == 'cos':
        difCrossSec_Ruth = (2*pi*D**2/(1-np.cos(theta))**2)

    if var == 'theta': 
        difCrossSec_Ruth = (D**2*pi*np.cos(theta/2)/(4*np.sin(theta/2)**3))
    
    if var == 'omega':        
        difCrossSec_Ruth = D**2/(16*np.sin(theta/2)**4)                                  

    return difCrossSec_Ruth



def scattering_differential_Mott(theta, difCrossSec_Ruth, D):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if var == 'cos':
        difCrossSec_Mott = difCrossSec_Ruth * ((1+np.cos(theta))/(2*(1+(((1-np.cos(theta))*kinEn)/(massTarget*c**2)))))

    if var == 'theta':
        difCrossSec_Mott = difCrossSec_Ruth * np.cos(theta/2)**2

    if var == 'omega':
        difCrossSec_Mott = difCrossSec_Ruth * np.cos(theta/2)**2   

    return difCrossSec_Mott 


# Calculations

theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.

# Convert theta to radians if necessary
if angUnit == 'degrees':
    theta_in = np.radians(theta_in)

b_out = impact_parameter(theta_in, D) # Impact parameter calculated.
difCrossSec_Ruth = scattering_differential_Ruth(theta_in, D) #Differential scattering cross section.


if mott == 'true':
    difCrossSec_Mott = scattering_differential_Mott(theta_in, difCrossSec_Ruth, D) # Mott correction cross section.


if var == 'cos':
    theta_in = np.cos(theta_in)
else:
    theta_in = np.degrees(theta_in)


# Write to file

file_path = "output.dat"
data = np.column_stack((theta_in, b_out, difCrossSec_Ruth))

if mott == "true":
    data = np.column_stack((data, difCrossSec_Mott))

np.savetxt(file_path, data, delimiter=",")


print("Data has been written to", file_path)