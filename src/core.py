from operator import ge
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
diracProton = parameters.get('diracProton')
formFactor = parameters.get('formFactor')
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



def scattering_differential_Dirac_Proton(theta, difCrossSec_Recoil, kinEn, massTarget):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    # Q = q**2 
    Q = -(2*massTarget*kinEn**2*(1-np.cos(theta)))/(massTarget+kinEn*(1-np.cos(theta)))

    difCrossSec_diracProton = difCrossSec_Recoil * (1-((Q)/(2*massTarget)*np.tan(theta/2)**2))

    return difCrossSec_diracProton



def scattering_differential_Form_Factor(theta, difCrossSec_Recoil, kinEn, massTarget):

    # Q = q**2 
    Q = -(2*(massTarget*5.6175e26)*(kinEn/(e*1e9))**2*(1-np.cos(theta)))/((massTarget*5.6175e26)+(kinEn/(e*1e9))*(1-np.cos(theta))) #Natural units (GeV)^2

    a = 0.71 # Experimental constant 0.71 GeV
    
    Form_Factor = ( 1 / (1+(Q/a**2)) )**2 # dipole

    magneticMoment = 2.79

    difCrossSec_formFactor = difCrossSec_diracProton * Form_Factor

    Ge = Form_Factor #Electric form factor 

    Gm = magneticMoment * Form_Factor # Magnetic form factor
    
    # Lorentz Invariant quantity
    Tau = -Q/(4*(massTarget*5.6175e26)**2) # Natural Units (mass in GeV)
    
    # Rosenbluth Formula
    #difCrossSec_formFactor = difCrossSec_diracProton * (( Ge**2 + Tau*Gm**2 )/(1+Tau) + 2*Tau*Gm**2*np.tan(theta/2)**2)
    
    return difCrossSec_formFactor



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

    if mott:
        difCrossSec_Mott = scattering_differential_Mott(theta_in, difCrossSec_Ruth, kinEn, massTarget) # Mott correction cross section.

    if recoil:   
         difCrossSec_Recoil = scattering_differential_Recoil(theta_in, difCrossSec_Mott, kinEn, massTarget) # Target recoil correction cross section.

    if diracProton:
        difCrossSec_diracProton = scattering_differential_Dirac_Proton(theta_in, difCrossSec_Recoil, kinEn, massTarget) # Dirac Proton correction cross section.

    if formFactor:
        difCrossSec_formFactor = scattering_differential_Form_Factor(theta_in, difCrossSec_Recoil, kinEn, massTarget) # Form Factor correction cross section.
        


# Write to file

file_path = "output.dat"

header = 'theta'
data = np.degrees(theta_in)

if impactParameter:
    header += ',b_out'
    data = np.column_stack((data, b_out))

if cross_section_variable in ['cos', 'theta', 'omega']:
    header += ',difCrossSec_Ruth'
    data = np.column_stack((data, difCrossSec_Ruth))

    if mott:
        header += ',difCrossSec_Mott'
        data = np.column_stack((data, difCrossSec_Mott))

    if recoil:
        header += ',difCrossSec_Recoil'
        data = np.column_stack((data, difCrossSec_Recoil))
        
    if diracProton:
        header += ',difCrossSec_diracProton'
        data = np.column_stack((data, difCrossSec_diracProton))
        
    if formFactor:
        header += ',difCrossSec_formFactor'
        data = np.column_stack((data, difCrossSec_formFactor))


np.savetxt(file_path, data, delimiter=",", header = header, comments="")

print("Data has been written to", file_path)