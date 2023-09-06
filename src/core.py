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
procedure = parameters.get('proc')
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
rosenbluth = parameters.get('rosenbluth')
cross_section_variable = parameters.get('CrossSecVariable')
bMin = parameters.get('bMin')
bMax = parameters.get('bMax')
nProj = parameters.get('nProj')
detectorDistance = parameters.get('detectorDistance')

# Outputs Nuclear mass in atomic mass units (u).
element = periodictable.elements[zTarget]
massTarget = element.mass
massTarget = massTarget * 0.932808457 # Converts mass to GeV


# If projectile is electron
if zProj != 0:
    massProj = periodictable.elements[zProj].mass * 0.932808457
else:
    zProj = 1
    massProj = electron_mass


kinEn = kinEn*1e9 # Converts energy of incoming particles from eV to GeV.
D = (zProj*zTarget*alpha/kinEn) # Minimum distance between incident particles and target


# Functions

def impact_parameter(theta, D):
    #Returns impact parameter when given the scattering angle.
    
    return D/(2*np.tan(theta/2))



def scattering_angle(bparam, D):
    # Returns scattering angle when given the impact parameter.

    return 2*np.arctan(D/(2*bparam))



def closest_distance(theta, D):
    # Returns closest distance that projectile reaches in relation to the target
    
    return D/2 * ( 1 + 1/(np.sin(theta/2)) )



def scattering_differential_Ruth(theta, D, cross_section_variable):
    # Returns differential scattering impact when given the scattering angle.

    if cross_section_variable == 'cos':
        difCrossSec_Ruth = (2*pi*D**2/(1-np.cos(theta))**2)

    if cross_section_variable == 'theta':
        difCrossSec_Ruth = (D**2*pi*np.cos(theta/2)/(4*np.sin(theta/2)**3))
    
    if cross_section_variable == 'omega':
        difCrossSec_Ruth = D**2/(16*np.sin(theta/2)**4)

    return difCrossSec_Ruth



def scattering_differential_Mott(theta, difCrossSec_Ruth, kinEn, massTarget):
    # Applies mott correction factor

    difCrossSec_Mott = difCrossSec_Ruth * np.cos(theta/2)**2 

    return difCrossSec_Mott 



def scattering_differential_Recoil(theta, difCrossSec_Mott, kinEn, massTarget):
    # Applies Recoil correction factor

    difCrossSec_Recoil = difCrossSec_Mott * (1/(1+(((1-np.cos(theta))*kinEn)/(massTarget*c**2))))

    return difCrossSec_Recoil



def scattering_differential_Dirac_Proton(theta, difCrossSec_Recoil, kinEn, massTarget):
    # Applies Dirac proton correction factor

    # Q = q**2 
    Q = -(2*massTarget*kinEn**2*(1-np.cos(theta)))/(massTarget+kinEn*(1-np.cos(theta)))
    
    difCrossSec_diracProton = difCrossSec_Recoil * (1-((Q)/(2*massTarget)*np.tan(theta/2)**2))

    return difCrossSec_diracProton 



def scattering_differential_Form_Factor(theta, difCrossSec_Recoil, kinEn, massTarget):
    # Applies Form Factor correction factor

    # Q = q**2
    Q = -(2*massTarget*kinEn**2*(1-np.cos(theta)))/(massTarget+kinEn*(1-np.cos(theta))) 
    a = 0.71 # Experimental constant 0.71 GeV 

    Form_Factor = ( 1 / (1+(Q/a)) )**2 # dipole
 
    difCrossSec_formFactor = difCrossSec_Recoil * Form_Factor

    return difCrossSec_formFactor



def scattering_differential_Rosenbluth(theta, difCrossSec_Recoil, kinEn, massTarget):
    # Applies Form Factor correction factor

    # Q = q**2
    Q = -(2*massTarget*kinEn**2*(1-np.cos(theta)))/(massTarget+kinEn*(1-np.cos(theta))) 
    a = 0.71 # Experimental constant 0.71 GeV 

    Form_Factor = ( 1 / (1+(Q/a)) )**2 # dipole
    magneticMoment = 2.79
    Ge = Form_Factor #Electric form factor 
    Gm = magneticMoment * Form_Factor # Magnetic form factor
    
    # Lorentz Invariant quantity
    Tau = -Q/(4*massTarget**2)

    # Rosenbluth Formula
    Rosenbluth = (( Ge**2 + Tau*Gm**2 )/(1+Tau) + 2*Tau*Gm**2*np.tan(theta/2)**2)
    difCrossSec_Rosenbluth = difCrossSec_Recoil * Rosenbluth

    return difCrossSec_Rosenbluth



# Calculations

# Calculates Impact parameter
if procedure == 'thvsb':
    theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.
    
    # Converts theta to radians if necessary
    if angUnit == 'degrees':
        theta_in = np.radians(theta_in)

    data = np.degrees(theta_in)
    
    b_out = impact_parameter(theta_in, D) 
    
    header = 'theta,b_out'
    data = np.column_stack((data, b_out))


# Calculates diferential cross section
if procedure == 'xsec':
    header = 'theta'
    theta_in = np.linspace(angStart,angEnd,1000)[1:] # Scattering angle input.
    
    # Converts theta to radians if necessary
    if angUnit == 'degrees':
        theta_in = np.radians(theta_in)

    data = np.degrees(theta_in)
    
    header += ',difCrossSec_Ruth'
    difCrossSec_Ruth = scattering_differential_Ruth(theta_in, D, cross_section_variable) #Differential scattering cross section.
    data = np.column_stack((data, difCrossSec_Ruth))
    
    if mott:
        header += ',difCrossSec_Mott'
        difCrossSec_Mott = scattering_differential_Mott(theta_in, difCrossSec_Ruth, kinEn, massTarget) # Mott correction cross section.
        data = np.column_stack((data, difCrossSec_Mott))
        
    if recoil:
        header += ',difCrossSec_Recoil' 
        difCrossSec_Recoil = scattering_differential_Recoil(theta_in, difCrossSec_Mott, kinEn, massTarget) # Target recoil correction cross section.
        data = np.column_stack((data, difCrossSec_Recoil))
        
    if diracProton:
        header += ',difCrossSec_diracProton'
        difCrossSec_diracProton = scattering_differential_Dirac_Proton(theta_in, difCrossSec_Recoil, kinEn, massTarget) # Dirac Proton correction cross section.
        data = np.column_stack((data, difCrossSec_diracProton))

    if formFactor:
        header += ',difCrossSec_formFactor'
        difCrossSec_formFactor = scattering_differential_Form_Factor(theta_in, difCrossSec_Recoil, kinEn, massTarget) # Form Factor correction cross section.
        data = np.column_stack((data, difCrossSec_formFactor))

    if rosenbluth:
        header += ',difCrossSec_rosenbluth'
        difCrossSec_rosenbluth = scattering_differential_Rosenbluth(theta_in, difCrossSec_Recoil, kinEn, massTarget) # Rosenbluth correction cross section.
        data = np.column_stack((data, difCrossSec_rosenbluth))


# Calculates 3D outgoing angle and detector position given beam of particles
if procedure == 'beam':
    b_random = np.random.uniform(bMin, bMax, nProj)
    angle = np.array([scattering_angle(b, D) for b in b_random])
    position = detectorDistance * np.sin(angle)
    rotation = np.random.uniform(0, 2 * np.pi, nProj)
    data = np.column_stack((np.degrees(angle), position, np.degrees(rotation)))
    header = 'theta,position,phi'
    

if procedure == 'bvsd':
    b_in = np.linspace(bMin,bMax,1000)[1:] # Impact parameter input converted to theta
    theta_in = scattering_angle(b_in, D)
    closestDistance = closest_distance(theta_in, D)
    header = 'b_in,closestDistance'
    data = np.column_stack((b_in,closestDistance))

file_path = "output.dat"

# Write to file
np.savetxt(file_path, data, delimiter=",", header = header, comments="")


print("Data has been written to", file_path)