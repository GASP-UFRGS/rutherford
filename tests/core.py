import matplotlib.pyplot as plt
import numpy as np
import sys
from card_reader import read_card, _raise_missing_card_error()
from scipy.constants import epsilon_0, pi, e

try:
    card_name = sys.argv[1]
except IndexError:
    _raise_missing_card_error()

# Constants

T, Z, z = read_card(card_name)

T = T*e # Energy of the particles in J.
k = 1/(4*pi*epsilon_0)
r0 = (k*z*Z*e**2/T)*1e15 # Minimum distance between incident particles and target in fm.

# Functions

def impact_parameter(theta, r0, angle_unit='degrees'):
    """
    Returns impact parameter when given the scattering angle.
    """

    if angle_unit == 'radians':
        b = r0/(2*np.tan(theta/2))
        return b

    elif angle_unit == 'degrees':
        b = r0/(2*np.tan(np.radians(theta)/2))
        return b

    return

def scattering_angle(b, r0, angle_unit='degrees'):
    """
    Returns scattering angle when given the impact parameter.
    """

    if angle_unit == 'radians':
        return 2*np.arctan(r0/(2*b))
    
    elif angle_unit == 'degrees':
        return np.degrees(2*np.arctan(r0/(2*b)))

    return
def scattering_differential(theta, r0, angle_unit ='degrees'):
    """
    Returns differential scattering impact when given the scattering angle.
    """
    if angle_unit == 'radians':
        dsig = np.pi*r0**2*np.cos(theta/2)
        dsigdtheta = dsig/(4*np.sin(theta/2)**3)
        return dsigdtheta
    
    if angle_unit == 'degrees':
        theta = np.radians(theta)
        dsig = np.pi*r0**2*np.cos(theta/2)
        dsigdtheta = dsig/(4*np.sin(theta/2)**3)
        return dsigdtheta
    return

# Calculations

theta_in = np.linspace(0,180,100)[1:] # Scattering angle input.
b_out = impact_parameter(theta_in, r0) # Impact parameter calculated.

b_in = np.linspace(0,30*r0,10000)[1:] # Impact parameter input.
theta_out = scattering_angle(b_in, r0) # Scattering angle calculated.

dsig_dtheta = scattering_differential(theta_in, r0) #Differential scattering cross section

# Plots

# b vs theta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(theta_in,b_out)
plt.ylabel(r'$b$ [fm]',fontsize=14)
plt.xlabel(r'$\theta$ [degrees]',fontsize=14)
plt.title(r'Impact parameter as function of scattering angle',fontsize=16)

plt.savefig('plot_b_vs_theta.png', dpi=300, bbox_inches='tight')

# theta vs b

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(b_in, theta_out)
plt.ylabel(r'$\theta$ [degrees]',fontsize=14)
plt.xlabel(r'$b$ [fm]',fontsize=14)
plt.title(r'Scattering angle as function of impact parameter',fontsize=16)

plt.savefig('plot_theta_vs_b.png', dpi=300, bbox_inches='tight')

# theta vs dsig_dtheta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(theta_in, dsig_dtheta)
plt.yscale("log")
plt.xlabel(r'$\theta [degrees]$',fontsize=14)
plt.ylabel(r'$d\sigma/d\theta$',fontsize=14)
plt.title(r'Distribution of $d\sigma/d\theta$ in function $\theta$',fontsize=16)

plt.savefig('plot_dsig_dtheta_vs_theta.png', dpi=300, bbox_inches='tight')

