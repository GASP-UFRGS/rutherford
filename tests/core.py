import matplotlib.pyplot as plt
import numpy as np
import sys
from card_reader import read_card, _raise_missing_card_error
from scipy.constants import epsilon_0, pi, e

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

# Calculations

theta_in = np.linspace(angStart,angEnd,100)[1:] # Scattering angle input.
b_out = impact_parameter(theta_in,r0,angUnit) # Impact parameter calculated.

dsig_dtheta = scattering_differential(theta_in,r0,angUnit) #Differential scattering cross section

# Plots

# b vs theta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(theta_in, b_out)
plt.ylabel(r'$b$ [fm]',fontsize=14)
plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
plt.title('Impact parameter as function of scattering angle',fontsize=16)

plt.savefig('plot_b_vs_theta.png', dpi=300, bbox_inches='tight')

# theta vs b

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(b_out, theta_in)
plt.ylabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
plt.xlabel(r'$b$ [fm]',fontsize=14)
plt.title('Scattering angle as function of impact parameter',fontsize=16)

plt.savefig('plot_theta_vs_b.png', dpi=300, bbox_inches='tight')

# theta vs dsig_dtheta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(theta_in, dsig_dtheta)
plt.yscale("log")
plt.xlabel(r'$\theta$ [{unit}]'.format(unit=angUnit),fontsize=14)
plt.ylabel(r'$d\sigma/d\theta$',fontsize=14)
plt.title(r'Distribution of $d\sigma/d\theta$ in function $\theta$',fontsize=16)

plt.savefig('plot_dsig_dtheta_vs_theta.png', dpi=300, bbox_inches='tight')

