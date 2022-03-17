#!/usr/bin/env python3

# import physical constants module
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import epsilon_0, pi, e

# Constants

T = 7.7e6*e # Energy of the alpha particles in J.
k = 1/(4*pi*epsilon_0)
r0 = (k*2*79*e**2/T)*1e15 # Minimum distance between incident particles and target in fm.

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

def scattering_differential(theta, r0):
    dsig = np.pi*r0**2*np.cos(theta/2)
    dsigdtheta = dsig/(4*np.sin(theta/2)**3)
    return dsigdtheta

# Calculations

theta_in = np.linspace(0,180,100)[1:] # Scattering angle input.
b_out = impact_parameter(theta_in, r0) # Impact parameter calculated.

b_in = np.linspace(0,30*r0,10000)[1:] # Impact parameter input.
theta_out = scattering_angle(b_in, r0) # Scattering angle calculated.

theta = np.linspace(-90, 90,100)[1:]
dsig_dtheta = scattering_differential(theta, r0) #Differential scattering impact

# Plots

# b vs theta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(theta_in,b_out)
plt.ylabel(r'b [$fm$]',fontsize=14)
plt.xlabel(r'$\theta [graus]$',fontsize=14)
plt.title(r'Distribuição de $b$ em função de $\theta$',fontsize=16)

plt.savefig('grafico_b_vs_theta.png', dpi=300, bbox_inches='tight')

# theta vs b

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(b_in, theta_out)
plt.ylabel(r'$\theta [graus]$',fontsize=14)
plt.xlabel(r'b [$fm$]',fontsize=14)
plt.title(r'Distribuição de $\theta$ em função de $b$',fontsize=16)

plt.savefig('grafico_theta_vs_b.png', dpi=300, bbox_inches='tight')

# dsig_dtheta vs theta

plt.figure(figsize=(8,6), facecolor='w')
plt.plot(dsig_dtheta,theta_in)
plt.ylabel(r'$\theta [graus]$',fontsize=14)
plt.xlabel(r'$d\sigma/d\theta$',fontsize=14)
plt.title(r'Distribution of $\theta$ in function $d\sigma/d\theta$',fontsize=16)

plt.savefig('grafico_theta_vs_dsgi_dtheta.png', dpi=300, bbox_inches='tight')