#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Author: Chris Jing 
Date: 2018-06-27

Description: this program produces a plot of the spectral irradiance 
             of the Sun’s radiation at the top of the Earth’s atmosphere
             (a) in W/(m2.nm)
             (b) in W/(m2.eV)
"""


import numpy as np
import math 
import matplotlib.pyplot as plt 


T_sun = 5800 # temperature of Sun's surface in Kelvin 
emissivity = 1 # assume perfect black body 

wavelength_min = 50 # in nm
wavelength_max = 4000 # in nm 
N = wavelength_max - wavelength_min # sample points 
irradiance = np.zeros(N+1)
wavelength = np.linspace(wavelength_min*(10**(-9)), wavelength_max*(10**(-9)), num = N+1) # in nanometers 
h = 6.626 * 10**(-34) # Planck's constant in J*s
c = 299792458 # speed of light in m/s 
k_b = 1.3806 * 10**(-23) # Boltzmann constant in J/K

radius_sun = (695700 * 10**3 + 696392 * 10**3) / 2 # radius of the Sun in m
radius_earth = 6371 * 10**3 # in meter 
radius_sun_earth_orbit = 149.6 * 10**9 - radius_earth # mean Earth-Sun Orbit radius in m
solid_angle = math.pi * (radius_sun**2)/(radius_sun_earth_orbit**2)

area = 0

for i in range(N+1):
    constant_term = 2 * h * c**2 / (wavelength[i]**5)
    irradiance[i] = solid_angle * (constant_term / (math.exp(h*c/(wavelength[i] * k_b * T_sun)) - 1)) 
    area += (irradiance[i] * 10**(-9))
  
surface_irradiance = np.zeros(N+1)
albedo = 0.25
area_surface = 0

for i in range(N+1):
    if (i+50 >= 400 and i+50 <= 900): 
        surface_irradiance[i] = irradiance[i]*(1-albedo)
        area_surface += (surface_irradiance[i] * 10**(-9))
        
    

# Plot the figure for the spectral irradiance of the Sun in unit [W m^-2 nm^-1]
plt.figure(1) 
plt.plot(wavelength*(10**(9)), irradiance / 10**(9), 'r', label='Top of Atmosphere')
plt.plot(wavelength*(10**(9)), surface_irradiance / 10**(9), 'b', label='Surface of Earth')
plt.ylabel('Spectral Irradiance (W m^-2 nm^-1)')
plt.xlabel('Wavelength (nm)') 
plt.title('Solar Irradiance of the Sun at the\ntop of Earth\'s atmospehere')
legend = plt.legend(loc='upper right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')
plt.savefig("irradiance_vs_wavelength_in_nm")

print("The solar irradiance at top of atmosphere is approximately %.2f W/(m^2)"%area)
print("The solar irradiance at surface of Earth is approximately %.2f W/(m^2)"%area_surface)


frequency = np.linspace(1, c / (200 * 10**(-9)), num = N + 1)
eVs = frequency * h / (1.6 * 10**(-19))

irradiance_eVs = np.zeros(N+1)
irradiance_f = np.zeros(N+1)

area_eVs = 0
area_f = 0
for i in range(N+1):
    constant_term = 2 * h * frequency[i]**3 / (c**2)
    irradiance_f[i] = solid_angle * (constant_term / (math.exp(h*frequency[i]/(k_b * T_sun)) - 1)) 
    irradiance_eVs[i] = (1.6 * 10**(-19)) * irradiance_f[i] / h
    area_eVs += irradiance_eVs[i] * (eVs[1] - eVs[0])
    area_f += irradiance_f[i] * (frequency[1] - frequency[0])


surface_irradiance_eVs = np.zeros(N+1)
area_surface_eVs = 0

for i in range(N+1):
    if (eVs[i] > 1.37 and eVs[i] < 3.1): 
        surface_irradiance_eVs[i] = irradiance_eVs[i]*(1-albedo)
        area_surface_eVs += surface_irradiance_eVs[i] * (eVs[1] - eVs[0]) 
          

# Plot the figure for the spectral irradiance of the Sun in unit [W m^-2 eV^-1]
plt.figure(2) 
plt.plot(eVs, irradiance_eVs, 'b', label='Solar Irradiance of the Sun at the\ntop of Earth\'s atmospehere')
plt.plot(eVs, surface_irradiance_eVs, 'r', label='Solar Irradiance of the Sun at the\n surface of Earth')
plt.ylabel('Spectral Irradiance (W m^-2 eV^-1)')
plt.xlabel('Energy (eV)') 
plt.title('Solar Irradiance of the Sun')
legend = plt.legend(loc='upper right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')
plt.savefig("irradiance_vs_energy_in_eV")

print("The solar irradiance at top of atmosphere is approximately %.2f W/(m^2)"%area_eVs)
print("The solar irradiance at surface of Earth is approximately %.2f W/(m^2)"%area_surface_eVs)

