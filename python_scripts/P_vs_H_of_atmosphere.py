#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Author: Chris Jing 
Date: 2018-05-29

Description: this program plots a graph of pressure vs. altitude, and reads off the scale height
"""


import numpy as np 
import math 
import matplotlib.pyplot as plt 

# Constants: (i = sea level)
T_i = 273 + 15 
P_i = 101000  # atmospheric pressure in Pascal
rho_i = 1.225 # density of air in kg/m^3
m = 0.02897   # mean molar mass of air in kg/mol 
g = 9.81      # gravity in m/s^2 
R = 8.314     # gas constant in J / mol / K 
gamma = 7/5   # heat capacity ratio for dry air (C_p/C_v)  
scale_height_pressure = P_i / math.e 
# height at which the pressure is 1/e of the sea level value

N = 10000 # maximum height in meter
dz = 1  # change in height is 0.01 meter 
scale_height = 0

"""
+----------------------------------------------------------------------+
| Q1. 1. Using the barymetric equation and the dry adiabatic lapse     |
| rate, make a simple finite difference calculation to find the        |
| rate of the change of pressure and density with altitude. Start at   |
| 15 C, 101 kPa. Plot a graph of pressure vs. altitude, and read off   |
| the scale height (where the pressure is 1/e of the sea level value)? |
+----------------------------------------------------------------------+

"""

# Create arrays to store temperature, pressure, density 
# rate of change of pressure with altitude, rate of change of density with altitude 
# data at different altitude and make initializations
T = np.zeros(int(N/dz));    T[0] = T_i 
P = np.zeros(int(N/dz));    P[0] = P_i 
rho = np.zeros(int(N/dz));  rho[0] = rho_i
 
dP_per_dz = np.zeros(int(N/dz)) 
drho_per_dz = np.zeros(int(N/dz)) 
height = np.linspace(0, N, num = int(N/dz)) 

for i in range(int(N/dz - 1)): 
    dT = (-m * g * (1 - 1/gamma) / R) * dz
    dP = (-m * g * P[i] / (R * T[i])) * dz 
    
    T[i+1] = T[i] + dT
    P[i+1] = P[i] + dP 
    
    dP_per_dz[i+1] = dP / dz 
    rho[i+1] = dP_per_dz[i+1] / -g
    drho_per_dz[i] = (rho[i+1] - rho[i])/dz
    
    # check to see if the current altitude has reached the scale height 
    if(P[i+1] <= scale_height_pressure and scale_height == 0):
        scale_height = (i + 1) * dz 
        T_scale_height = T[i+1]

# Plot the figure for the pressure vs. altitude
plt.figure(1) 
plt.plot(height, P/1000, 'b', label='Atmospheric Pressure vs. Height')
plt.plot(height, np.ones(int(N/dz))*(scale_height_pressure/1000), 'k--', label='scale height')
plt.text(0, (scale_height_pressure/1000) + 2, r'scale height = %.2f m' %scale_height, fontsize=10)
plt.text(0, (scale_height_pressure / 1000) - 5, r'Pressure at scale height = P_0/e = %.2f kPa' %(scale_height_pressure / 1000), fontsize=10)
plt.ylabel('Pressure in kiloPascal')
plt.xlabel('Height in meter') 
plt.title('Pressure vs. Height Plot for the Atmosphere')

# Plot the figure for the rate of change in pressure vs. altitude 
plt.figure(2) 
plt.plot(height[1:int(N/dz)], dP_per_dz[1:int(N/dz)], 'r', label='Rate of change of Pressure with altitude vs. Altitude')

plt.axvline(x= scale_height, ls = '--', color = 'k')
plt.text(scale_height - 4500, dP_per_dz[1], r'Scale Height = %.2f m' %scale_height, fontsize=10)
plt.ylabel('Rate of change in Pressure with Altitude in Pa/m')
plt.xlabel('Height in meter') 
plt.title('Rate of change of Pressure with altitude vs. Altitude')

# Plot the figure for the rate of change in air density vs. altitude 
plt.figure(3) 
plt.plot(height[1:int(N/dz)-1], drho_per_dz[1:int(N/dz)-1]*1000, 'g', label='Rate of change of air density with altitude vs. Altitude')

plt.axvline(x= scale_height, ls = '--', color = 'k')
plt.text(scale_height - 4500, drho_per_dz[1]*1000, r'Scale Height = %.2f m' %scale_height, fontsize=10)
plt.ylabel('Rate of change in Air Density with Altitude in g/m^3 /m')
plt.xlabel('Height in meter') 
plt.title('Rate of change of Air Density with altitude vs. Altitude\n')




# Plot the figure for the altitude vs. pressure 
plt.figure(4) 
plt.plot(P/1000, height, 'b', label='Height vs. Pressure')
plt.plot(P/1000, np.ones(int(N/dz))*scale_height + 500, 'k--', label='scale height')
plt.text(scale_height_pressure / 1000, scale_height, r'scale height = %.2f m' %scale_height, fontsize=10)
plt.text(scale_height_pressure / 1000, scale_height + 1000, r'Pressure at scale height = P_0/e = %.2f kPa' %(scale_height_pressure / 1000), fontsize=10)
plt.xlabel('Pressure in kiloPascal')
plt.ylabel('Height in meter') 
plt.title('Height vs. Pressure Plot for the Atmosphere')

# Plot the figure for the altitude vs. temperature
plt.figure(5) 
plt.plot(T - 273, height, 'r', label='Height vs. Temperature')
plt.plot(T - 273, np.ones(int(N/dz))*scale_height + 500, 'k--', label='scale height')
plt.text(T_scale_height - 273, scale_height, r'scale height = %.2f m' %scale_height, fontsize=10)
plt.text(T_scale_height - 273, scale_height + 1000, r'Temperature at scale height = %.2f C' %(T_scale_height - 273), fontsize=10)
plt.xlabel('Temperature in Celsius')
plt.ylabel('Height in meter') 
plt.title('Height vs. Temperature Plot for the Atmosphere')

plt.savefig("Pressure_vs_Height_of_Atmosphere_Plot")
