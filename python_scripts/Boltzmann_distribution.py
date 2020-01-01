#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Author: Chris Jing 
Date Created: 2018-05-17

Instruction: To run this program, open this file in a Python IDE such as Spyder 
             and simply click F5 or the "Run file" icon. 
"""

"""
1. Recreate the velocity plots shown in class for air and helium at 20 C.
What does a “probability density” expressed in s/m mean?
Investigate the relationship between these distributions and the speed of sound. 
Calculate numerically the rms velocities for 20 C air and helium and find from 
these values the speed of sound in these gases.
Show your working.
"""

import numpy as np 
import math 
import matplotlib.pyplot as plt


# constants
k_b   = 1.38 * 10**(-23) # Boltzmann Constant 
m_air = 0.02897 # molar mass of air in kg
m_He  = 0.00402  # molar mass of Helium gas in kg 
m_Ne  = 0.02018  # molar mass of Neon gas in kg 
R     = 8.3145 # gas constant in J/(K*mol) 

 
speed_max = 2500 # maximum speed on the scale 
speed = np.linspace(0, speed_max, num = speed_max + 1) 
T = 293   # room temperature in K 
avogadro = 6.22 * 10**23 # Avogadro's Number 
C_v = (5/2) * R # heat capacity for diatomic gas under constant volume 
C_p = (7/2) * R # heat capacity for diatomic gas under constant pressure
gamma = 7/5 # heat capacity ratio for dry air (C_p/C_v)  

prob_density_air = np.zeros(speed_max + 1)
prob_density_He = np.zeros(speed_max + 1)
prob_density_Ne = np.zeros(speed_max + 1)


for i in range(speed_max + 1):
    
    term_1_air = math.sqrt(((m_air/avogadro)/(2*math.pi*k_b*T))**3) 
    term_2_air = 4*math.pi*(i**2)
    term_3_air = math.exp(-(m_air/avogadro)*(i**2)/(2*k_b*T))
    
    term_1_He = math.sqrt(((m_He/avogadro)/(2*math.pi*k_b*T))**3) 
    term_2_He = 4*math.pi*(i**2)
    term_3_He = math.exp(-(m_He/avogadro)*(i**2)/(2*k_b*T))
    
    term_1_Ne = math.sqrt(((m_Ne/avogadro)/(2*math.pi*k_b*T))**3) 
    term_2_Ne = 4*math.pi*(i**2)
    term_3_Ne = math.exp(-(m_Ne/avogadro)*(i**2)/(2*k_b*T))
    
    prob_density_air[i] = term_1_air * term_2_air * term_3_air
    prob_density_He[i] = term_1_He * term_2_He * term_3_He
    prob_density_Ne[i] = term_1_Ne * term_2_Ne * term_3_Ne

# Root mean square velocities
rms_air = (3 * k_b * T / (m_air/avogadro))**(0.5)
rms_He = (3 * k_b * T / (m_He/avogadro))**(0.5)
rms_Ne = (3 * k_b * T / (m_Ne/avogadro))**(0.5)
v_sound = rms_air * (gamma/3)**(0.5)

print("The root mean square speed of air = %.2f m/s" %rms_air)
print("The root mean square speed of He  = %.2f m/s" %rms_He)
print("The root mean square speed of Ne  = %.2f m/s" %rms_Ne)


plt.figure(1) 
plt.plot(speed, prob_density_air, 'r', label='Air')
plt.plot(speed, prob_density_He, 'b', label='Helium')
plt.plot(speed, prob_density_Ne, 'g', label='Neon')

plt.vlines(x=rms_air, ymin = 0, ymax= 0.002, linestyle = '-.', color = 'k') 
plt.text(rms_air + 10, 0.002, r'r.m.s. speed of air = %.f m/s'%rms_air, fontsize=10)
plt.vlines(x=rms_He, ymin = 0, ymax= 0.002, linestyle = '-.', color = 'k') 
plt.text(rms_He + 10, 0.001, r'r.m.s. speed of He = %.f m/s'%rms_He, fontsize=10)
plt.vlines(x=rms_Ne, ymin = 0, ymax= 0.002, linestyle = '-.', color = 'k') 
plt.text(rms_Ne + 10, 0.0005, r'r.m.s. speed of Ne = %.f m/s'%rms_Ne, fontsize=10)

plt.xlabel('Speed (m/s)')
plt.ylabel('Probability density (s/m)') 
plt.title('Maxwell-Boltzmann Distribution for Gases') 


legend = plt.legend(loc='upper right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')

plt.savefig("Boltzmann_Distribution_Plot")
