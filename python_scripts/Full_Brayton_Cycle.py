# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Brayton Cycle Demo 

Author: Chris Jing 
Date Created: 2018-05-17

Description: this program plots pressure vs volumne for the Brayton Cycle
             of an ideal gas with an isentropic pressure ratio of 20
"""
import numpy as np 
import math 
import matplotlib.pyplot as plt


# Constants for the Ideal Gas used in the Brayton Cycle (User change this)
P_1 = 100000.0 # Initial pressure in Pascal 
molar_mass = 28.97 # molar mass of air in 28.97 g/mol 
n = 1000.0/28.97 # gas molecules in mol
R = 8.314 # gas constant in J / mol / K 
T_1 = 293    # Initial temperature in Kelvin (ambient temp = 20 C)
V_1 = n * R * T_1 / P_1     # Initial volume in m^3 
pressure_ratio = 20 # isentropic pressure ratio 
#heat_added = 2000000 # heat added to the system in Stage 2 -> 3 in J


# Other constants 
C_v = (5/2) * R # heat capacity for diatomic gas under constant volume 
C_p = (7/2) * R # heat capacity for diatomic gas under constant pressure
gamma = 7/5 # heat capacity ratio for dry air (C_p/C_v)  

T_3 = 1000 + 273 # max operating temperature of the system in Kelvin 


"""
+---------------+
| Brayton Cycle |
+---------------+

"""

# Process 1 -> 2 (Isentropic Compression Stage) 

P_2 = P_1 * pressure_ratio 

P_12 = np.linspace(P_1, P_2, num = 50)
V_12 = ((P_1 * V_1**gamma) / P_12)**(1/gamma) 
V_2 = V_12[len(V_12) - 1]
T_2 = P_2 * V_2 / (n * R) 
W_12 = P_1 * V_1**gamma * ((V_2**(1-gamma) - V_1**(1-gamma))/(1-gamma)) 


# Process 2 -> 3 (Isobaric Combustion Stage)

dT_23 = T_3 - T_2 # net temperature change 
W_23  = n * R * dT_23     # work done via the heat added to the system 
dV_23 = W_23 / P_2        # net change in volume 

P_23 = np.linspace(P_2, P_2, num = 50)
V_23 = np.linspace(V_2, V_2 + dV_23, num = 50)

P_3 = P_2 
V_3 = V_23[len(V_23) - 1]
heat_added = n * C_p * (T_3 - T_2) # depends on T_3 (max operating temeprature)


# Process 3 -> 4 (Isentropic Expansion Stage)

P_4 = P_3 / pressure_ratio 

P_34 = np.linspace(P_3, P_4, num = 50)
V_34 = ((P_3 * V_3**gamma) / P_34)**(1/gamma) 
V_4 = V_34[len(V_34) - 1]

W_34 = P_3 * V_3**gamma * ((V_4**(1-gamma) - V_3**(1-gamma))/(1-gamma)) 


# Process 4 -> 1 (Isobaric Heat Rejection Stage)

V_41 = np.linspace(V_4, V_1, num = 50)
P_41 = np.linspace(P_4, P_4, num = 50)

dV_41 = V_1 - V_4 
P_4 = P_1 
W_41 = P_4 * dV_41     # work done by the system 
dT_41 = W_41 / (n * R) # net temperature change 
heat_rejected = n * (C_v + R) * dT_41

T_4 = P_4 * V_4 / (n * R) 


"""
+-------------------------+
| Efficiency Verification |
+-------------------------+

"""

theoretical_efficiency = 1 - (P_1/P_2)**(1-1/gamma) # or equivalent to 1 - T_1/T_2 
W = W_23 + W_34 + W_41 + W_12  # Net work done by the system = area inside PV locus
Q_H = n * C_p * (T_3 - T_2) # = heat added
calculated_efficency = W / Q_H 


print("theoretical_efficiency = %.2f %%" %(theoretical_efficiency*100))
print("calculated_efficiency  = %.2f %%"  %(calculated_efficency*100))

print("P_1 = %.2f kPa" %(P_1/1000))
print("V_1 = %.2f m^3" %(V_1))
print("T_1 = %.2f K\n" %(T_1))

print("P_2 = %.2f kPa" %(P_2/1000))
print("V_2 = %.2f m^3" %(V_2))
print("T_2 = %.2f K\n" %(T_2))

print("P_3 = %.2f kPa" %(P_3/1000))
print("V_3 = %.2f m^3" %(V_3))
print("T_3 = %.2f K\n" %(T_3))

print("P_4 = %.2f kPa" %(P_4/1000))
print("V_4 = %.2f m^3" %(V_4))
print("T_4 = %.2f K\n" %(T_4))

# Plot the P vs V figure for the Brayton Cycle 
plt.plot(V_12, P_12 / 1000, 'b', label='Isentropic Compression Stage 1->2')
plt.plot(V_23, P_23 / 1000, 'r', label='Isobaric Heat Addition Stage   2->3')
plt.plot(V_34, P_34 / 1000, 'g', label='Isentropic Expansion Stage     3->4')
plt.plot(V_41, P_41 / 1000, 'k', label='Isobaric Heat Rejection Stage  4->1')

plt.text(V_1 + 0.02, P_1 / 1000, r'1', fontsize=10)
plt.text(V_2 - 0.05, P_2 / 1000, r'2', fontsize=10)
plt.text(V_3 + 0.02, P_3 / 1000, r'3', fontsize=10)
plt.text(V_4 + 0.02, P_4 / 1000, r'4', fontsize=10)

plt.ylabel('Pressure (kiloPascal))')
plt.xlabel('Volume (cubic meter)')
plt.title('Pressure vs Volume Plot for the Brayton Cycle of Ideal Gas\n ' + 
          'with an isentropic pressure ratio of %.f' %pressure_ratio)

legend = plt.legend(loc='upper right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')

standard_entropy = 6.8484 # Specific Entropy of Air at Room temperature (20 C) in kJ / kg K

plt.savefig("Full_Brayton_Cycle_Plot")



