# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Brayton Cycle Demo 

Author: Chris Jing 
Date Created: 2018-05-17

Description: This program create a computer model of a Brayton Cycle and make 
             a fully-labelled PV plot for an isentropic pressure ratio of 15. 

    - Show that the area enclosed by the cycle is consistent with the theoretical efficiency formula. 
    - Make a fully-labelled plot of thermal efficiency vs. pressure ratio, 
      showing the results of the theoretical formula as a line, 
      and the results of your numerical calculation as “data points”.
"""
import numpy as np 
import matplotlib.pyplot as plt


# Constants for the Ideal Gas used in the Brayton Cycle (User change this)
P_1 = 100000 # Initial pressure in Pascal 
V_1 = 1      # Initial volume in m^3 
T_1 = 293    # Initial temperature in Kelvin (ambient temp = 20 C)
pressure_ratio = 15 
#heat_added = 2000000 # heat added to the system in Stage 2 -> 3 in J


# Other constants 
R = 8.314 # gas constant in J / mol / K 
n = P_1 * V_1 / (R * T_1) # gas molecules in mol
C_v = (5/2) * R # heat capacity for diatomic gas under constant volume 
C_p = (7/2) * R # heat capacity for diatomic gas under constant pressure
gamma = 7/5 # heat capacity ratio for dry air (C_p/C_v)  

T_3 = 1000 + 273 # max operating temperature of the system in Kelvin 

"""
+---------------+
| Brayton Cycle |
+---------------+

"""
N = 40
pressure_ratio = np.ones(N)

calculated_efficency = np.ones(N)
theoretical_efficiency = np.ones(N)

for i in range(N):
    pressure_ratio[i] = i+1
    
    

    # Process 1 -> 2 (Isentropic Compression Stage) 
    
    P_2 = P_1 * pressure_ratio[i] 
    
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
    
    P_4 = P_3 / pressure_ratio[i] 
    
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
    
    theoretical_efficiency[i] = 1 - (P_1/P_2)**(1-1/gamma) # or equivalent to 1 - T_1/T_2 
    W = W_23 + W_34 + W_41 + W_12  # Net work done by the system = area inside PV locus
    Q_H = n * C_p * (T_3 - T_2) # = heat added
    calculated_efficency[i] = W / Q_H 


plt.figure(1) 
plt.plot(pressure_ratio, calculated_efficency, 'ro', label='calculated efficency')
plt.plot(pressure_ratio, theoretical_efficiency, 'b', label='theoretical efficiency')

plt.xlabel('pressure ratio')
plt.ylabel('efficiency (%)') 
plt.title('Efficency vs. Pressure Ratio Plot') 

legend = plt.legend(loc='lower right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')

plt.savefig("Brayton_Cycle_Efficiency_Plot") 
