# -*- coding: utf-8 -*-
"""
ENPH 257 Thermodynamics Simulation Project

Author: Chris Jing 
Date Created: 2018-05-17

Desription: this program simulate the heat distribution of an aluminum rod when 
            10 Watts of thermal power is applied at one end
"""
import numpy as np 
import math 
import matplotlib.pyplot as plt


# Shared constants: 
length = 0.3   # length of the rod in meter 
radius = 0.01  # radius of the rod in m 
x_section_area   = math.pi * radius**2 

k_Al   = 235 # conductivity of aluminum in Watts / K / m 
c_Al   = 921.096 # heat capacity of aluminum in J / kg K 
rho_Al = 2700  # density of aluminum in kg / m^3 

sigma   = 5.67 * 10**(-8) # Stefan-Boltzmann constant in W/m^2/K^4 
epsilon = 0.09               # assumed emissivity for a perfect blackbody 
k_c     = 7.17               # convection coefficient in W / m^2 / K

"""
+----------------------------------------------------------------------+
| Q1. Put 10 W of thermal power into one end of a 30 cm long, 2 cm     |
| diameter aluminum rod, painted matt black. Assume the ambient        |
| temperature is a uniform 20 C, and the convection constant kc = 5    |
| W/m2/K. In your submission, include your code and a properly labelled| 
| graph of steady-state temperature versus position on the rod.        |
+----------------------------------------------------------------------+

"""

# Constants:
T_amb    = 293      # ambient temperature in Kelvin 
Power_in = 10      # Thermal power into the heated rod in W / m^2
Power_state = 1   # power turned on

# Variables (user change this) 
dx = 0.01             # slice length in m 
dt = 0.01              # time step in second 
N  = int(length / dx) # number of slices 
runtime = 200        # in second (user change this to estimate steady-state)
power_modulation_period = 60*25 # the time period to turn on/off the power supply 

surface_area = 2 * math.pi * radius * dx 

# Set the initial temperatures for the rod as the ambient temperature
T = np.ones(N) * T_amb
 
T_0 = np.zeros(int(runtime / dt) + 1); T_0[0] = T_amb # 1.25 cm 
T_1 = np.zeros(int(runtime / dt) + 1); T_1[0] = T_amb # 8.25 cm
T_2 = np.zeros(int(runtime / dt) + 1); T_2[0] = T_amb # 15.25 cm
T_3 = np.zeros(int(runtime / dt) + 1); T_3[0] = T_amb# 22.25 cm
T_4 = np.zeros(int(runtime / dt) + 1); T_4[0] = T_amb # 29.25 cm

time = np.linspace(0, runtime, num = int(runtime/dt)+1) 

# Create arrays to store position and new temperature data 
position = np.linspace(0, length, num=N)
new_T = np.zeros(N)
cycle = 0

# Loop through the rod over time to update the temperatures 
for t in range(int(runtime / dt) + 1):
    
    # confirm whether to turn on/off the power supply 
    if (cycle >= int(power_modulation_period / dt)):
        cycle = 0
        if(Power_state ==1): 
            Power_in = 0
            Power_state = 0

        else:
            Power_in = 10
            Power_state = 1
    
        
        
    for i in range(N): # i th slice of the rod 
        # temperature change that applies to all slices 
        
        # heat loss via thermal radiation 
        dT_radi = (-2 * epsilon * sigma * (T[i]**4 - T_amb**4) / (c_Al * rho_Al * radius)) * dt 
            
        # heat loss via convection 
        dT_conv = (- 2 * k_c * (T[i] - T_amb) / (c_Al * rho_Al * radius)) * dt 
            
        if i == 0: # first slice (heated end)
            
            # heat loss via conduction for the heated end 
            dT_cond = (-k_Al * (T[i] - T[i+1]) / (c_Al * rho_Al * dx**2)) * dt
            
            # heat gain via power addition
            dT_power = Power_in  * dt / (c_Al * rho_Al * x_section_area * dx) 
            
            # effect of net temperature change 
            new_T[i] = T[i] + dT_cond + dT_power + dT_radi + dT_conv 
                
        elif i == N - 1: # last slice (cool end)
            
            # heat loss via conduction for the cool end 
            dT_cond = (k_Al * (T[i-1] - T[i]) / (c_Al * rho_Al * dx**2)) * dt
            
            # effect of net temperature change 
            new_T[i] = T[i] + dT_cond + dT_radi + dT_conv 
                
        else: # middle slices
            # heat change from thermal diffusion 
            double_differential = (T[i-1] - 2*T[i] + T[i+1])/(dx**2) 
            dT_diff = (double_differential * k_Al / (c_Al * rho_Al)) * dt
    
            # effect of net temperature change 
            new_T[i] = T[i] + dT_diff + dT_radi + dT_conv 
       
    T = np.array(new_T) 
    
    T_0[t] = T[int((0.0125/length)* N)] 
    T_1[t] = T[int((0.0825/length)* N)] 
    T_2[t] = T[int((0.1525/length)* N)] 
    T_3[t] = T[int((0.2225/length)* N)] 
    T_4[t] = T[int((0.2925/length)* N)] 
    
    cycle += 1
    


# Plot the figure for the rod temperature distribution after a long time 
plt.figure(1) 
plt.plot(time, T_0 - 273, 'r', label='T @ 1.25 cm')
plt.plot(time, T_1 - 273, 'b', label='T @ 8.25 cm')
plt.plot(time, T_2 - 273, 'g', label='T @ 15.25 cm')
plt.plot(time, T_3 - 273, 'k', label='T @ 22.25 cm')
plt.plot(time, T_4 - 273, 'y', label='T @ 29.25 cm')

plt.xlabel('time elasped (second)')
plt.ylabel('temperature at the position in rod (Celsius)') 
plt.title('Temperature vs Time Plot for the Rod\nwith Power Modulation Period of %.f seconds' %power_modulation_period) 
plt.savefig("Al_Rod_Heat_Distribution_Diagram")

#plt.ylim(T_amb - 273, )
#plt.xlim(0, runtime)

legend = plt.legend(loc='upper right', shadow=False)
frame = legend.get_frame()
frame.set_facecolor('0.90')

print("The average temperature of the rod after %.f seconds " %runtime + 
      " = %.2f degree C\n" %(np.average(new_T) - 273)) 
