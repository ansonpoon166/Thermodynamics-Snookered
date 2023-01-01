#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 15:32:14 2020

@author: ansonpoon
"""
import numpy as np
import scipy as sp
import Ball as b
import pylab as pl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import simulation as S
#%%
"To run the simulation"
sim = S.Simulation(10, 1, 0.05, [0,0], [0,0] , con_R = 10)
sim.run(1000, True, False, False, False)

#%%
"to create a collection of graphs,including the kinetic energy, distance from centre and the velocity distribution"
sim = S.Simulation(10, 1, 0.05, [0,0], [0,0] , con_R = 10)
sim.run(1000, False, True, False, False)

#%%
"Task 11 pressure_area graph"
pressure = []
area_inverse = []
container_radius = [30, 40, 50, 60, 70, 80]
for ab in container_radius:
    sim = S.Simulation(10, 1, 0.05, [0,0], [0,0] , con_R = ab)
    each_pressure, temperature = sim.run(1000, False, False, False, True)
    pressure.append(each_pressure)

for i in range(0, len(container_radius), 1):
    area_inverse.append(1/(sp.pi*(container_radius[i]**2)))

print(pressure)
    
plt.plot(area_inverse, pressure, 'rx')
def func(v, a):
    return a*v

guesses = 0.1
para, cov =sp.optimize.curve_fit(func, area_inverse, pressure, p0 = guesses)
x_func = np.linspace(0, 4e-4, num=1000, endpoint=True)
print("the param is",para)
plt.plot(x_func, func(x_func, *para), 'r')
plt.xlabel("inverse area (m2)")
plt.ylabel("pressure(Pa)")
plt.title("Pressure against area")
print ("the temperature is", temperature, "K")
plt.grid()
plt.show()

#%%
"pressure_temp graph"
pressure=[]
temperature_list=[]
velocity_factor = [1, 2, 3 ,4 ,5]
for i in velocity_factor:
    sim = S.Simulation(10, 1, 0.05, [0,0], [0,0], con_R = 10, v_factor=i)
    each_pressure, temperature = sim.run(1000, False, False, False, True)
    pressure.append(each_pressure)
    temperature_list.append(temperature)

plt.plot(temperature_list, pressure, 'rx')
def func(v, a):
    return a*v
guesses=0.1
para, cov =sp.optimize.curve_fit(func, temperature_list, pressure, p0 = guesses)
x_func = np.linspace(0, 2.5e24, num=100*100, endpoint=True)
print("the param is",para)
plt.plot(x_func, func(x_func, *para), 'r')
plt.xlabel("Temperature(K)")
plt.ylabel("Pressure(Pa)")
plt.title("Pressure against temperature")
plt.grid()
plt.show()

#%%
"pressure against number"
number_root = [5, 7, 9, 11, 13, 15]
pressure=[]
number = []
for i in number_root:
    sim = S.Simulation(i, 1, 0.05, [0,0], [0,0], con_R = 10, v_factor=1)
    each_pressure, temperature = sim.run(1000, False, False, False, True)
    pressure.append(each_pressure)

def func(v, a):
    return a*v

for i in number_root:
    number.append(i**2)

plt.plot(number, pressure, 'rx')
guesses=0.1
para, cov =sp.optimize.curve_fit(func, number, pressure, p0 = guesses)
x_func = np.linspace(0, 256, num=100*10, endpoint=True)
print("the param is",para)
plt.plot(x_func, func(x_func, *para), 'r')
plt.xlabel("Number of balls")
plt.ylabel("Pressure(Pa)")
plt.title("Pressure against number of balls")
print ("The temperature is", temperature, "K.")
plt.grid()
plt.show()

#%%
"Temperature against area"
temperature_list=[]
area_inverse = []
container_radius = [30, 40, 50, 60, 70, 80]
for ab in container_radius:
    sim = S.Simulation(10, 1, 0.05, [0,0], [0,0], con_R = ab)
    each_pressure, temperature = sim.run(1000, False, False, False, True)

    temperature_list.append(temperature)

for i in range(0, len(container_radius), 1):
    area_inverse.append(1/(sp.pi*(container_radius[i]**2)))

plt.plot(area_inverse, temperature_list, 'rx')
def func(v, a):
    return a*v

plt.xlabel("inverse area (m2)")
plt.ylabel("Temperature(K)")
plt.title("Temperature against area")

plt.grid()
plt.show()

#%%

"Temperature against number"
number_root = [5, 7, 9, 11, 13, 15]

temperature_list=[]
number = []
for i in number_root:
    sim = S.Simulation(i, 1, 0.05, [0,0], [0,0] , con_R = 10, v_factor=1)
    each_pressure, temperature = sim.run(1000, False, False, False, True)

    temperature_list.append(temperature)

def func(v, a):
    return a*v

for i in number_root:
    number.append(i**2)

plt.plot(number, temperature_list, 'rx')

plt.xlabel("Number of balls")
plt.ylabel("Temperature(K)")
plt.title("Temperature against number of balls")
plt.grid()
plt.show()

#%%
"Velocity distribution"
sim = S.Simulation(10, 1, 0.05, [0,0], [0,0] , con_R = 10, v_factor=1)
sim.run(10, False, False, True, False)



    