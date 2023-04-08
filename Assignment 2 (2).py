#!/usr/bin/env python
# coding: utf-8

# Kashish Bhatt 
# 
# 
# Assignment 2 

# ### Question 1: 
# 

# Calculate the equilibrium (reversal) potentials for the following ions at a temperature of 20◦C.
# 
# 
# (a) K+, [K+]out= 5 mmol/L, [K+]in= 150 mmol/L.
# 
# 
# (b) Na+, [Na+]out= 150 mmol/L, [Na+]in= 15 mmol/L.
# 
# 
# (c) Cl−, [Cl−]out= 125 mmol/L, [Cl−]in= 10 mmol/L.
# 
# 
# (d) Ca2+, [Ca2+]out= 2 mmol/L, [Ca2+]in= 0.0002 mmol/L.
# 
# 
# (e) What is the effect of increasing the temperature to 25◦ C on these equilibrium potentials.
# 
# 
# (f) Plot graphs of ENa, EK, ECa and ECl as a function of temperature in the range 10◦ C
# to 40◦ C.

# <hr style="border:0.7px solid blue"> 
#     Here we will employ the Nerst Equation: 
# 
# $$
#     E_i=\frac{RT}{zF}\,ln\frac{X_{out}}{X_{in}}
# $$
# 
# 

# <hr style="border:0.7px solid blue">
# At 25◦ C
#                $$ \frac{RT}{zF}\ $$ 
# 
# 
# will be slightly different depending on the z value which is the valency of ion. 
# 
#  
# has a valency of +1  $${K^+}$$ 
# 
# has a valency of +1  $${Na^+}$$       
# 
# has a valency of 2+  $${Ca^{2+}}$$   
# 
# has a valency of -1.  $${Cl^-}$$  

# <hr style="border:0.7px solid blue">
# 
# So using the Formula we get: 
# 
# | Molecule      | Forumla                         |
# | --------------|---------------------------------|
# | $${Na^+}$$    | $${58}log\frac{X_{out}}{X_{in}}$$|
# | $${K^+}$$     | $${58}log\frac{X_{out}}{X_{in}}$$|
# | $${Cl^-}$$    |$${-58}log\frac{X_{out}}{X_{in}}$$|
# | $${Ca^2+}$$   |$${29}log\frac{X_{out}}{X_{in}}$$ |
# 
# 
# 

# <hr style="border:0.9px solid red">
# 
# (a)   K<sup>+</sup> 
# 
# 
# $$
#           {K^+}_{out} = 5 mmol/L
# $$          
# $$          {K^+}_{in} = 150 mmol/L.
# $$

# Formula: 
#     
# $$
#      E_i={58}log\frac{X_{out}}{X_{in}}
# $$
# 
# 
# 
# $$
#      E_i={58}*log(\frac{5 mmol/L}{150 mmol/L})
# $$
# 
# Answer: 
# 
# $$
#      E_i={58}*log({0.0333})
# $$
# 
# $$
#      E_i = -{86 mV}
# $$
# 
#     

# <hr style="border:0.5px solid red">
# 
# (b)   Na<sup>+</sup> 
# 
# 
# $$
#           {Na^+}_{out} = 150 mmol/L
# $$          
# $$          {Na^+}_{in} = 15 mmol/L.
# $$
# 

# Formula: 
#     
# $$
#      E_i={58}log\frac{X_{out}}{X_{in}}
# $$
# 
# 
# $$
#      E_i={58}*log(\frac{150 mmol/L}{15 mmol/L})
# $$
# 
# Answer: 
# 
# $$
#      E_i={58}*log({10})
# $$
# 
# $$
#      E_i = + {58 mV}
# $$
# 

# <hr style="border:0.5px solid red">
# 
# (c)   Cl<sup>-</sup> 
# 
# 
# $$
#           {Cl^-}_{out} = 125 mmol/L
# $$          
# $$          {Cl^-}_{in} = 10 mmol/L.
# $$
# 

# Formula: 
#     
# $$
#      E_i={-58}log\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#      E_i={-58}*log(\frac{125 mmol/L}{10 mmol/L})
# $$
# 
# Answer: 
# 
# $$
#      E_i= {-58}*log({12.5})
# $$
# 
# $$
#      E_i = -{64mV}
# $$
# 

# <hr style="border:0.5px solid red">
# 
# (d)   Ca<sup>2+</sup> 
# 
# 
# $$
#           {Ca^2+}_{out} = 2 mmol/L
# $$          
# $$          {Ca^2+}_{in} = 0.0002 mmol/L.
# $$
# 

# Formula: 
#     
# $$
#      E_i={29}log\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#      E_i={29}*log(\frac{2 mmol/L}{0.0002 mmol/L})
# $$
# 
# Answer: 
# 
# $$
#      E_i= {29}*log({10000})
# $$
# 
# $$
#      E_i = -{116}
# $$
# 

# <hr style="border:0.5px solid red">
# (f.) E<sub>Na</sub>, E<sub>K</sub>, E<sub>Ca</sub> and E<sub>Cl</sub> as a function of temperature in the range 10◦ C to 40◦ C

# The Nerst equation like stated before is 

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import pandas as pd
from math import log
from math import exp


# In[2]:


F = 96485
R = 8.314
T = np.linspace(start=283,stop=313)


# In[3]:


Na = (R*T)/(1*F)*log(150/15)*1000
K = (R*T)/(1*F)*log(5/150)*1000
Cl = (R*T)/((-1)*F)*log(125/10)*1000
Ca = (R*T)/(2*F)*log(2/0.0002)*1000

plot(T-273, K, c="maroon", label = 'K+')
plot(T-273, Ca, c="blue", label = 'Ca2+')
plot(T-273, Cl, c="purple",label = 'Cl-' )
plot(T-273, Na, c="green", label ='Na+')
plt.legend()
xlabel('Temperature')
ylabel('Equilibrium potential (mV)')
plt.show()


# <hr style="border:1px solid red">
# 
# ### Question 2

# (a) Consider an ion X+ at 20◦ C. What is the concentration relation [X+]out/ [X+]in necessary
# to maintain a resting membrane potential V = -60 mV?

# <hr style="border:0.5px solid red">
# 
# Remember: Where we are looking for X+ ion, meaning z = 1. 
# 
# Using the Nernst Equation: 
# 
# $$
#     E_i=\frac{RT}{zF}\,ln\frac{X_{out}}{X_{in}}
# $$
# 
# 

# $$
#     -60mv=\frac{{8.314 {J/K*mol}}*{293.15K}}{{1}*{96485{C/mol}}}\,ln\frac{X_{out}}{X_{in}}
# $$
# 

# $$
#     -0.06v=\frac{{8.314 {J/K*mol}}*{293.15K}}{{1}*{96485 C/mol}}\,ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     -0.06v={0.0253}*ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     \frac{-0.06v}{0.0253v}\,=ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     -2.372 =ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#    e^{ -2.372} =\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#    0.0933 =\frac{X_{out}}{X_{in}}
# $$
# 

# <hr style="border:0.5px solid red">
# 
# (b) How is this value affected when the temperature is increased to 25◦ C?

# $$
#     -60mv=\frac{{8.314 {J/K*mol}}*{298.15K}}{{1}*{96485{C/mol}}}\,ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     -0.06v=\frac{{8.314 {J/K*mol}}*{298.15K}}{{1}*{96485 C/mol}}\,ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     -0.06v={0.02569v}*ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     \frac{-0.06v}{0.02569v}\,=ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#     -2.333 =ln\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#    e^{ -2.333} =\frac{X_{out}}{X_{in}}
# $$
# 
# $$
#    0.0967 =\frac{X_{out}}{X_{in}}
# $$
# 

# <hr style="border:0.5px solid red">
# 
# (c) Plot a graph of the concentration relation [X+]out/ [X+]in necessary to maintain a resting
# membrane potential V = -60 mV as a function of the temperature in the range 10◦ C to 40◦
# C.

# In[4]:


T = np.linspace(start=283,stop=313)


# In[5]:


X = np.exp((-578.91)/(T*8.314))
plt.plot(T-273,X)
xlabel('Temperature (celsius)');
ylabel('[X+]out / [X+]in');


# <hr style="border:0.5px solid red">
# 
# ### Question 3

# Consider the following passive membrane equation: 
# 
# $$ 
#     τ \frac{dV}{dT}\,= −(V − EL) + R(I_{app})
# $$

# with:
# $$
#     V_0= {EL}
# $$
# $$
#     R = {100}MΩ
# $$
# $$
#    C = 100 pF
# $$
# $$
#     I_{app} = 0.25 nA 
# $$
# $$
#     E_L = {−60} mV 
# $$

# (a) Calculate the value of τ

# $$
# τ = R(C) 
# $$

# Before we continue lets convert picoFarads to Microfarad. 
# $$
# 1𝑝𝐹 =  1.0E-6 µF
# $$
# 
# Using that conversion: 
# 
# $$
# 100 pF = 0.0001 µF
# $$

# Therefore: 
# $$
# τ=𝑅(𝐶)
# $$
# <hr style="border:0.5px solid red">
# $$
# τ=0.0001µ𝐹(100𝑀Ω) = 0.01 Seconds
# $$

# <hr style="border:0.5px solid red">
# (b) Calculate the value of the term $R_{Iapp}$

# Formula for Resistance: 
# $$
# R_{app}= R(I_{app})
# $$    

# So here we will use the $I_{app}$: 
# $$
# R_{app}= \frac{(100)(0.25)}{1000}\ = 0.025
# $$

# <hr style="border:0.5px solid red">
# (c) Calculate, if possible, the time it takes the voltage to reach V = −50 mV .

# Given: 
# $$
# V(t) = -50 mV
# $$
# 
# We will use the simplified version of this:
# $$
# V(t) = V_∞ + E_l + (V_0+ V_∞ + E_l)e^\frac{t}{τ}\
# $$
# 
# Simplified version:
# $$
# V(t) = V_∞ + E_l + (V_∞)e^\frac{t}{τ}\
# $$
# 
# $$
# 50E^{-3}= 0.025 + (-60E^-3) + (-0.025)e^\frac{t}{τ}\
# $$
# 
# $$
# -t = τ*ln(0,6)= -0.01*(-0.51) = 5.1E^{-3}
# $$
# 
# $$
# t =5.1 msec
# $$
# 

# <hr style="border:0.5px solid red">
# (d) Calculate, if possible, the time it takes the voltage to reach V = −30 mV .

# $$
# V(t) = -30mV
# $$
# 
# We will use the simplified version of this:
# $$
# V(t) = V_∞ + E_l + (V_0+ V_∞ + E_l)e^\frac{t}{τ}\
# $$
# 
# Simplified version:
# $$
# V(t) = V_∞ + E_l + (V_∞)e^\frac{t}{τ}\
# $$
# 
# $$
# 50E^{-3}= 0.025 + (-60E^-3) + (-0.025)e^\frac{t}{τ}\
# $$
# 
# $$
# -0.2 = e^\frac{t}{τ}\
# $$
# 
# $$
# -t/τ = ln(0,6)= -τln(-0.2)
# $$

# Here log of -0.2 is not possible, and therefore the time calculation is not accesible, because the line would continue on to go to infinity but never cross the x-axis. 

# ### Question 4

# Write a code to solve numerically eq. (1) or adapt the template code provided in the course
# website. Simulate eq. (1) for the parameter provided in Question 3.

# In[58]:


from PIL import Image, ImageDraw, ImageFilter
im1 = Image.open('/Users/kashish/Desktop/Screen Shot 2022-12-10 at 11.42.48 PM.png')
im1


# The graph depicts that visual reasoning why the log function from the previous question will
# not be possible. The function of voltage and time wiill continue on to infinity but never reach
# 0
