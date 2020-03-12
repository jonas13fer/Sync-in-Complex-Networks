# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 00:58:04 2019

@author: Oliveira, Jonas
"""
from scipy.integrate import odeint
import numpy as np 

a=0.6 #alpha
b=2.18 #beta
g=470 # gamma
N=25
y0 = np.random.normal(1,0.1,N) #initial conditions
x0 = np.random.normal(1,0.1,N)
z0 = np.random.normal(1,0.1,N)
i0 = [x0,y0,z0] # initial conditions

t = np.linspace(0,10000,50000) #time array 

def circuit(y,t,a,b):
    w1 = a*y[0] + y[2]
    w2 = y[2] - f(y[1])
    w3 = - y[0] - b*y[1]
    return [w1,w2,w3]

def f(y): 
    '''linear piece-wise diode resistance with slope 'g' '''
    return (g/2)*(abs(y+(1/g)) - abs(y-(1/g)))    

sol = odeint(circuit,i0,t,args=(a,b)) #odeint integration

import matplotlib.pyplot as plt

plt.subplot(311)
plt.plot(t,sol[:,0],color="red",label="x")
plt.legend(loc=2)

plt.subplot(312)
plt.plot(t,sol[:,1],color="orange",label="y")
plt.legend(loc=2)

plt.subplot(313)
plt.plot(t,sol[:,2],color="green",label="z")
plt.legend(loc=2)

plt.show()

plt.plot(sol[:,2],sol[:,0])

plt.show()
