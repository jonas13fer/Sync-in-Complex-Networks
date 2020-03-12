# -*- coding: utf-8 -*-
"""
Algorithm wrtiten to simulate the dynamic of a cat-cortex using the Kuramoto model. 
Author: Oliveira, Jonas. 
Federal University of Goiás, Catalão - Goiás. 
Graduate Program in Modelling and Optimization 
"""

from numpy import*
import matplotlib.pyplot as plt
import networkx as nx 
from scipy.integrate import odeint

M = loadtxt("cat-cortex.csv", unpack="True") #wheighted connection matrix
N, mu,sigma, K, ispace = len(M),0.0,0.5,0.01,1000
t = linspace(0,100,ispace)          #time
theta0 = random.uniform(-pi,pi,N)  #initial phases
omega = random.uniform(-0.5,0.5,N)

def f(teta,t,omega,K):
    A, B = sin(teta), cos(teta)
    C, D = dot(M,sin(teta)), dot(M,cos(teta))    	
    return omega + (K)*((B*C) - (A*D))

theta = odeint(f, theta0, t, args=(omega,K))

freq=[]
for j in range(ispace):
    freq.append(f(theta[j],t,omega,K))  
    
freq=array(freq)

M = loadtxt("cat-cortex-adj.csv", unpack="True")

theta_T = theta.transpose() 

S1=dot(M,cos(theta_T))
S2=dot(M,sin(theta_T))

r_i = sqrt((S1**2)+(S2**2))

r = (sum(r_i,axis=0))/(sum(M))

osc = random.randint(N, size = 65)  

plt.subplot(121)
for l in osc: 
    plt.plot (t,theta[:,l],label="O_{0}".format(l))
plt.title('Phases Cat Cortex')
plt.xlim(0,10)
plt.xlabel('t')
plt.ylabel('$\Theta$')

plt.subplot(122)
for m in osc: 
    plt.plot (t,freq[:,m],label="O_{0}".format(m))
plt.title('Freqs. Cat Cortex')
plt.xlim(0,10)
plt.xlabel('t')
plt.ylabel('$\omega$')

plt.show()

plt.plot(t,r)

plt.show()
