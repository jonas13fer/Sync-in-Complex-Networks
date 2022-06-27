# -*- coding: utf-8 -*-
'''
Vinte e sete de junho de dois mil e vinte dois. 

Solução do modelo de Kuramoto utilizando Runge-Kutta de quarta ordem. Via de regra, utilizávamos o ODEINT do scipy.integrate.
Contudo, o ODEINT não permitia, por exemplo, insirir uma perturbação em um dado tempo da integração.
Não há nada de revolucionário nesta abordagem, use como quiser, e caso encontre algum erro, comunique-me. 

Cordialmente,

Jonas Oliveira (jonas.oliveira@inpe.br)
Computação Aplicada, INPE
'''

import numpy as np 

N, K = 100, 1.1 #number of nodes and coupling strength
theta0 = np.random.uniform(-np.pi,np.pi,N)  #initial phases
omega = np.random.uniform(-0.5,0.5,N) #intrinsic frequencies

t = np.linspace(0,50,1000) #time vector

'''Kuramoto Model all-to-all '''
def kuramoto_model(theta,t,omega,K): 			
    A,B = np.sin(theta), np.cos(theta)
    return omega + (K/N)*(B*sum(A)-A*sum(B))

'''Fourth Order Runge-Kutta'''
def rungekutta4th(f, y0, t, args=()):
    t_n = len(t)
    y = np.zeros((t_n,len(y0))) #matrix with the times series of each oscillator 
    y[0] = y0
    
    '''It randomly sets the step where the pertubation will be inserted'''
    perturbation_step = np.random.randint(t_n)
    
    for step in range(t_n - 1):
        h = t[step + 1] - t[step]
        k1 = f(y[step], t[step], *args)
        k2 = f(y[step] + k1 * h / 2., t[step] + h / 2., *args)
        k3 = f(y[step] + k2 * h / 2., t[step] + h / 2., *args)
        k4 = f(y[step] + k3 * h / 1., t[step] + h / 1., *args)
        
        y[step + 1] = y[step] + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)  
        
        '''Pertubation insertion'''
        if step == perturbation_step:
           print("The perturbation was inserted at the {0}th step.".format(step+1))      
           y[step + 1] = y[step + 1] + np.random.normal(0,1.0,N)
           
    return y 
    
from scipy.integrate import odeint 		

'''Runge-Kutta 4h order against ODEINT'''

theta_rk4th = rungekutta4th(kuramoto_model,theta0,t,args=(omega,K))
theta_rk4th = np.mod(theta_rk4th,2*np.pi)

'''If you want to test Runge-Kutta 4th precision, just compare it to ODEINT'''
theta_odeint = odeint(kuramoto_model,theta0,t,args=(omega,K))
theta_odeint = np.mod(theta_odeint,2*np.pi)

"""Frequency's calculation """ 
freq=[]
for j in range(len(t)):
    freq.append(kuramoto_model(theta_rk4th[j],t,omega,K))      
freq=np.array(freq)

import matplotlib.pyplot as plt

osc = np.random.randint(N, size = 50)  #It chooses 50 oscillators out of N to plot the phases and frequencies

markersize=7
markeredgewidth=0.5

'''Frequencies and Phases Plot'''
plt.subplot(121)
for m in osc: 
    plt.plot(t,freq[:,m])
plt.title("Frequencies")
plt.xlim(0,50)
plt.xlabel('t')
plt.ylabel('$\dot{\Theta}$')

plt.subplot(122)
for m in osc: 
    plt.plot(t,theta_rk4th[:,m])
plt.title("Phases")
plt.xlim(0,50)
plt.xlabel('t')
plt.ylabel('$\Theta$')

plt.show()

''' Order parameter r calc'''
S1 = [sum(np.cos(theta_rk4th[i])) for i in range(len(t))]
d1 = np.array([i**2 for i in S1])

S2 = [sum(np.sin(theta_rk4th[i])) for i in range(len(t))]
d2 = np.array([i**2 for i in S2])

r = (1.0/N)*np.sqrt(d1 + d2)

#Visualização:	

plt.plot(t, r)
plt.xlabel('t')
plt.xlim(0,50)
plt.ylabel('r (t)')
plt.ylim(0,1)
plt.show()

