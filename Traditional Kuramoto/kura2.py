#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 01:05:09 2017

@author: jonas
"""
#Implementação para calcular r(t): 

from numpy import *
N = 100
K = 2.5
ispace = 1000

theta0 = random.uniform(0,2*pi,N)
mu, sigma = 0, 0.5  # Média e Desvio Padrão, respectivamente;
omega = random.normal(mu,sigma,N)


def theta_dot(theta,t,omega,K,N): # Função que define o modelo de Kuramoto; 
    A,B = sin(theta), cos(theta)
    return omega + (K/N)*(B*sum(A)-A*sum(B))
    
from scipy.integrate import odeint # Chamada do integrador numérico;

t = linspace(0,100,ispace) # Vetor do tempo 

theta = odeint(theta_dot, theta0, t, args=(omega,K,N)) #Matriz gerada pelo integrador númerico, contendo os tetas.


#r(t) calculation:

S1 = [sum(cos(theta[i])) for i in range(ispace)]
d1 = array([i**2 for i in S1])

S2 = [sum(sin(theta[i])) for i in range(ispace)]
d2 = array([i**2 for i in S2])

r = (1.0/N)*sqrt(d1 + d2)


#Visualização:
	
import matplotlib.pyplot as plt
plt.scatter(t, r, s=5.0, marker='o', color='blue')
plt.xlabel('t')
plt.xlim(0,30)
plt.ylabel('r (t)')
plt.ylim(0,1.25)
plt.show()


