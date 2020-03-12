  #!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 01:21:03 2017

@author: jonas
"""
from numpy import *

#Implementação para calcular theta(t): 

N = 1000
K = 3.5
ispace = 1000
theta0 = random.uniform(0,2*pi,N)
mu, sigma = 0, 1.0					# Média e Desvio Padrão, respectivamente;
omega = random.normal(mu,sigma,N)

def theta_dot(theta,t,omega,K,N): 			# Função que define o modelo de Kuramoto; 
    A,B = sin(theta), cos(theta)
    return omega + (K/N)*(B*sum(A)-A*sum(B))

    
from scipy.integrate import odeint 			# Chamada do integrador numérico;

t = linspace(0,100,ispace)				 # Vetor do tempo 

theta = odeint(theta_dot, theta0, t, args=(omega,K,N)) #Matriz gerada pelo integrador númerico, contendo os tetas.

#Visualização;
	
import matplotlib.pyplot as plt

osc = random.randint(N, size = 25) #Vetor contendo 25 osciladores escolhidos aleatoriamente entre os N osciladores. 

for i in osc:
	plt.plot(t, theta[:,i])
	#plt.scatter(t, theta[:,1], s= .75, color = 'blue')
	#plt.scatter(t, r, s= .75, color='red')
	#plt.xlim(0,30)
	#plt.ylim(0,1)

plt.title('theta versus t ')
plt.xlabel('t')
plt.xlim(0,100)
plt.ylabel('theta')
plt.savefig('theta-t.pdf')
plt.show()





