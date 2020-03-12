#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Created on Monday Oct 23 2017
@author: jonas
"""
#Implementação para calcular o 'parâmetro de ordem r saturado' vs. o 'fator de acoplamento k' 

from numpy import *
from scipy.stats import cauchy
from scipy.integrate import odeint

r_sat = []
k_list = linspace(0,10,100)
mean, gamma, N = 0, 2, 5000
kc = 2*gamma 							#That's the critical coupling for the Cauchy-Lorentz distribution. 
sticks = 200
t = linspace(0,50,sticks)					#Time vector

# Error parameters for ODEINT:
atol=1e-4
rtol=1e-4		 

def Jacobian(theta,t,omega,K,N):				#The Jacobian
	A,B = sin(theta), cos(theta)
	jac = diagflat((K/N)*(B*sum(B)+A*sum(A)))
	return jac

def theta_dot(theta,t,omega,K,N):			#Kuramoto Model Equation
	    A,B = sin(theta), cos(theta)
	    return omega + (K/N)*(B*sum(A)-A*sum(B))    				

for it in range(len(k_list)): 
								#Número de Osciladores;	
	K = k_list[it]						#Coupling Factor    		

	theta0 = random.uniform(0,2*pi,N)			#Phases

	omega = cauchy.rvs(mean, gamma, N)			#Frequencies

	if it == 10: 
		print Jacobian(theta0,t,omega,K,N)
	
	theta = odeint(theta_dot, theta0, t, args=(omega,K,N)
			,Dfun= Jacobian,rtol=rtol,atol=atol)	#Matriz gerada pelo integrador númerico, contendo os tetas.

	S1 = [sum(cos(theta[i])) for i in range(sticks)]
	d1 = array([i**2 for i in S1])

	S2 = [sum(sin(theta[i])) for i in range(sticks)]
	d2 = array([i**2 for i in S2])
#r(t) calculation:
	r = (1.0/N)*sqrt(d1 + d2)
	#Calcula a média dos últimos 1000
	r_med = [(1./25.)*(sum(r[j] for j in range(sticks-25,sticks)))]
	r_sat.append(r_med)
	fh = open('kr-jac-cauchy', 'a')
	fh.write('{}\n'.format(r_med[0]))
	print 'iteração',it, 'K', K , 'r' , r_med

f = piecewise(k_list, [k_list < 4 , k_list >= 4], [lambda k_list: 0, lambda k_list: (1-(kc/k_list))**(0.5) ])

#Visualização:
import matplotlib.pyplot as plt
plt.title('Order Parameter vs Coupling Factor - (Jacobian)')
plt.plot( k_list,f, k_list, r_sat,'ro')
plt.xlabel('k')
plt.xlim(0,10)
plt.ylabel('r')
plt.ylim(0,1.25)
plt.savefig('rk-jac-cauchy')
plt.show()



