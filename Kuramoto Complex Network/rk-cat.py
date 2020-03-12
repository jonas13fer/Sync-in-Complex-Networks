# -*- coding: utf-8 -*-

from numpy import*
import matplotlib.pyplot as plt
import networkx as nx 
from scipy.integrate import odeint

M = loadtxt("cat-cortex.csv", unpack="True") #wheighted connection matrix
N, mu,sigma, ispace = len(M),0.0,0.5,700
K_list = linspace(-0.001,0.2,50)
t = linspace(0,700,ispace)          #time
theta0 = random.uniform(-pi,pi,N)  #initial phases
omega = random.uniform(-0.5,0.5,N)

def f(teta,t,omega,K,M):
    A, B = sin(teta), cos(teta)
    C, D = dot(M,sin(teta)), dot(M,cos(teta))    	
    return omega + (K)*((B*C) - (A*D))

T = 300
delta = len(t)-T

def C(theta):
    theta = theta[T:delta+T]
    thetaT = theta.transpose()

    A,B = sin(theta),cos(theta)
    C,D = dot(cos(thetaT),cos(theta)),dot(sin(thetaT),sin(theta))
    E,F = dot(sin(thetaT),cos(theta)),dot(sin(thetaT),cos(theta))
    return (sqrt((C+D)**2+(E-F)**2))/delta

r_link = []
r_med = []

for K in K_list: 

    M = loadtxt("cat-cortex.csv", unpack="True") #wheighted connection matrix    

    theta = odeint(f, theta0, t, args=(omega,K,M))
    
    """
    #complex network order parameter calculation

    M = loadtxt("cat-adj.csv", unpack="True")

    theta_T = theta.transpose() 

    S1=dot(M,cos(theta_T))
    S2=dot(M,sin(theta_T))

    r_i = sqrt((S1**2)+(S2**2))

    r = (sum(r_i,axis=0))/(sum(M))

    r_med.append(sum(r)/len(r))
    """

    #r_link calculation
    r_link.append(sum(C(theta))/(N*(N-1)))

    #traditional order parameter calculation
    
    S1 = [sum(cos(theta[i])) for i in range(ispace)]
    d1 = array([i**2 for i in S1])

    S2 = [sum(sin(theta[i])) for i in range(ispace)]
    d2 = array([i**2 for i in S2])

    r = (1.0/N)*sqrt(d1 + d2)
    
    r_med.append(sum(r)/len(r))


plt.semilogx(K_list,r_med,'o')

plt.semilogx(K_list,r_link,'d')

plt.show()

    
