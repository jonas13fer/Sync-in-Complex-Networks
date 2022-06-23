import numpy as np 

N, K = 50, 0.65
theta0 = np.random.uniform(-np.pi,np.pi,N)  #initial phases
omega = np.random.uniform(0.0,np.pi,N)
omega = 0.0
t = np.linspace(0,50,1000)

def kuramoto_model(theta,t,omega,K): 			# Função que define o modelo de Kuramoto; 
    A,B = np.sin(theta), np.cos(theta)
    return omega + (K/N)*(B*sum(A)-A*sum(B))

def rungekutta4th(f, y0, t, args=()):
    t_n = len(t)
    y = np.zeros((t_n,len(y0)))
    y[0] = y0
    for step in range(t_n - 1):
        h = t[step + 1] - t[step]
        k1 = f(y[step], t[step], *args)
        k2 = f(y[step] + k1 * h / 2., t[step] + h / 2., *args)
        k3 = f(y[step] + k2 * h / 2., t[step] + h / 2., *args)
        k4 = f(y[step] + k3 * h / 1., t[step] + h / 1., *args)
        y[step + 1] = y[step] + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
    return y 
    
from scipy.integrate import odeint 		

'''Runge-Kutta 4h order against ODEINT'''

theta_rk4th = rungekutta4th(kuramoto_model,theta0,t,args=(omega,K))
theta_rk4th = theta_rk4th % 2*pi
theta_odeint = odeint(kuramoto_model,theta0,t,args=(omega,K))
theta_odeint = theta_odeint % 2*pi

import matplotlib.pyplot as plt

osc = np.random.randint(N, size = 50)  #Vetor contendo osciladores escolhidos aleatoriamente entre os N osciladores. 

markersize=7
markeredgewidth=0.5

'''Error plot'''
plt.subplot(121)
for m in osc: 
    plt.scatter(t,abs(theta_odeint[:,m]-theta_rk4th[:,m]))
plt.title("ODEINT vs. Runge-Kutta4")
plt.xlim(0,50)
plt.xlabel('t')
plt.ylabel('$Error$')

plt.subplot(122)
for m in osc: 
    plt.scatter(t,theta_rk4th[:,m])
plt.title("Runge-Kutta4")
plt.xlim(0,50)
plt.xlabel('t')
plt.ylabel('$\Theta$')

plt.show()

