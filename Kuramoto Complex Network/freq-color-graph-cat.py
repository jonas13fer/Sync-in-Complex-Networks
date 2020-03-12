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
N, K, ispace = len(M),0.5,50
t = linspace(0,50,ispace)          #time
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

cat_cortex = loadtxt("cat-adj.csv", unpack="True")
G = nx.DiGraph(cat_cortex)

degrees = [val for (node, val) in G.degree()]

ent = dict(G.in_degree(G.nodes()))

sai = dict(G.out_degree(G.nodes()))

node_size = []

for i in range(len(degrees)):
    node_size.append(degrees[i]*10)


import matplotlib.animation as animation

cmap=plt.cm.YlOrRd
fig, ax = plt.subplots()
#fig = plt.figure()
pos = nx.spring_layout(G)

it=0

def update(it):
    fig.clear()
    nx.draw(G,pos=pos,node_size=node_size,node_color = freq[it],   with_labels=False,width=0.15,alpha=0.75,cmap=cmap)
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=-pi, vmax=pi))
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    it=it+1
    print("Iteration ",it)
    return

ani = animation.FuncAnimation(fig, update, frames=ispace, interval=200, repeat=False)

ani.save('color-freq.mp4', writer='ffmpeg')

