#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 14:38:55 2017

@author: henriette
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from operator import add



n = 100;          # number of penguin speed
r = 0.5;         # radius of penguin
temp_init = 1;   # temperture of penguin
T0 = 1;

L = 20;         # domain size
noise = 0.01;   # noise amplitude
dt = 0.1;       # time step size

v0 = 0.1;       # penguin speed
a = 1;          # weighting of prior orientation on heading
f = 10000;         # weighting of soft repulsion on orientation
t = 0.1;        # weighting of temperature gradient on orientation
k = 0.5;        # temp decay
tau = 0.01;     # relaxation constant for angle
s = 1;
tau1 = 0.1;     # temperature relaxation

#pos = 2*L*rand(n,2) - 0.5*L; # position of penguins     rand(n,2) creates nx2 matrix
pos = np.array(2*L*np.random.random((n,2)) - 0.5*L)
R = np.zeros([n,n])      # array for penguin penguin distance;
F = np.zeros([n,2])     # orrientation vector of penguins
S = np.zeros([n,2])            # visual for each penguin to follow



################################################################################
#                                                                              #
#                Initial equilibration (huddling) phase                        #
#                                                                              #
################################################################################

fig, ax = plt.subplots()

peng, = ax.plot(pos[:,0], pos[:,1], 'ok', markersize= 11)

ax.set_xlim(0, L)
ax.set_ylim(0, L)
plt.title('time = 0')

    
def animate(t):
    
    # penguin-penguin distance
    for i in range(1,n):
        R[:, i] = np.sqrt(np.square(pos[:, 0] - pos[i, 0]) + np.square(pos[:, 1] - pos[i, 1]))
        R[i, i] = 4*r

#    # soft core repulsion
    for i in range(1,n):
        F[i,:] = [0, 0];
        rr = R[:,i];
        I = np.where(rr < 2*r);
        I = []
        for (i,item) in enumerate(rr):
            if item < 2*r:
                I.append(i)
        if I:
            for j in I:
                F[i,:] = F[i,:] - (2*r - rr[j])**(1/2)*(pos[j,:] - pos[i,:])/rr[j]
       
    
    # huddling
    for i  in range(1,n):
        S[i,:] = [0, 0]
        rr = R[:,i]
        I = np.arange(n)
        for j in range(n):
            if I[j] == i: 
                I[j] = 0
#        I[I == i] = []
        for j in I:
            pij = pos[j, :] - pos[i, :]
            S[i,:] = S[i,:] + pij/np.linalg.norm(pij)
       
    
    # recalulcate orientation of penguins
    N_temp = S + f*F
    
    for i in range(1,n):
        N_temp[i,:] = N_temp[i,:]/np.linalg.norm(N_temp[i,:])
    
#    N_temp(isnan(N_temp)) = 0;
    
    # update position
    pos = pos + dt*v0*10*N_temp #+ 0.01*randn(size(pos))
    peng.set_data([pos[i][0] for i in range(n)], [pos[i][1] for i in range(n)])  # update the data
    return peng,

ani = animation.FuncAnimation(fig,animate, np.arange(1,500))#, init_func=init)#, interval = 1, blit = True)
plt.show()








