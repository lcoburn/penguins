# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 10:55:35 2016

@author: silke
"""


import numpy as np

import matplotlib.pyplot as plt
#import matplotlib.patches as ptch
#import matplotlib.lines as lne
import matplotlib.animation as ani



def PlotCurve(i,N,mat):
    t=i
    
    k = np.random.randint(0,N)
    l = np.random.randint(0,N)
    mat[k,l] = -mat[k,l]
    ax.imshow(mat)

    return ax



fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
N = 10
mat = np.random.choice([-1,1], [N,N])
#class matplotlib.animation.FuncAnimation(fig, func, frames=None, init_func=None, fargs=None, save_count=None, **kwargs)
anim = ani.FuncAnimation(fig, PlotCurve,fargs=(N,mat))
plt.show()
