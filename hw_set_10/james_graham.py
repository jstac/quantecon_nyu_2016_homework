# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 12:38:14 2016

@author: James
"""

import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from numpy.linalg import norm

# parameter
s1 = 0.5
theta = 2.5
delta = 0.7
rho = 0.4

# Create the h function, which will be defined on a grid
def h_fun(n, rho):
    def g_fun(h):
        g = 2 - 1/(h + rho*n) - 1/(h + n/rho)
        return g
    h = bisect(g_fun, 0, 1)
    return h

Npoints = 1000

n_grid = np.linspace(1e-10,1,Npoints)
h_grid = np.empty([Npoints,1])

for n in range(Npoints):
    h_grid[n] = h_fun(n_grid[n],rho)



#%%  Approximate the h function with a spline

h_spl = UnivariateSpline(n_grid, h_grid,k=5)

xs = np.linspace(0, 1, 1000)
plt.plot(xs, h_spl(xs), 'g', lw=3)
plt.plot(n_grid, h_grid, 'r', lw=3)
plt.show()

#%% Compute the paths

def law_motion(n1,n2):    
    if n1 <= 0.5 and n2 <= 0.5:
        n1_next = delta*(theta/2 + (1-theta)*n1)
        n2_next = delta*(theta/2 + (1-theta)*n2)
    elif n1 >= h_spl(n2) and n2 >= h_spl(n1):
        n1_next = delta*n1
        n2_next = delta*n2
    elif n1 >= 0.5 and n2 <= h_spl(n1):
        n1_next = delta*n1
        n2_next = delta*(theta*h_spl(n1) + (1-theta)*n2)
    elif n1 <= h_spl(n2) and n2 >= 0.5:
        n1_next = delta*(theta*h_spl(n2) + (1-theta)*n1)        
        n2_next = delta*n2
    return n1_next, n2_next


def check_conv(n1_first, n2_first, n1_third, n2_third):
    first = np.array([n1_first, n2_first])     
    third = np.array([n1_third, n2_third])
    if norm(first - third) < 1e-10:
        conv = 1
    else:
        conv = 0
    return conv
         

N_iter = 100
n1 = np.empty([N_iter,1])
n2 = np.empty([N_iter,1])
n1[0] = 0.2
n2[0] = 0.3
    
# Test the simulations
for t in range(1,N_iter):
    n1[t], n2[t] = law_motion(n1[t-1],n2[t-1])
    if t >=2 and check_conv(n1[t-2],n2[t-2],n1[t],n2[t]) == 1:
        if np.isclose(n1[t],n2[t], 1e-4):
            sync = 1
        else:
            sync = 0
        break    
plt.plot(n1, 'r', lw=2,label='n1', alpha = 0.2)
plt.plot(n2, 'b', lw=2,label='n2', alpha = 0.2)
plt.legend()
plt.show()



#%%
# Run across all grid points 
Npoints = 100
N_iter = 1000
n_grid = np.linspace(1e-10,1,Npoints)
sync = np.empty(([Npoints,Npoints]))
synctime = np.empty(([Npoints,Npoints]))
for nn in range(Npoints):
    for mm in range(Npoints):
        n1 = np.empty(([N_iter,1]))
        n2 = np.empty(([N_iter,1]))
        n1[0] = n_grid[nn]
        n2[0] = n_grid[mm]        
        for t in range(1,N_iter):
            n1[t], n2[t] = law_motion(n1[t-1],n2[t-1])
            if t >=2 and check_conv(n1[t-2],n2[t-2],n1[t],n2[t]) == 1:
                if np.isclose(n1[t],n2[t], 1e-4):
                    sync[nn,mm] = 1 
                    synctime[nn,mm] = t
                else:
                    sync[nn,mm] = 0
                    synctime[nn,mm] = 0
                break
        


print('Im done!')

#%% Plot the contour stuff


# Note to self: plt.contour = contour lines; plt.contourf = filled contours 
levels = [0.5, 0.9]
c = plt.contour(n_grid,n_grid,sync, levels, linewidths=2, colors='k')
colorscheme = ('r','w') 
c = plt.contourf(n_grid,n_grid,sync, levels, colors=colorscheme) 
plt.xlabel("$n_1$")
plt.ylabel("$n_2$")
plt.show()

#%%
# Plot with time to synch
levels = np.linspace(120,350,30);
c = plt.contour(n_grid,n_grid,synctime, levels, linewidths=2)   # , levels
# colorscheme = ('r','w') 
c = plt.contourf(n_grid,n_grid,sync,levels)   #, levels, colors=colorscheme) 
plt.xlabel("$n_1$")
plt.ylabel("$n_2$")
plt.show()
