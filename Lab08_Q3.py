__authors__ = "Zachary Klouchnikov and Hannah Semple"

# This code answers Q3 of Lab8 for PHY407, where we solve Burgers' equation using a modified leapfrog method
# for this nonlinear wave equation. We show the results at t=0s, 0.5s, 1s, and 1.5s

"""
IMPORTS
"""
import numpy as np
import matplotlib.pyplot as plt

"""
SOLVING
"""
#constants and initialisation
e = 1
dx = 0.02
dt = 0.005
Lx = 2 * np.pi  #m
Tf = 2  #s

beta = e * dt / dx
Nx = int(Lx / dx)  # needs to be integer value
Nt = int(Tf / dt)

burger = np.zeros((Nt, Nx))  # Nt many time rows, Nx many distance columns

# initial + boundary conds
burger[0] = np.sin(xs)
burger[0][0] = 0  # x=0 condition
burger[0][-1] = 0  # x=Lx condition
xs = np.linspace(0,Lx,Nx)  #array of x values in m

# one forward euler step in the time domain
for i in range(1, Nx-1):
    burger[1][i] = burger[0][i] - (beta/4) * ((burger[0][i+1])**2 - (burger[0][i-1])**2)
burger[1][0] = 0  #initial conds
burger[1][-1] = 0

# leapfrog using eq 8 from lab8 manual
for j in range(1, Nt-1):
    for i in range(1, Nx-1):
        burger[j+1][i] = burger[j-1][i] - (beta/2) * (burger[j][i+1]**2 - burger[j][i-1]**2)
    
    # Apply boundary conditions at each time step
    burger[j+1][0] = 0
    burger[j+1][-1] = 0

"""
PLOTTING
"""
t = [0, 0.5, 1, 1.5]  #time slices
inds = [0,100,200,300]  #j indices that correspond to t=[0, 0.5, 1, 1.5]

plt.figure(figsize=(8,7))
for p in range(1,5):
    plt.subplot(2,2,p)
    plt.plot(xs,burger[inds[p-1]], color='teal')
    time = t[p-1]
    plt.title('Burgers\' Solution at t={}s'.format(time), fontsize=14)
    plt.xlabel('Distance [m]', fontsize=14)
    plt.ylabel('Velocity Disturbance [m/s]', fontsize=14)
    plt.grid()

plt.tight_layout()
# plt.savefig('burger.pdf', bbox_inches='tight')
plt.show()
