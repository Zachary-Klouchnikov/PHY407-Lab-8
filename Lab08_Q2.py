__authors__ = "Zachary Klouchnikov and Hannah Semple"

# This program simulates the vibration of a string fixed at both ends
# using the finite difference time domain (FTCS) method. It visualizes the
# displacement of the string over time.

"""
IMPORTS
"""
import numpy as np
import matplotlib.pyplot as plt

from collections.abc import Callable

plt.rc('text', usetex = True)
plt.rc('font', family = 'serif')  

"""
PART A) AND B)
"""
"Constants"
C = 1.0 # ms^-1
L = 1.0 # m
D = 0.1 # m
SIGMA = 0.3 # m
V = 100.0 # ms^-1

"Initializing Arrays"
h = 10e-6 # Step size
x = np.linspace(0, L, 100, dtype = float) # Position array
a = L / (len(x) - 1) # Spatial step size
phi = np.zeros_like(x, dtype = float) # Displacement array

# Velocity array
psi = ((C * x * (L - x)) / (L ** 2)) * np.exp((-(x - D) ** 2) / (2 * 
                                                                 SIGMA ** 2)) 

"Main FTCS Loop"
new_phi = np.zeros_like(phi, dtype = float)
new_psi = np.zeros_like(psi, dtype = float)

plt.figure()
t = 0.0
frame_counter = 0
while t < 0.02:
    for i in range(1, len(x) - 1):
        new_phi[i] = phi[i] + h * psi[i]
        new_psi[i] = psi[i] + h * (V ** 2 / a ** 2) * (phi[i + 1] + 
                                                       phi[i - 1] - 2 * phi[i])
    
    # Update boundary conditions
    new_phi[0] = 0.0
    new_phi[-1] = 0.0
    new_psi[0] = 0.0
    new_psi[-1] = 0.0

    phi = new_phi.copy()
    psi = new_psi.copy()

    t += h
    frame_counter += 1

    if frame_counter % 20 == 0:
        # Plotting the string displacement
        plt.clf()
        plt.plot(x, phi, color = 'Teal')

        # Labels
        plt.title("Displacement of a Vibrating String", fontsize = 12)
        plt.xlabel("Position Along String $(m)$", fontsize = 12)
        plt.ylabel("Displacement $(m)$", fontsize = 12)

        # Limits
        plt.xlim(0, L)
        plt.ylim(-0.001, 0.001)

        plt.grid()

        # if frame_counter in [200, 600, 800, 900, 1000]:
            # plt.savefig(f'Figures\\Displacement of a Vibrating String at t = {t:.3f}s.pdf')
        plt.pause(0.01)
