__authors__ = "Zachary Klouchnikov and Hannah Semple"

# HEADER

"""
IMPORTS
"""
import numpy as np
import matplotlib.pyplot as plt

"""
FUNCTIONS
"""

"""
MAIN
"""
"Constants"
M = 200  # Number of grid points
TARGET = 1e-6  # Target accuracy

"Initializing Arrays"
t = np.zeros((M + 1, M + 1), dtype = float)
t_prime = np.zeros((M + 1, M + 1), dtype = float)

"Boundary Conditions"
t[0 : 50, 0] = np.linspace(0, 5, 50) #AB boundary
t[50, 0 : 30] = np.linspace(5, 7, 30) #BC boundary
t[50 : 150, 30] = 7.0 #CD boundary
t[150, 0 : 30] = np.linspace(5, 7, 30) #DE boundary
t[150 : 200, 0] = np.linspace(5, 0, 50) #EF boundary
t[200, 0 : 80] = np.linspace(0, 10, 80) #FG boundary
t[0 : 200, 80] = 10.0 #GH boundary
t[0, 0 : 80] = np.linspace(0, 10, 80) #HA boundary

omega = 0.6 # Relaxation factor

plt.figure()
# Main loop
delta = 1.0
while delta > TARGET:

    t_prime = np.copy(t)
    # Calculate new values of the potential
    for x in range(1, M):
        for y in range(1, M):
            t[x, y] = ((1 + omega) / 4.0) * (t[x + 1, y] + t[x - 1, y] + t[x, y + 1] + t[x, y - 1]) - (omega * t[x, y])

    t[0 : 51, 0] = np.linspace(0, 5, 51) #AB boundary
    t[50, 0 : 31] = np.linspace(5, 7, 31) #BC boundary
    t[50 : 151, 30] = 7.0 #CD boundary
    t[150, 0 : 31] = np.linspace(5, 7, 31) #DE boundary
    t[150 : 201, 0] = np.linspace(5, 0, 51) #EF boundary
    t[200, 0 : 81] = np.linspace(0, 10, 81) #FG boundary
    t[0 : 201, 80] = 10.0 #GH boundary
    t[0, 0 : 81] = np.linspace(0, 10, 81) #HA boundary
    t[51 : 150, 0 : 30] = 0.0
    t[:, 81 : ] = 0.0

    # Calculate maximum difference from old values
    delta = np.max(np.abs(t - t_prime))
    print(delta)

    plt.clf()
    plt.contourf(np.transpose(t, axes = (1, 0)), origin = 'lower', cmap = 'gray')
    plt.colorbar(label = 'Temperature ($C^\circ$)')
    plt.draw()
    plt.pause(0.01)

print("Complete!")
