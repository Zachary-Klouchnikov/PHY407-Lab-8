__authors__ = "Zachary Klouchnikov and Hannah Semple"

# This program solves for the steady-state temperature distribution on a 2D
# plate with specified boundary conditions using the finite difference method.
# It implements the Gauss-Seidel iterative method with over-relaxation to
# accelerate convergence.

"""
IMPORTS
"""
import numpy as np
import matplotlib.pyplot as plt

"""
PART B)
"""
"Constants"
M = 200  # Number of grid points

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

omega = 0.0 # Relaxation factor

# Gauss-Seidel Iteration with replacement and over-relaxation
i = 0 # Iteration counter
while i != 100:

    t_prime = np.copy(t) # Store previous iteration
    
    # Calculate new values of the temperature
    for x in range(1, M):
        for y in range(1, M):
            t[x, y] = ((1 + omega) / 4.0) * (t[x + 1, y] + t[x - 1, y] + 
                                             t[x, y + 1] + t[x, y - 1]) - (
                                                 omega * t[x, y])

    # Reapply Boundary Conditions
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

    i += 1 # Increment iteration counter

"Plotting Temperature Distribution With Omega = 0.0"
plt.figure()

# Plotting temperature distribution with omega = 0.0
plt.contourf(np.transpose(t, axes = (1, 0)), origin = 'lower', cmap = 'gray')

# Labels
plt.title("Temperature Distribution With $\omega = 0.0$", fontsize = 12)
plt.xlabel("X Position $(cm)$", fontsize = 12)
plt.ylabel("Y Position $(cm)$", fontsize = 12)

plt.savefig('Figures\\Temperature Distribution With Omega = 0.0.pdf')
plt.show()

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

omega = 0.9 # Relaxation factor

# Gauss-Seidel Iteration with replacement and over-relaxation
i = 0 # Iteration counter
while i != 100:

    t_prime = np.copy(t) # Store previous iteration

    # Calculate new values of the temperature
    for x in range(1, M):
        for y in range(1, M):
            t[x, y] = ((1 + omega) / 4.0) * (t[x + 1, y] + t[x - 1, y] + 
                                             t[x, y + 1] + t[x, y - 1]) - (
                                                 omega * t[x, y])

    # Reapply Boundary Conditions
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

    i += 1 # Increment iteration counter

"Plotting Temperature Distribution With Omega = 0.9"
plt.figure()

# Plotting temperature distribution with omega = 0.9
plt.contourf(np.transpose(t, axes = (1, 0)), origin = 'lower', cmap = 'gray')

# Labels
plt.title("Temperature Distribution With $\omega = 0.9$", fontsize = 12)
plt.xlabel("X Position $(cm)$", fontsize = 12)
plt.ylabel("Y Position $(cm)$", fontsize = 12)

plt.savefig('Figures\\Temperature Distribution With Omega = 0.9.pdf')
plt.show()

"""
PART C)
"""
"Constants"
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

omega = 0.9 # Relaxation factor

plt.figure()
# Gauss-Seidel Iteration with replacement and over-relaxation
delta = 1.0 # Initial difference
while delta > TARGET:

    t_prime = np.copy(t) # Store previous iteration

    # Calculate new values of the temperature
    for x in range(1, M):
        for y in range(1, M):
            t[x, y] = ((1 + omega) / 4.0) * (t[x + 1, y] + t[x - 1, y] + 
                                             t[x, y + 1] + t[x, y - 1]) - (
                                                 omega * t[x, y])

    # Reapply Boundary Conditions
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
    plt.contourf(np.transpose(t, axes = (1, 0)), cmap = 'gray')
    plt.colorbar(label = "Temperature ($C^\circ$)")
    plt.draw()
    plt.pause(0.01)

print("Temperature at (2.5cm, 1cm): ", t[25, 10])
