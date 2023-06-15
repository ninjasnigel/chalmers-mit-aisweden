import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

list1 = [random.randint(1, 4) for _ in range(2)]
list2 = [random.randint(1, 10) for _ in range(2)]

print(list1)

# Lengths and masses
L1, L2 = list1
m1, m2 = list2

# Gravity
g = 9.81

# Initial conditions
theta1_0 = np.pi/2
theta2_0 = np.pi
theta1_dot_0 = 0.0
theta2_dot_0 = 0.0

# Time array
t = np.linspace(0, 1000, 50000)

# The system of ODEs
def dy_dt(y, t, L1, L2, m1, m2):
    theta1, theta2, theta1_dot, theta2_dot = y

    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)

    theta1_dot2 = (m2*g*np.sin(theta2)*c - m2*s*(L1*theta1_dot**2*c + L2*theta2_dot**2) - (m1+m2)*g*np.sin(theta1)) / (L1 * (m1 + m2*s**2))
    theta2_dot2 = ((m1+m2)*(L1*theta1_dot**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + m2*L2*theta2_dot**2*s*c) / (L2 * (m1 + m2*s**2))
    return theta1_dot, theta2_dot, theta1_dot2, theta2_dot2

# Integrate the ODE
y0 = np.array([theta1_0, theta2_0, theta1_dot_0, theta2_dot_0])
sol = odeint(dy_dt, y0, t, args=(L1, L2, m1, m2))
theta1, theta2 = sol[:, 0], sol[:, 1]

# Convert to Cartesian coordinates of the centers of mass
x1 = L1/2 * np.sin(theta1)
y1 = -L1/2 * np.cos(theta1)
x2 = L1 * np.sin(theta1) + L2/2 * np.sin(theta2)
y2 = -L1 * np.cos(theta1) - L2/2 * np.cos(theta2)

# Rod ends' coordinates
x1_end = L1 * np.sin(theta1)
y1_end = -L1 * np.cos(theta1)
x2_end = L1 * np.sin(theta1) + L2 * np.sin(theta2)
y2_end = -L1 * np.cos(theta1) - L2 * np.cos(theta2)

fig, ax = plt.subplots()

# Create line objects for the rods and the path of the mass centers
rod1, = ax.plot([], [], 'k-', lw=2)
rod2, = ax.plot([], [], 'k-', lw=2)
trace, = ax.plot([], [], 'r-', lw=1)

def init():
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    rod1.set_data([], [])
    rod2.set_data([], [])
    trace.set_data([], [])
    return rod1, rod2, trace,

def update(i):
    # Update rod positions
    rod1.set_data([0, x1_end[i]], [0, y1_end[i]])
    rod2.set_data([x1_end[i], x2_end[i]], [y1_end[i], y2_end[i]])
    
    # Update trace line
    trace.set_data([x1[:i+1], x2[:i+1]], [y1[:i+1], y2[:i+1]])
    
    return rod1, rod2, trace,

ani = animation.FuncAnimation(fig, update, frames=len(t), init_func=init, blit=True, interval=1)

plt.show()
