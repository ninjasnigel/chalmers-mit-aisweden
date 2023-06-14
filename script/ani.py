import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define the system parameters
m1, m2, L1, L2, g = 1.0, 1.0, 1.0, 1.0, 9.81

def equations(Y, t, m1, m2, L1, L2, g):
    """The function that describes the ODE system."""
    theta1, z1, theta2, z2 = Y
    c, s = np.cos(theta1-theta2), np.sin(theta1-theta2)
    theta1_dot = z1
    z1_dot = (m2*g*np.sin(theta2)*c - m2*s*(L1*z1**2*c + L2*z2**2) - (m1+m2)*g*np.sin(theta1)) / L1 / (m1 + m2*s**2)
    theta2_dot = z2
    z2_dot = ((m1+m2)*(L1*z1**2*s - g*np.sin(theta2) + g*np.sin(theta1)*c) + m2*L2*z2**2*s*c) / L2 / (m1 + m2*s**2)
    return theta1_dot, z1_dot, theta2_dot, z2_dot

# Set initial conditions
Y0 = [np.pi/2, 0.0, np.pi/2, 0.0]

# Create an array of time values to solve at
t = np.linspace(0, 20, 1000)

# Solve the ODE system
Y = odeint(equations, Y0, t, args=(m1, m2, L1, L2, g))

# Create a figure and axis to plot on
fig, ax = plt.subplots()

# Create two line objects for the rods of the pendulum
line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], 'r-', lw=1)

# Define the x and y limits of the plot
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

# This function initializes the animation
def init():
    line.set_data([], [])
    trace.set_data([], [])
    return line, trace

# This function is called for each frame of the animation
def animate(i):
    thisx = [0, L1*np.sin(Y[i, 0]), L1*np.sin(Y[i, 0]) + L2*np.sin(Y[i, 2])]
    thisy = [0, -L1*np.cos(Y[i, 0]), -L1*np.cos(Y[i, 0]) - L2*np.cos(Y[i, 2])]
    line.set_data(thisx, thisy)
    trace.set_data(thisx[2:], thisy[2:])
    return line, trace

# Create the animation
ani = animation.FuncAnimation(fig, animate, np.arange(1, len(Y)), interval=25, blit=True, init_func=init)

# Show the animation
plt.show()
