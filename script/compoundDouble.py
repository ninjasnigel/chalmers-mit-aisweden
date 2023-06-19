import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
L1 = 1.0  # Length of the first arm
L2 = 1.0  # Length of the second arm
m1 = 1.0  # Mass of the first arm
m2 = 1.0  # Mass of the second arm
g = 9.8   # Acceleration due to gravity

# Initial conditions (theta1, theta2, omega1, omega2)
theta1_0 = np.pi / 2
theta2_0 = np.pi / 4
omega1_0 = 0.0
omega2_0 = 0.0

# Time
dt = 0.01   # Time step
t = np.arange(0, 10, dt)

# Function to calculate the derivatives of the state variables
def derivatives(state, t):
    theta1, theta2, omega1, omega2 = state

    # Equations of motion
    alpha1 = (-g * (2 * m1 + m2) * np.sin(theta1) - m2 * g * np.sin(theta1 - 2 * theta2)
              - 2 * np.sin(theta1 - theta2) * m2 * (omega2**2 * L2 + omega1**2 * L1 * np.cos(theta1 - theta2))) \
             / (L1 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

    alpha2 = (2 * np.sin(theta1 - theta2) * (omega1**2 * L1 * (m1 + m2) + g * (m1 + m2) * np.cos(theta1)
              + omega2**2 * L2 * m2 * np.cos(theta1 - theta2))) \
             / (L2 * (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2)))

    return [omega1, omega2, alpha1, alpha2]

# Integrate the equations of motion using the Runge-Kutta method
state0 = [theta1_0, theta2_0, omega1_0, omega2_0]
state = np.zeros((len(t), 4))
state[0] = state0

for i in range(len(t) - 1):
    k1 = dt * np.array(derivatives(state[i], t[i]))
    k2 = dt * np.array(derivatives(state[i] + 0.5 * k1, t[i] + 0.5 * dt))
    k3 = dt * np.array(derivatives(state[i] + 0.5 * k2, t[i] + 0.5 * dt))
    k4 = dt * np.array(derivatives(state[i] + k3, t[i] + dt))
    state[i + 1] = state[i] + (1 / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

# Calculate the positions of the center of mass
x1 = L1 * np.sin(state[:, 0])
y1 = -L1 * np.cos(state[:, 0])
x2 = x1 + L2 * np.sin(state[:, 1])
y2 = y1 - L2 * np.cos(state[:, 1])

# Create the figure and axes
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.grid()

# Initialize the lines representing the pendulum rods and paths
line1, = ax.plot([], [], lw=2)
line2, = ax.plot([], [], lw=2)
trace1, = ax.plot([], [], 'r-', alpha=0.5)
trace2, = ax.plot([], [], 'b-', alpha=0.5)

# Initialize the center of mass positions and traces
center_mass_positions_1 = []
center_mass_positions_2 = []
trace_length = 200  # Number of frames for the trace to fade out
trace_lines_1 = []
trace_lines_2 = []

# Function to update the animation
def update(frame):
    # Update the positions of the pendulum rods
    line1.set_data([0, x1[frame]], [0, y1[frame]])
    line2.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])

    # Update the center of mass positions
    center_mass_positions_1.append([(0.5 * x1[frame]), (0.5 * y1[frame])])
    center_mass_positions_2.append([((x1[frame] + x2[frame]) / 2), ((y1[frame] + y2[frame]) / 2)])
    center_mass_xs_1, center_mass_ys_1 = zip(*center_mass_positions_1)
    center_mass_xs_2, center_mass_ys_2 = zip(*center_mass_positions_2)

    # Update the center of mass trace paths
    trace1.set_data(center_mass_xs_1[max(0, frame - trace_length):frame], center_mass_ys_1[max(0, frame - trace_length):frame])
    trace2.set_data(center_mass_xs_2[max(0, frame - trace_length):frame], center_mass_ys_2[max(0, frame - trace_length):frame])

    # Fade out the trace by reducing alpha gradually
    trace_alpha = np.linspace(1.0, 0.0, min(frame, trace_length))
    for i in range(len(trace_lines_1)):
        trace_lines_1[i].set_alpha(trace_alpha[i])
        trace_lines_2[i].set_alpha(trace_alpha[i])

    return line1, line2, trace1, trace2

# Animate the pendulum
ani = FuncAnimation(fig, update, frames=len(t), interval=dt * 1000, blit=True)

# Display the animation
plt.show()
