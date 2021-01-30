from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

"""
This program uses odeint() to solve the given differential equation of a
launched projectile. The trajectory of this numerical solution is then plotted
alongside the analytical solution for comparison. Also, a series of
questions are answered using the data returned by odeint().

"""

# Define initial conditions
m = 0.1
g = 9.8
c = 0.65  # atmospheric friction constant
v0 = 10


def dydt(ya, ta):
    x = ya[0]  # x position
    y = ya[1]  # y position
    vx = ya[2]  # x velocity
    vy = ya[3]  # y velocity
    ax = -(c*vx)/m  # x acceleration from given equation
    ay = (-g - (c*vy/m))  # y acceleration from given equation
    return np.array([vx, vy, ax, ay])


t0 = 0.0  # start time
tmax = 2  # end time
steps = 100

# Time Array
ta = np.linspace(t0, tmax, steps+1)
# Initial conditions for position and velocity
theta = (5*np.pi/18)  # converted degrees to radians
y0a = [0, 0, 10*np.cos(theta), 10*np.sin(theta)]
# Solve differential equation
yM = odeint(dydt, y0a, ta)

# Plot numerical solution
plt.subplot(211)
plt.plot(yM[0:100:2, 0], yM[0:100:2, 1], 'o', label='Numerical Solution')
# Plot analytical solution
t = ta
vT = m*g/c
x = (v0*vT/g)*np.cos(theta)*(1-np.exp(-g*t/vT))
y = (vT/g)*((v0*np.sin(theta))+vT)*(1-np.exp(-g*t/vT))-vT*t
plt.plot(x, y, 'r', label='Analytical Solution')
plt.title('Launched Projectile')
plt.xlabel('Distance (m)')
plt.ylabel('Height (m)')
plt.legend()

# Find point of impact
x = yM[46, :]  # Used np.where to find the last row where the value of y>=0
d = x[0]  # Distance to impact
# Final velocity
vx = x[2]  # Final x velocity
vy = x[3]  # Final y velocity
vf = np.sqrt(vx**2+vy**2)
# Find max height reached
h = np.max(yM[:, 1])
# Find time of flight
t = (vT/-g)*(np.log(1-((d*g)/(v0*vT*np.cos(theta)))))

print('The distance d to impact is: {:.3} m'.format(d))
print('The maximum height reached is: {:.3} m'.format(h))
print('The time of flight is: {:.3} s'.format(t))
print('The velocity at the impact point is: {:.3} m/s'.format(vf))

# How does d change if the angle is increased to 70 degrees?
theta2 = (7*np.pi/18)
y0a = [0, 0, 10*np.cos(theta2), 10*np.sin(theta2)]
yM = odeint(dydt, y0a, ta)
p = yM[46, :]
d2 = p[0]

# Plot angles and distances
x = [theta, theta2]
y = [d, d2]
plt.subplot(212)
plt.plot(x, y, 'r')
plt.title('Angle vs Distance')
plt.xlabel('Launch angle (rad)')
plt.ylabel('Distance traveled (m)')
plt.tight_layout()
plt.savefig('ballistic.pdf')


