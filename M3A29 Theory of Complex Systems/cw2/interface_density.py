import numpy as np
import matplotlib.pyplot as plt

def interface_density(u, v, N, t):
    x = 0.5
    xlist=[x]
    for i in range(t):
        x1 = x + 4/N
        x2 = x - 4/N
        x3 = x
        p1 = (1 - x)**2*(1 - u)*v/2
        p2 = x*(x - 2/N)*(1-v)*u/2
        p3 = 1 - p2 - p1
        x = np.random.choice([x1, x2, x3], p=[p1, p2, p3])
        xlist.append(x)
    return xlist

u = 0.5
v = 0.25
alpha = 2*(1 - u)*v
beta = 2*(1 - v)*u
eq9 = 1/(1+np.sqrt(beta/alpha))
t = 200000

x_axis = range(t+1)
y_axis1 = np.asarray(interface_density(u, v, 10000, t))
y_axis2 = np.asarray(interface_density(u, v, 20000, t))
y_axis3 = np.asarray(interface_density(u, v, 30000, t))
for i in range(19):
    y_axis1 += np.asarray(interface_density(u, v, 10000, t))
    y_axis2 += np.asarray(interface_density(u, v, 20000, t))
    y_axis3 += np.asarray(interface_density(u, v, 30000, t))
y_axis1 = y_axis1/20
y_axis2 = y_axis2/20
y_axis3 = y_axis3/20
y_axis = [eq9 for i in x_axis]

plt.figure()
plt.plot(x_axis, y_axis1, label="<x(t)>, N=10000")
plt.plot(x_axis, y_axis2, label="<x(t)>, N=20000")
plt.plot(x_axis, y_axis3, label="<x(t)>, N=30000")
plt.plot(x_axis, y_axis, label="Eq. (9)")
plt.title("Plot of x(t) averaged over 20 times and Eq. (9) against Time")
plt.xlabel("Time")
plt.ylabel("")
plt.legend()
plt.grid()
plt.show()