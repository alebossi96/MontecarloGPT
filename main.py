import matplotlib.pyplot as plt
import numpy as np
import math
"""
values = []
def p(theta):
    g = 0.  9
    return 1/(4*math.pi)*(1-g**2)/(1+g**2-2*g*np.cos(theta))**(1.5)
with open("output.txt", "r") as f:
  for line in f:
    values.append(float(line.strip()))
values = np.array(values)
bins = 1000
hist, bin_edges = np.histogram(values, bins=bins, density=True)

#values = values[values<10]
plt.plot(bin_edges[:bins], hist/np.max(hist), label="Measured")
step = 1e-3
x = np.arange(step,math.pi,step = step)
plt.plot(x,p(x)*np.sin(x)*2*math.pi/np.max(p(x)*np.sin(x)*2*math.pi), 'k--', label = "Theoretical")
print(np.sum(p(x)*np.sin(x))*step*2*math.pi)
print(np.sum(hist)*(bin_edges[1]- bin_edges[0]))
plt.show()
"""
import montecarlomodule as mc
for mu_s in [20]:
    res = mc.mc(0.5,mu_s)
    plt.plot(res, label = str(mu_s))
plt.legend()
plt.yscale('log')
plt.show()


