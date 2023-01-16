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
class obj_det:
    def __init__(self, x,y,z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
d = obj_det(0,0,0,0.1)
import montecarlomodule as mc

for mu_s in [1, 5, 10, 20, 40, 80]:
    res = mc.mc(0.5,mu_s, d)
    time = np.arange(0,len(res))/1e3
    plt.plot(time, res/np.sum(res[100:]), label = "$mu_s$" +str(mu_s))

v = 29.9
mu_s = 20*(0.5)
def de_0(time):
    return 1./time**(3/2)
plt.plot(time,de_0(time)/np.sum(de_0(time)[100:]), label = "teo")
plt.xlabel("time [ns]")
plt.ylabel("amplitude")
plt.legend()
plt.xlim(-0.01, 0.15)
#plt.savefig("imgs/mu_s_lin.png")
plt.yscale('log')
plt.xlim(-0.01, time[-1])
#plt.savefig("imgs/mu_s_log.png")
plt.show()
#compara g
plt.close()
for g in [0, 0.1,0.2, 0.4,0.8,0.9]:
    res = mc.mc(g, 10, d)
    time = np.arange(0,len(res))/1e3
    plt.plot(time, res/np.sum(res[100:]), label = "g" + str(g))
plt.legend()
v = 29.9
def de_0(time):
    return 1./time**(3/2)
plt.plot(time,de_0(time)/np.sum(de_0(time)[100:]), label = "teo")
plt.xlabel("time [ns]")
plt.ylabel("amplitude")
plt.legend()
plt.xlim(-0.01, 0.15)
#plt.savefig("imgs/g_lin.png")
plt.yscale('log')
plt.xlim(-0.01, time[-1])
#plt.savefig("imgs/g_log.png")
plt.show()
plt.close()
res = mc.mc(0, 10, d)
time = np.arange(0,len(res))/1e3
plt.plot(time, res/np.sum(res[100:]), label = "g = 0 mu_s = 10")
res = mc.mc(0.9, 100, d)
time = np.arange(0,len(res))/1e3
plt.plot(time, res/np.sum(res[100:]), label = "g = 0.9 mu_s = 100")
plt.xlabel("time [ns]")
plt.ylabel("amplitude")
plt.legend()
plt.xlim(-0.01, 0.15)
#plt.savefig("imgs/comp_lin.png")
plt.yscale('log')
plt.xlim(-0.01, time[-1])
#plt.savefig("imgs/comp_log.png")
plt.show()
plt.close()

