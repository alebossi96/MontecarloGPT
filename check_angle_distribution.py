import matplotlib.pyplot as plt
import numpy as np
import math
class obj_det:
    def __init__(self, x,y,z, r):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
def p(cos_theta, g):
    return 1/(4*math.pi)*(1-g**2)/(1+g**2-2*g*cos_theta)**(1.5)
d = obj_det(0,0,0,0.1)

import montecarlomodule as mc
time = np.arange(0,1e3)/1e3
for g in [0.9]:
    ang_dist = mc.angle_distribution(g)
    angle = []
    for _ in range(100000):
        idx = np.random.randint(0,len(ang_dist))
        angle.append(ang_dist[idx])
    bins = 100
    values = np.array(angle)
    hist, bin_edges = np.histogram(values, bins=bins, density=True)

    #values = values[values<10]
    plt.plot(np.cos(bin_edges[:bins]), hist/np.max(hist), label="from_dists")
    step = 1e-3
    x = np.arange(step,math.pi,step = step)
    plt.plot(np.cos(x),p(np.cos(x), g)*np.sin(x)*2*math.pi/np.max(p(np.cos(x), g)*np.sin(x)*2*math.pi), 'k--', label = "Theoretical")
    print(np.sum(p(np.cos(x),g)*np.sin(x))*step*2*math.pi)
    print(np.sum(hist)*(bin_edges[1]- bin_edges[0]))
    values = mc.test_angle (int(1e4),g);
    bins = 1000
    hist, bin_edges = np.histogram(values, bins=bins, density=True)
    #values = values[values<10]
    plt.plot(bin_edges[:bins], hist/np.max(hist), label="Measured")
    plt.show()  
    """
    res = mc.mc(g,20, d)
    values = []
    
    with open("cos_ang.txt", "r") as f:
      for line in f:
        values.append(float(line.strip()))
    values = np.array(values)
    bins = 1000
    hist, bin_edges = np.histogram(values, bins=bins, density=True)

    #values = values[values<10]
    plt.plot(bin_edges[:bins], hist/np.max(hist), label="Measured")
    plt.show()  
    """
    """
    step = 1e-3
    
    x = np.arange(step,math.pi,step = step)
    plt.plot(np.cos(x),p(np.cos(x), g)*np.sin(x)*2*math.pi/np.max(p(np.cos(x), g)*np.sin(x)*2*math.pi), 'k--', label = "Theoretical")
    print(np.sum(p(np.cos(x),g)*np.sin(x))*step*2*math.pi)
    print(np.sum(hist)*(bin_edges[1]- bin_edges[0]))
    plt.show()
    """



