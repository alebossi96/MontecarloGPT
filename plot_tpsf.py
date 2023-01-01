import matplotlib.pyplot as plt
import numpy as np
import math
values = []
start = 1000
with open("TPSF_3.txt", "r") as f:
  for line in f:
    values.append(int(line.strip()))
values = np.array(values)
time = np.arange(0,len(values))
plt.plot(time,values/np.sum(values[start:]))

def de_0(time):
    return 1/time**(3/2)
plt.plot(time,de_0(time)/np.sum(de_0(time)[start:]))
plt.yscale('log')
plt.show()


