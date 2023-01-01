import matplotlib.pyplot as plt
import numpy as np
import math
values = []
start = 100
with open("TPSF_2.txt", "r") as f:
  for line in f:
    values.append(int(line.strip()))
values = np.array(values)
time = np.arange(0,len(values))
plt.plot(time[:200],values[:200]/np.sum(values[start:]))

def de_0(time):
    return 1/time**(3/2)
plt.plot(time[:200],de_0(time)[:200]/np.sum(de_0(time)[start:]))

plt.show()


