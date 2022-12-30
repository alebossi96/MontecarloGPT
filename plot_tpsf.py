import matplotlib.pyplot as plt
import numpy as np
import math
values = []
with open("TPSF.txt", "r") as f:
  for line in f:
    values.append(int(line.strip()))
values = np.array(values)
plt.plot(values)
plt.show()
