import numpy as np
from matplotlib import pyplot as plt

x, y, z = np.loadtxt("AVG_TIME.csv", unpack=True, usecols=(0,1,2), delimiter=',')
plt.plot(x, y, 'r--', x, z, 'b-*')
plt.xlabel("loop length")
plt.ylabel("time cost")
plt.legend(("length-time","length-step"))
plt.show()