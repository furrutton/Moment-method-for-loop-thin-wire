import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('sigma.txt', usecols=(0,1), unpack=True, dtype=np.float)
plt.plot(x, y, 'b-',linewidth=0.7)
plt.show()
