import numpy as np
from matplotlib import pyplot as plt

x = np.loadtxt("An_dB.txt", unpack=True, usecols=(0), dtype=np.float)
plt.plot(x)
plt.show()
