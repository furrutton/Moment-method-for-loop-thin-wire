import numpy as np
from matplotlib import pyplot as plt

y1, y2 = np.loadtxt("rb.txt", unpack=True, usecols=(0,1), dtype=np.float)
plt.plot(y1, 'r-')
plt.plot(y2, 'b-')
plt.show()
