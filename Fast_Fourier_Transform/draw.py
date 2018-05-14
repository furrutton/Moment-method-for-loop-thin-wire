import numpy as np
from matplotlib import pyplot as plt

y1 = np.loadtxt("An_dB100.csv", unpack=True, usecols=(0), dtype=np.float)
y2 = np.loadtxt("An_dB200.csv", unpack=True, usecols=(0), dtype=np.float)
y3 = np.loadtxt("An_dB300.csv", unpack=True, usecols=(0), dtype=np.float)
y4 = np.loadtxt("An_dB400.csv", unpack=True, usecols=(0), dtype=np.float)

x = [2*i/len(y1) for i in range(len(y1))]

plt.plot(x, y1, 'r-', x, y2, 'g-', x, y3, 'b-', x, y4, 'y-')
plt.title("Current distribution along the loop presented by dB")
plt.xlabel("rad")
plt.ylabel("dB")
plt.legend(("100 lambda", "200 lambda", "300 lambda", "400 lambda"))
plt.show()
