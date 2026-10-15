import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("wyniki/minimalizacja.dat")

Epot = data[:, 2]
it = data[:, 0]

fig, ax = plt.subplots(figsize=(6,4))
ax.plot(it, Epot)
plt.show()
