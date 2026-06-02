import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("wyniki_male/xrd.dat")

two_theta = data[:,0]
intensity = data[:,1]

# Normalize intensity
intensity /= np.max(intensity)

plt.figure(figsize=(10,5))

plt.plot(two_theta, intensity, linewidth=1.5)

plt.xlabel(r"$2\theta$ [deg]")
plt.ylabel("Normalized Intensity")
plt.title("Simulated XRD Pattern")
plt.grid(True, linestyle="--", alpha=0.6)

plt.xlim(two_theta.min(), two_theta.max())

plt.tight_layout()

plt.savefig("wyniki/xrd_pocz.png", dpi=300)

plt.show()