import numpy as np
import matplotlib.pyplot as plt
import astropy.units as unit
import astropy.constants as const

#plt.style.use("seaborn-darkgrid")

fonts = {
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}

plt.rcParams.update(fonts)

data = np.loadtxt("recombination.txt")

x = data[:, 0]
Xe_of_x = data[:, 1]
ne_of_x = data[:, 2]
tau_of_x = data[:, 3]
dtaudx_of_x = data[:, 4]
ddtauddx_of_x = data[:, 5]

plt.semilogy(x, Xe_of_x)
plt.figure()
plt.semilogy(x, ne_of_x)
plt.figure()
plt.semilogy(x, tau_of_x)
plt.semilogy(x, dtaudx_of_x)
plt.semilogy(x, ddtauddx_of_x)
plt.show()