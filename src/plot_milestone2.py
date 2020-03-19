import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import astropy.units as unit
import astropy.constants as const

def x_to_z(x):
    """
    Returns redshift z as a function of x=log(a).
    Parameters:
    -------------
    x:
        x = log(a)
    """
    return(np.exp(-x) - 1)


fonts = {
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}

plt.rcParams.update(fonts)

save_directory = "../doc/milestoneII/figures/"
data = np.loadtxt("recombination.txt")

x = data[:, 0]
Xe_of_x = data[:, 1]
ne_of_x = data[:, 2]
tau_of_x = data[:, 3]
dtaudx_of_x = data[:, 4]
ddtauddx_of_x = data[:, 5]
g_of_x = data[:, 6]
dgdx_of_x = data[:, 7]
ddgddx_of_x = data[:, 8]


x_lss = x[np.argmax(g_of_x)] # x=log(a) coordinate value of last scattering surface
x_rec = x[np.argmin(np.abs(Xe_of_x - 0.5))]
print(f"Time of last scattering surface x_* = {x_lss}")
print(f"Time of last scattering surface z_* = {x_to_z(x_lss)}")
print(f"Time of recombination x_rec = {x_rec}")
print(f"Time of recombination z_rec = {x_to_z(x_rec)}")


fig, ax = plt.subplots() 
ax.semilogy(x, Xe_of_x)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$X_e$")
ax.set_title("Free electron fraction")
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "xe.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.semilogy(x, ne_of_x)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$kg/m^3$")
ax.set_title("Free electron density")
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "ne.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.semilogy(x, tau_of_x, label=r"$\tau(x)$")
ax.semilogy(x, np.abs(dtaudx_of_x), label=r"$\vert\frac{d \tau(x)}{dx}\vert$")
ax.semilogy(x[10: -10], ddtauddx_of_x[10:-10], label=r"$\frac{d^2\tau(x)}{dx^2}$")
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_title("Optical depth and derivatives")
ax.legend()
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "tau.pdf", dpi=1000)

fig, ax = plt.subplots()
ax.plot(x, g_of_x)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$g(x)$")
ax.set_title("Visibility function")
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "g.pdf", dpi=1000)

fig, ax = plt.subplots()
ax.plot(x, dgdx_of_x)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\frac{dg(x)}{dx}$")
ax.set_title("visibility function derivative")
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "dgdx.pdf", dpi=1000)

fig, ax = plt.subplots()
ax.plot(x, ddgddx_of_x)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\frac{d^2g(x)}{dx^2}$")
ax.set_title("visibility function double derivative")
ax.grid()
plt.tight_layout()
fig.savefig(save_directory + "ddgdxx.pdf", dpi=1000)
#plt.show()

"""
Run example:
(master)> python plot_milestone2.py
Time of last scattering surface x_* = -6.98376
Time of last scattering surface z_* = 1077.9676679708941
Time of recombination x_rec = -7.16616
Time of recombination z_rec = 1293.8627707684411
"""