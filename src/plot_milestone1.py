import numpy as np
import matplotlib.pyplot as plt
import astropy.units as unit
import astropy.constants as const

plt.style.use("seaborn-darkgrid")

fonts = {
    "font.family": "serif",
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
}

plt.rcParams.update(fonts)

conversion_factor = 3.08567758e22 / 1e3 # factor to convert s^-1 to km/s/Mpc

def x_to_z(x):
    """
    Returns redshift z as a function of x=log(a).

    Parameters:
    -------------
    x:
        x = log(a)
    """
    return(np.exp(-x) - 1)



data = np.loadtxt("cosmology.txt") 


x = data[:, 0]
eta = data[:, 1] * unit.m
eta = eta.decompose()
eta = eta.to(unit.Gpc)
Hp_of_x = data[:, 2] * conversion_factor
H_of_x = data[:, 3] * conversion_factor
Omega_B = data[:, 5]
Omega_CDM = data[:, 6]
Omega_Lambda = data[:, 7]
Omega_R = data[:, 8]
z = x_to_z(x)



fig, ax = plt.subplots(2, 2)

ax[0, 0].semilogy(x, eta.value)
ax[0, 0].set_title(r"$\eta(x)$")
ax[0, 0].set_xlabel(r"$x=log(a)$")
ax[0, 0].set_ylabel("Gpc")
ax[0, 0].grid()

ax[0, 1].semilogy(x, Hp_of_x)
ax[0, 1].set_title(r"$\mathcal{H}(x)$")
ax[0, 1].set_xlabel(r"$x=log(a)$")
ax[0, 1].set_ylabel("km/s/Mpc")
ax[0, 1].grid()

ax[1, 0].semilogy(x, H_of_x)
ax[1, 0].set_title(r"H(x)")
ax[1, 0].set_xlabel(r"$x=log(a)$")
ax[1, 0].set_ylabel("km/s/Mpc")
ax[1, 0].grid()

ax[1, 1].loglog(z[np.where(z>0)], H_of_x[np.where(z>0)])
ax[1, 1].set_title(r"H(z)")
ax[1, 1].set_xlabel(r"Redshift $z$")
ax[1, 1].set_ylabel("Km/s/Mpc")
ax[1, 1].grid()
plt.tight_layout()
fig.savefig("../doc/milestoneI/figures/evolution.pdf", dpi=1000)


fig, ax = plt.subplots()
ax.plot(x, Omega_B + Omega_CDM, label=r"$\Omega_b(x) + \Omega_{CDM}(x)$")
ax.plot(x, Omega_Lambda, label=r"$\Omega_\Lambda$(x)")
ax.plot(x, Omega_R, label=r"$\Omega_r(x)$")
ax.set_xlabel(r"x=log(a)")
ax.set_ylabel("Fraction of total energy content")
ax.legend()
ax.grid()
fig.savefig("../doc/milestoneI/figures/densityparams.pdf", dpi=1000)


plt.show()