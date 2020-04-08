import numpy as np
import matplotlib.pyplot as plt

save_directory = "../doc/milestoneIII/figures/"
data_k0 = np.loadtxt("perturbations_k0.01.txt")

x         = data_k0[:, 0]
Phi       = data_k0[:, 1]
delta_b   = data_k0[:, 2]
delta_cdm = data_k0[:, 3]
v_b       = data_k0[:, 4]
v_cdm     = data_k0[:, 5]

fig, ax = plt.subplots() 
ax.plot(x, Phi)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\Phi$")
ax.set_title(r"$\Phi$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "Phi.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.plot(x, delta_b)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\delta_b$")
ax.set_title(r"$\delta_b$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "deltab.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.plot(x, delta_cdm)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\delta_{cdm}$")
ax.set_title(r"$\delta_{cdm}$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "deltacdm.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.plot(x, v_b, label="Combined solution")
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$v_b$")
ax.set_title(r"$v_b$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "vb.pdf", dpi=1000)

fig, ax = plt.subplots() 
ax.plot(x, v_cdm, label="Combined solution")
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$v_{cdm}$")
ax.set_title(r"$v_{cdm}$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "vcdm.pdf", dpi=1000)

