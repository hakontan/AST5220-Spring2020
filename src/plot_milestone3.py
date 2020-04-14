import numpy as np
import matplotlib.pyplot as plt

save_directory = "../doc/milestoneIII/figures/"
data_k0 = np.loadtxt("perturbations_k0.01.txt")

x         = data_k0[:, 0]
Psi       = data_k0[:, 1]
Phi       = data_k0[:, 2]
delta_b   = data_k0[:, 3]
delta_cdm = data_k0[:, 4]
v_b       = data_k0[:, 5]
v_cdm     = data_k0[:, 6]

fig, ax = plt.subplots() 
ax.plot(x, Psi)
ax.set_xlabel(r"$x=$log$(a)$")
ax.set_ylabel(r"$\Psi$")
ax.set_title(r"$\Psi$")
ax.grid()
#ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "Psi.pdf", dpi=1000)


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
ax.semilogy(x, np.abs(delta_b), "b--", label=r"$\delta_b$, k=0.01/Mpc")
ax.semilogy(x, np.abs(delta_cdm), "r--", label=r"$\delta_{cdm}$, k=0.01/Mpc")
ax.set_xlabel(r"$x=$log$(a)$")
ax.grid()
ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "delta.pdf", dpi=1000)



fig, ax = plt.subplots() 
ax.semilogy(x, np.abs(v_b), "b--", label=r"$v_b$, k=0.01")
ax.semilogy(x, np.abs(v_cdm), "r--", label=r"$v_{cdm}$, k=0.01")
ax.set_xlabel(r"$x=$log$(a)$")
ax.grid()
ax.legend()
plt.tight_layout()
fig.savefig(save_directory + "v.pdf", dpi=1000)

