import numpy as np
import matplotlib.pyplot as plt
import re

class PlotPerturbations():
        def __init__(self, filenames):
            self.filenames=filenames
            self.save_directory = "../doc/milestoneIII/figures/"

        def store_data(self, filenames):
            data_k0 = np.loadtxt(filenames)

            self.x         = data_k0[:, 0]
            self.Psi       = data_k0[:, 1]
            self.Phi       = data_k0[:, 2]
            self.delta_b   = data_k0[:, 3]
            self.delta_cdm = data_k0[:, 4]
            self.v_b       = data_k0[:, 5]
            self.v_cdm     = data_k0[:, 6]
            self.theta0    = data_k0[:, 7]
            self.theta1    = data_k0[:, 8]
            self.theta2    = data_k0[:, 9]

        def plot(self):
            figPsi, axPsi     = plt.subplots() 
            figPhi, axPhi     = plt.subplots()
            figdelta, axdelta = plt.subplots()
            figv, axv         = plt.subplots()
            fig0, ax0         = plt.subplots()
            fig1, ax1         = plt.subplots()
            fig2, ax2         = plt.subplots()
            for files in self.filenames:
                self.store_data(files)
                
                kvalue = re.findall("(\d+)", files)      

                axPsi.plot(self.x, self.Psi, label=f"k={kvalue[0]}.{kvalue[1]}/Mpc")
                axPsi.set_xlabel(r"$x=$log$(a)$")
                axPsi.set_ylabel(r"$\Psi$")
                axPsi.set_title(r"$\Psi$")
                axPsi.grid()
                axPsi.legend()
                             
                axPhi.plot(self.x, self.Phi, label=f"k={kvalue[0]}.{kvalue[1]}/Mpc")
                axPhi.set_xlabel(r"$x=$log$(a)$")
                axPhi.set_ylabel(r"$\Phi$")
            
                axPhi.grid()
                axPhi.legend()
                       
                axdelta.semilogy(self.x, np.abs(self.delta_b), label=fr"$\delta_b$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                axdelta.semilogy(self.x, np.abs(self.delta_cdm), "--", label=fr"$\delta_{{cdm}}$, k={kvalue}/Mpc")
                axdelta.set_xlabel(r"$x=$log$(a)$")
                axdelta.grid()
                axdelta.legend()
                
                axv.semilogy(self.x, np.abs(self.v_b), label=fr"$v_b$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                axv.semilogy(self.x, np.abs(self.v_cdm), "--", label=fr"$v_{{cdm}}$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                axv.set_xlabel(r"$x=$log$(a)$")
                axv.grid()
                axv.legend()

                ax0.plot(self.x, self.theta0, label=fr"$\theta_0$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                ax0.set_xlabel(r"$x=$log$(a)$")
                ax0.grid()
                ax0.legend()

                ax1.plot(self.x, self.theta2, label=fr"$\theta_1$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                ax1.set_xlabel(r"$x=$log$(a)$")
                ax1.grid()
                ax1.legend()

                ax2.plot(self.x, self.theta2, label=fr"$\theta_2$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                ax2.set_xlabel(r"$x=$log$(a)$")
                ax2.grid()
                ax2.legend()
            

            axPhi.set_title(r"$\Phi$")
            axPsi.set_title(r"$\Psi$")
            ax0.set_title(r"$\theta_0$")
            ax1.set_title(r"$\theta_1$")
            ax1.set_xlabel(r"$x=$log$(a)$")
            ax1.grid()
            ax1.legend()
            ax2.set_title(r"$\theta_2$")

                
            figPsi.set_tight_layout(True)
            figPhi.set_tight_layout(True)
            figdelta.set_tight_layout(True)
            figv.set_tight_layout(True)
            figv.savefig(self.save_directory + "v.pdf", dpi=1000)
            figdelta.savefig(self.save_directory + "delta.pdf", dpi=1000)
            figPhi.savefig(self.save_directory + "Phi.pdf", dpi=1000)
            figPsi.savefig(self.save_directory + "Psi.pdf", dpi=1000)
            fig0.savefig(self.save_directory + "t0.pdf", dpi=1000)
            fig1.savefig(self.save_directory + "t1.pdf", dpi=1000)
            fig2.savefig(self.save_directory + "t2.pdf", dpi=1000)

if __name__=="__main__":
    filenames = ["perturbations_k0.001.txt", "perturbations_k0.01.txt", "perturbations_k0.1.txt"]
    plots = PlotPerturbations(filenames)
    plots.plot()