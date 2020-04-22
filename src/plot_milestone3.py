import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
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

        def plot_comparison(self, filename, nr_vert = 2):
            self.store_data(filename)
            kvalue = re.findall("(\d+)", filename)

            maxima0 = argrelextrema(self.theta0, np.greater)
            maxima1 = argrelextrema(self.theta1, np.greater)

            fig0, ax0 = plt.subplots(2)
            ax0[0].plot(self.x, self.theta0,
                     label=fr"$\theta_0$, k={kvalue[0]}.{kvalue[1]}/Mpc")
            ax0[1].semilogy(self.x, np.abs(self.delta_b),
                             label=fr"$\delta_b$, k={kvalue[0]}.{kvalue[1]}/Mpc")
            


            fig1, ax1 = plt.subplots(2)
            ax1[0].plot(self.x, self.theta1,
                     label=fr"$\theta_1$, k={kvalue[0]}.{kvalue[1]}/Mpc")
            ax1[1].semilogy(self.x, np.abs(self.v_b),
                             label=fr"$\delta_b$, k={kvalue[0]}.{kvalue[1]}/Mpc")

            print(maxima0)
            for i in range(nr_vert):
                ax0[0].axvline(x = self.x[maxima0[0][i]])
                ax0[1].axvline(x = self.x[maxima0[0][i]])
                ax1[0].axvline(x = self.x[maxima1[0][i]])
                ax1[1].axvline(x = self.x[maxima1[0][i]])

            fig0.savefig(self.save_directory + "comparison0.pdf", dpi=1000)
            fig1.savefig(self.save_directory + "comparison1.pdf", dpi=1000)
        


        def plot(self):
            a_eq = -8.68 # Output from milestone I
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
                axPsi.plot(self.x, self.Phi, label=f"k={kvalue[0]}.{kvalue[1]}/Mpc")
                       
                axdelta.semilogy(self.x, np.abs(self.delta_b), label=fr"$\vert\delta_b\vert$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                axdelta.semilogy(self.x, np.abs(self.delta_cdm), "--", label=fr"$\delta_{{cdm}}$, k={kvalue[0]}.{kvalue[1]}/Mpc")

                
                axv.semilogy(self.x, np.abs(self.v_b), label=fr"$\vert v_b \vert$, k={kvalue[0]}.{kvalue[1]}/Mpc")
                axv.semilogy(self.x, np.abs(self.v_cdm), "--", label=fr"$v_{{cdm}}$, k={kvalue[0]}.{kvalue[1]}/Mpc")

                ax0.plot(self.x, self.theta0, label=fr"$\theta_0$, k={kvalue[0]}.{kvalue[1]}/Mpc")

                ax1.plot(self.x, self.theta1, label=fr"$\theta_1$, k={kvalue[0]}.{kvalue[1]}/Mpc")
 
                ax2.plot(self.x, self.theta2, label=fr"$\theta_2$, k={kvalue[0]}.{kvalue[1]}/Mpc")
            

            axPhi.set_title(r"$\Phi$")
            axPsi.set_title(r"$\Psi$, $\Phi$")
            ax0.set_title(r"$\theta_0$")
            ax1.set_title(r"$\theta_1$")
            ax1.set_xlabel(r"$x=$log$(a)$")
            ax1.grid()
            ax1.legend()
            ax2.set_title(r"$\theta_2$")

            axPsi.set_xlabel(r"$x=$log$(a)$")
            axPsi.set_ylabel(r"$\Psi$ , $\Phi$")
    
            axPsi.grid()
            axPsi.legend()

            axPhi.set_xlabel(r"$x=$log$(a)$")
            axPhi.set_ylabel(r"$\Phi$")
            
            axPhi.grid()
            axPhi.legend()


            axdelta.set_xlabel(r"$x=$log$(a)$")
            axdelta.grid()
            axdelta.legend()


            ax0.set_xlabel(r"$x=$log$(a)$")
            ax0.grid()
            ax0.legend()


            ax1.set_xlabel(r"$x=$log$(a)$")
            ax1.grid()
            ax1.legend()

            ax2.set_xlabel(r"$x=$log$(a)$")
            ax2.grid()
            ax2.legend()

            axv.set_xlabel(r"$x=$log$(a)$")
            axv.grid()
            axv.legend()

            axPhi.axvline(x=a_eq)
            axPsi.axvline(x=a_eq)
            ax0.axvline(x=a_eq)
            ax1.axvline(x=a_eq)
            ax2.axvline(x=a_eq)
            axdelta.axvline(x=a_eq)
            axv.axvline(x=a_eq)


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
    filenames = ["perturbations_k0.001.txt",
                 "perturbations_k0.05.txt",
                 "perturbations_k0.3.txt"] 
                 #"perturbations_k0.1.txt",
                 #"perturbations_k0.3.txt"]
    plots = PlotPerturbations(filenames)
    plots.plot()

    plots.plot_comparison("perturbations_k0.1.txt",)