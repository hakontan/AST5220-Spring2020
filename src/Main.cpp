#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.25;
  double OmegaLambda = 0.7;
  double Neff        = 3.046;
  double TCMB        = 2.7255;

  // Recombination parameters
  double Yp          = 0.0;
  
  //=========================================================================
  // Module I
  //=========================================================================


  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaLambda, Neff, TCMB);
  cosmo.solve();
  //std::cout << "solved" << std::endl;
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");


  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed


  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue1 = 0.0001 / Constants.Mpc;
  double kvalue2 = 0.001 / Constants.Mpc;
  double kvalue3 = 0.05 / Constants.Mpc;
  double kvalue4 = 0.1 / Constants.Mpc;
  double kvalue5 = 0.3 / Constants.Mpc;
  pert.output(kvalue1, "perturbations_k0.0001.txt");
  pert.output(kvalue2, "perturbations_k0.001.txt");
  pert.output(kvalue3, "perturbations_k0.05.txt");
  pert.output(kvalue4, "perturbations_k0.1.txt");
  pert.output(kvalue5, "perturbations_k0.3.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
