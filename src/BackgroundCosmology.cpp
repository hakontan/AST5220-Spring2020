#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaLambda,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaLambda(OmegaLambda),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaK, H0, ...
  //=============================================================================
  //...
  //...
  //...
  //...
  H0 = Constants.H0_over_h * h; // s^-1
  OmegaR = 8 * Constants.G * std::pow(M_PI, 3) * std::pow(Constants.k_b * TCMB, 4)
                  / (45 * H0 * H0 * std::pow(Constants.hbar, 3) * std::pow(Constants.c, 5));
  double OmegaNu = 0.0;
  double OmegaK = 0.0;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  //double x_start = 1e-7;
  //double x_end   = 1.0;
  int npts       = 100;

  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector eta(npts);
  
  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    //std::cout << "entered ODE function detadx" << std::endl;
    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...
    
    detadx[0] = Constants.c / Hp_of_x(x);
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...
  
  ODESolver ode;
  Vector eta_ic{0.0}; // initial conitions for eta
  ode.solve(detadx, x_array, eta_ic);
  auto all_data = ode.get_data();

  for (int i = 0; i<npts; i++) {
    eta[i] = all_data[i][0];
  }
  
  eta_of_x_spline.create(x_array, eta, "Test Spline");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================


double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double H = get_H0() * std::sqrt((OmegaB + OmegaCDM) * std::exp(-3 * x)
                                   + OmegaR * std::exp(-4 * x)
                                   + OmegaLambda);
  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double Hp_of_x = get_H0() * std::sqrt((OmegaB + OmegaCDM) * std::exp(-x)
                                   + OmegaR * std::exp(-2 * x)
                                   + OmegaLambda * std::exp(2 * x));

  return Hp_of_x;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double factor = (- (OmegaB + OmegaCDM) * std::exp(-x) 
                   - 2 * OmegaR * std::exp(-2*x)
                   + 2 * OmegaLambda * std::exp(2*x));

  return get_H0() * factor / (2 * std::sqrt(Hp_of_x(x) / get_H0()));
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double factor    = (- (OmegaB + OmegaCDM) * std::exp(-x) 
                      - 2 * OmegaR * std::exp(-2*x)
                      + 2 * OmegaLambda * std::exp(2*x));

  double ddxfactor = ((OmegaB + OmegaCDM) * std::exp(-x) 
                  - 4 * OmegaR * std::exp(-2*x)
                  + 4 * OmegaLambda * std::exp(2*x));

  return  0.5 * get_H0() * (ddxfactor * Hp_of_x(x) - dHpdx_of_x(x) * factor) / Hp_of_x(x);
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...


  return std::pow(get_H0() / H_of_x(x), 2) * std::exp(- 3 * x) * OmegaB;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return std::pow(get_H0() / H_of_x(x), 2) * std::exp(- 4 * x) * OmegaR;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return std::pow(get_H0() / H_of_x(x), 2) * std::exp(- 3 * x) * OmegaCDM;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return std::pow(get_H0() / H_of_x(x), 2) * OmegaLambda;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  //std::cout << "OmegaK:      " << OmegaK      << "\n";
  //std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = Constants.x_start;
  const double x_max =  Constants.x_end;
  const int    n_pts =  npts;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << H_of_x(x)          << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}



