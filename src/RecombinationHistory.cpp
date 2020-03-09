#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector log_Xe_arr(npts_rec_arrays);
  Vector log_ne_arr(npts_rec_arrays);
  double Xe_current_peebles;
  double a;
  double n_b;
  const double OmegaB_0 = cosmo->get_OmegaB(0);
  bool finished = false;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      
      log_Xe_arr[i] = std::log(Xe_current);
      log_ne_arr[i] = std::log(ne_current);

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
        //dXedx[0] = rhs_peebles_ode(x, Xe, dXedx); //return rhs_peebles_ode(x, Xe, dXedx);
        //return GSL_SUCCESS;
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
      //Solving the Peebles ODE from the transition between Saha untill today
      Vector peeblesIC{Xe_current};
    
      Vector peebles_interval = Utils::linspace(x_array[i], x_array[npts_rec_arrays-1], npts_rec_arrays-i + 1);
      
      peebles_Xe_ode.solve(dXedx, peebles_interval, peeblesIC, gsl_odeiv2_step_rkf45);
      auto Xe_res = peebles_Xe_ode.get_data_by_component(0);
      double Xe_res_peebles;


      //updating Xe and ne arrays
      for (int k=0; k<npts_rec_arrays-i; k++) {
        a = std::exp(x_array[k+i]);
        n_b = OmegaB_0 * Constants.rhoc0 / (Constants.m_H * a * a * a);
        Xe_current_peebles = Xe_res[k];
        log_Xe_arr[k+i] = std::log(Xe_current_peebles);
        log_ne_arr[k+i] = std::log(n_b * Xe_current_peebles);
      }
      finished = true;
    }
    //Break the loop when Peebles equation is solved from onset untill today
    if(finished){
      break;
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "log_Xe_Spline");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "log_ne_Spline");


  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;
  const double rhoc0       = Constants.rhoc0;
  const double c           = Constants.c;
  const double OmegaB      = cosmo->get_OmegaB(0);
 


  double n_b = OmegaB * rhoc0 / (m_H * a * a * a);
  double F = std::pow(k_b * m_e * 2.725 / (2 * M_PI * a), 3.0/2) * std::exp(-epsilon_0 * a / (2.725 * k_b)) 
             / (n_b * hbar * hbar * hbar);
  // Electron fraction and number density
  double Xe;
  
  if (F > 1e5) {
    Xe = 1;
  }
  else {
    Xe = 0.5 * (-F + std::sqrt(F * F + 4 * F));
  }
  double ne = Xe * n_b;
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...
  
  return std::pair<double,double>(Xe, ne);
}



//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double rhoc0       = Constants.rhoc0;

  // calculating parameters for the RHS of the peebles equation.
  const double n_b = cosmo->get_OmegaB(0) * rhoc0 / (m_H * a * a * a);

  double n1s = (1 - X_e) * n_b;

  double T_b = 2.725 * k_b / a;

  double lambda_alpha = cosmo->H_of_x(x) * (27 * epsilon_0 * epsilon_0 * epsilon_0) 
                        / (64 * M_PI * M_PI * n1s * c * c * c * hbar * hbar * hbar);

  double phi2 = 0.448 * std::log(epsilon_0 / T_b);

  double alpha = std::sqrt(sigma_T * 3 / (8 * M_PI)) * m_e * c / hbar;
  double alpha2 = (64 * M_PI / std::sqrt(27 * M_PI)) * (alpha * alpha / (m_e * m_e))
                  * std::sqrt(epsilon_0 / T_b) * phi2
                  * hbar * hbar / c;

  double beta1 = alpha2 * std::pow(2.725 * m_e / (2 * M_PI * a), 3.0 / 2) * std::exp(-epsilon_0 / T_b)
                 * std::pow(k_b, 3.0/2) /(hbar * hbar * hbar);
  
  double factor_beta2 = 3 * epsilon_0 / (4 * T_b);
  double beta2;

  //Avoiding numerical error due to exponent becoming too large
  if (factor_beta2 < 200.0) {
    beta2 = beta1 * std::exp(factor_beta2);
  }
  else {
    beta2 = 0;
  }

  double Cr = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta2);
  
  //=============================================================================
  // Write the expression for dXedx
  //=============================================================================


  dXedx[0] = Cr * (beta1 * (1.0 - *Xe) - n_b * alpha2 * *Xe * *Xe) / cosmo->H_of_x(x);
  //std::cout << dXedx[0] << std::endl;
  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);
  Vector tau(npts);
  Vector g(npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // Set the derivative for photon optical depth
    //=============================================================================

    
    dtaudx[0] = - Constants.c * Constants.sigma_T * ne_of_x(x) / cosmo->H_of_x(x);
    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  
  ODESolver tau_ode;
  Vector tau_ic{1e5};
  tau_ode.solve(dtaudx, x_array, tau_ic, gsl_odeiv2_step_rkf45);

  auto tau_data = tau_ode.get_data_by_component(0);

  for (int i = 0; i<npts; i++) {
    tau[i] = tau_data[i] - tau_data[npts - 1];

  }
  
  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================

  for (int i = 0; i<npts; i++) {
    g[i] =  std::exp(-tau[i]) * Constants.c * Constants.sigma_T * ne_of_x(x_array[i]) / cosmo->H_of_x(x_array[i]);
  }

  tau_of_x_spline.create(x_array, tau, "tau Spline");
  g_tilde_of_x_spline.create(x_array, g, "g");
  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...
  return std::exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...
  return std::exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

