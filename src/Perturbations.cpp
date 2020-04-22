#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  //compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  Vector k_array(n_k);
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, n_x);
  Vector Psi(n_k * n_x);
  Vector Phi(n_k * n_x);
  Vector v_cdm(n_k * n_x);
  Vector delta_cdm(n_k * n_x);
  Vector v_b(n_k * n_x);
  Vector delta_b(n_k * n_x);
  Vector2D Theta(Constants.n_ell_theta, Vector(n_k * n_x));
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);

  //Initializing k_array with logarithmic spacing
  double delta_k = std::log(k_max / k_min) / double(n_k - 1);
  for(int i = 0; i < n_k; i++){
      k_array[i] = std::log(k_min) + delta_k * i;
  }
  k_array[n_k - 1] = std::log(k_max);

  for(int i = 0; i < n_k; i++){
      k_array[i] = std::exp(k_array[i]);
  }

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    std::pair<double,int> x_end_data = get_tight_coupling_time(k);
    const double x_end_tight = x_end_data.first;
    const int x_end_index = x_end_data.second;
    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================
  
    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    Vector x_tc = Utils::linspace(x_start, x_end_tight, x_end_index);
    ODESolver tc_ODE;
    //tc_ODE.set_verbose(true);

    tc_ODE.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
  
    auto Phi_tc         = tc_ODE.get_data_by_component(Constants.ind_Phi_tc);
    auto delta_b_tc     = tc_ODE.get_data_by_component(Constants.ind_deltab_tc);
    auto delta_cdm_tc   = tc_ODE.get_data_by_component(Constants.ind_deltacdm_tc);
    auto v_cdm_tc       = tc_ODE.get_data_by_component(Constants.ind_vcdm_tc);
    auto v_b_tc         = tc_ODE.get_data_by_component(Constants.ind_vb_tc);
    auto theta0_tc       = tc_ODE.get_data_by_component(Constants.ind_start_theta_tc);
    auto theta1_tc       = tc_ODE.get_data_by_component(Constants.ind_start_theta_tc + 1);
    



    for (int ix = 0; ix < x_end_index; ix++) {
        Theta[0][ix + n_x * ik]   = theta0_tc[ix];
        Theta[1][ix + n_x * ik]   = theta1_tc[ix];
        Phi[ix + n_x * ik]        = Phi_tc[ix];
        delta_b[ix + n_x * ik]    = delta_b_tc[ix];
        delta_cdm[ix + n_x * ik]  = delta_cdm_tc[ix];
        v_b[ix + n_x * ik]        = v_b_tc[ix];
        v_cdm[ix + n_x * ik]      = v_cdm_tc[ix];
        Theta[2][ix + n_x * ik]   = 20.0 * Constants.c * k * Theta[1][ix + n_x * ik]
                                    / (45.0 * cosmo->Hp_of_x(x_array[ix]) * rec->dtaudx_of_x(x_array[ix]));
    }
    
    for (int ell = 3; ell < Constants.n_ell_theta - 1; ell++) {
      for (int ix = 0; ix < x_end_index; ix++) {
        Theta[ell][ix + n_x * ik] = ell * Constants.c * k * Theta[ell-1][ix + n_x * ik]
                                    / ((2.0 * ell + 1.0) * cosmo->Hp_of_x(x_array[ix])
                                    * rec->dtaudx_of_x(x_array[ix]));
        
      }
    }
    
    
    for (int ix = 0; ix < x_end_index; ix++) {
      Psi[ix + n_x * ik]       = - Phi_tc[ix]
                                  - 12.0 * cosmo->get_H0() * cosmo->get_H0() * cosmo->get_OmegaR(0)
                                  * Theta[2][ix + n_x * ik]
                                  / (Constants.c * Constants.c * k * k * std::exp(x_array[ix]) * std::exp(x_array[ix]));
    }
  
    //std::cout << "tc_ode solution stored" << std::endl;
    //===================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(tc_ODE.get_data_by_xindex(x_end_index-1), x_end_tight, k);
    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_full = Utils::linspace(x_end_tight, x_end, n_x - x_end_index);
    ODESolver full_ODE;
    //full_ODE.set_verbose(true);
    full_ODE.solve(dydx_full, x_full, y_full_ini);
    
    auto y_full = full_ODE.get_data(); 

    auto Phi_after         = full_ODE.get_data_by_component(Constants.ind_Phi);
    auto delta_b_after     = full_ODE.get_data_by_component(Constants.ind_deltab);
    auto delta_cdm_after   = full_ODE.get_data_by_component(Constants.ind_deltacdm);
    auto v_cdm_after       = full_ODE.get_data_by_component(Constants.ind_vcdm);
    auto v_b_after         = full_ODE.get_data_by_component(Constants.ind_vb);
    auto theta_after       = full_ODE.get_data_by_component(Constants.ind_start_theta);
    
      
    for (int ix = 0; ix < n_x - x_end_index; ix ++) {
        //std::cout << ix << std::endl;
        Phi[x_end_index + ix + n_x * ik]        = Phi_after[ix];
        delta_b[x_end_index + ix + n_x * ik]    = delta_b_after[ix];
        delta_cdm[x_end_index + ix + n_x * ik]  = delta_cdm_after[ix];
        v_b[x_end_index + ix + n_x * ik]        = v_b_after[ix];
        v_cdm[x_end_index + ix + n_x * ik]      = v_cdm_after[ix];
    }

    for (int ell = 0; ell < Constants.n_ell_theta; ell++) {
      auto theta_after       = full_ODE.get_data_by_component(Constants.ind_start_theta + ell);
      for (int ix = 0; ix < n_x - x_end_index; ix ++) {
        Theta[ell][x_end_index + ix + n_x * ik] = theta_after[ix];
        //std::cout << Theta[ell][x_end_index + ix + n_x * ik] << std::endl;
      }
    }

    //std::cout << "Full ODE solution stored" << std::endl;
    
    for (int ix = 0; ix < n_x - x_end_index; ix ++) {
       Psi[x_end_index + ix + n_x * ik] = - Phi_after[ix] 
                                            - 12.0 * cosmo->get_H0() * cosmo->get_H0() * cosmo->get_OmegaR(0)
                                            * Theta[2][x_end_index + ix + n_x * ik] 
                                            / (Constants.c * Constants.c * k * k 
                                            * std::exp(2 * x_array[ix + x_end_index]));

    }
    
    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    //===================================================================
    //...
    //...

  }
  Utils::StartTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
  Psi_spline.create(x_array, k_array, Psi, "Psi_spline");
  Phi_spline.create(x_array, k_array, Phi, "Phi_spline");
  delta_cdm_spline.create(x_array, k_array, delta_cdm, "delta_cdm_spline");
  delta_b_spline.create(x_array, k_array, delta_b, "delta_b_spline");
  v_cdm_spline.create(x_array, k_array, v_cdm, "v_cdm_spline");
  v_b_spline.create(x_array, k_array, v_b, "v_b_spline");


  for (int ell = 0; ell < Constants.n_ell_theta; ell++) {
    Theta_spline[ell].create(x_array, k_array, Theta[ell]);
  }
  
  std::cout << "Splines created" << std::endl;
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  //double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  double Psi_init       = - 2.0 / 3.0;
  double Phi_init       = - Psi_init;
  double delta_cdm_init = - Psi_init * 3.0 / 2.0;
  double v_cdm_init     = - Psi_init * Constants.c * k / (2.0 * cosmo->Hp_of_x(x));
  double delta_b_init   = delta_cdm_init;
  double v_b_init       = v_cdm_init;

  Phi                   = Phi_init;
  delta_cdm             = delta_cdm_init;
  v_cdm                 = v_cdm_init;
  delta_b               = delta_b_init;
  v_b                   = v_b_init;
  
  // SET: Photon temperature perturbations (Theta_ell)

  double theta_0_init   = - Psi_init / 2.0;
  double theta_1_init   = Psi_init * Constants.c * k / (6.0 * cosmo->Hp_of_x(x));


  Theta[0] = theta_0_init;
  Theta[1] = theta_1_init;

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  //Constants and variables used in calculation
  const double c              = Constants.c;
  const double Hp             = cosmo->Hp_of_x(x);
  const double dHpdx          = cosmo->dHpdx_of_x(x);
  const double H0             = cosmo->get_H0();
  const double Omega_CDM      = cosmo->get_OmegaCDM(0);
  const double Omega_B        = cosmo->get_OmegaB(0);
  const double Omega_R        = cosmo->get_OmegaR(0);
  const double a              = std::exp(x);
  const double R              = 4.0 * Omega_R / (3.0 * Omega_B * a);
  const double dtaudx         = rec->dtaudx_of_x(x);  
  const double ddtauddx       = rec->ddtauddx_of_x(x);

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  //const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  //double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  //double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)

  double Phi_init       = Phi_tc;
  double delta_cdm_init = delta_cdm_tc;
  double v_cdm_init     = v_cdm_tc;
  double delta_b_init   = delta_b_tc;
  double v_b_init       = v_cdm_tc;

  Phi                   = Phi_init;
  delta_cdm             = delta_cdm_init;
  v_cdm                 = v_cdm_init;
  delta_b               = delta_b_init;
  v_b                   = v_b_init;

  // SET: Photon temperature perturbations (Theta_ell)

  double theta_0_init   = Theta_tc[0];
  double theta_1_init   = Theta_tc[1];
  double theta_2_init   = - 20.0 * Constants.c * k * Theta[1]
                        / (45.0 * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x));
  double Psi       = - Phi - 12.0 * H0 * H0 * Omega_R * Theta[2] / (c * c * k * k * a * a);

  Theta[0] = theta_0_init;
  Theta[1] = theta_1_init;
  Theta[2] = theta_2_init;

  for (int ell = 3; ell < n_ell_theta; ell++) {
    Theta[ell] = - ell * Constants.c * k * Theta[ell-1] 
                 / ((2.0 * ell + 1.0) * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x));
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double,int> Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  double dtaudx;
  double ck_over_Hp;
  double Xe;
  Vector x_arr = Utils::linspace(x_start, x_end, n_x);
  int c = 0;
  for (int i = 0; i < n_x; i++) {
    dtaudx     = rec->dtaudx_of_x(x_arr[i]);
    ck_over_Hp = Constants.c * k / cosmo->Hp_of_x(x_arr[i]);
    Xe = rec->Xe_of_x(x_arr[i]);
    x_tight_coupling_end = x_arr[i];
    c++;
    if (std::fabs(dtaudx) < 10 * std::min(1.0, ck_over_Hp)) {   
      break;
    } 
    if (rec->Xe_of_x(x_arr[i]) < 0.99) {
      break;
    }
  }
  //std::cout << x_tight_coupling_end << std::endl;
  return std::pair<double,int>(x_tight_coupling_end, c);
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  Vector k_array;
  Vector x_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  //Constants used in the calculation
  const double c              = Constants.c;
  const double Hp             = cosmo->Hp_of_x(x);
  const double dHpdx          = cosmo->dHpdx_of_x(x);
  const double H0             = cosmo->get_H0();
  const double Omega_CDM      = cosmo->get_OmegaCDM(0);
  const double Omega_B        = cosmo->get_OmegaB(0);
  const double Omega_R        = cosmo->get_OmegaR(0);
  const double a              = std::exp(x);
  const double R              = 4.0 * Omega_R / (3.0 * Omega_B * a);
  const double dtaudx         = rec->dtaudx_of_x(x);  
  const double ddtauddx       = rec->ddtauddx_of_x(x);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  //double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Theta_2   = - 20.0 * Constants.c * k * Theta[1]
                     / (45.0 * cosmo->Hp_of_x(x) * rec->dtaudx_of_x(x));

  double Psi       = - Phi - 12.0 * H0 * H0 * Omega_R * Theta_2 / (c * c * k * k * a * a);
  
  dPhidx           = Psi - Phi * c * c * k * k / (3.0 * Hp * Hp) + H0 * H0 / (2.0 * Hp * Hp)
                     * (Omega_CDM * delta_cdm / a + Omega_B * delta_b / a + 4.0 * Omega_R * Theta[0] / (a * a));

  ddelta_cdmdx     = c * k * v_cdm / Hp - 3.0 * dPhidx;

  dv_cdmdx         = - v_cdm - c * k * Psi / Hp;

  ddelta_bdx       = c * k * v_b / Hp - 3.0 * dPhidx;


  // SET: Photon multipoles (Theta_ell)
  double q = (- ((1.0 - R)*dtaudx + (1.0 + R)*ddtauddx) * (3.0 * Theta[1] + v_b) - (c * k / Hp) * Psi
             + (1.0 - dHpdx / Hp) * c * k * (-Theta[0] + 2.0 * Theta_2) / Hp - (c * k / Hp) * dThetadx[0])
             / ((1.0 + R)*dtaudx + dHpdx / Hp - 1.0);

  dv_bdx = (1.0 / (1.0 + R)) * (- v_b - c * k * Psi / Hp + R * (q + (c * k / Hp) * (-Theta[0] + 2 * Theta_2) - c * k * Psi / Hp));

  dThetadx[0] = - c * k * Theta[1] / Hp - dPhidx;
  dThetadx[1] = (1.0 / 3) * (q - dv_bdx);  

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  //Constants and variables used in the calculation
  const double c              = Constants.c;
  const double Hp             = cosmo->Hp_of_x(x);
  const double dHpdx          = cosmo->dHpdx_of_x(x);
  const double H0             = cosmo->get_H0();
  const double Omega_CDM      = cosmo->get_OmegaCDM(0);
  const double Omega_B        = cosmo->get_OmegaB(0);
  const double Omega_R        = cosmo->get_OmegaR(0);
  const double a              = std::exp(x);
  const double R              = 4 * Omega_R / (3 * Omega_B * a);
  const double dtaudx         = rec->dtaudx_of_x(x);  
  const double ddtauddx       = rec->ddtauddx_of_x(x);
  const double ck_over_Hp     = c * k / Hp;
  const double eta            = cosmo->eta_of_x(x);
  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi = - Phi - 12.0 * H0 * H0 * Omega_R * Theta[2] / (c * c * k * k * a * a);
  

  dPhidx = Psi - Phi * c * c * k * k / (3.0 * Hp * Hp) + H0 * H0 / (2.0 * Hp * Hp)
           * (Omega_CDM * delta_cdm / a + Omega_B * delta_b / a + 4.0 * Omega_R * Theta[0] / (a * a));
  
  ddelta_cdmdx   = c * k * v_cdm / Hp - 3.0 * dPhidx;
  dv_cdmdx       = - v_cdm - c * k * Psi /Hp;

  ddelta_bdx     = c * k * v_b / Hp - 3 * dPhidx;
  dv_bdx         = - v_b - c * k * Psi / Hp + dtaudx * R * (3 * Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  dThetadx[0] = - c * k * Theta[1] / Hp - dPhidx;

  dThetadx[1] = ck_over_Hp * Theta[0] / 3.0 - ck_over_Hp * 2.0 * Theta[2] / 3.0
                + ck_over_Hp * Psi / 3.0 + dtaudx * (Theta[1] + v_b / 3.0);

  dThetadx[2] = ck_over_Hp * 2 * Theta[2-1] / (2 * 2 + 1) 
                  - ck_over_Hp * (2 + 1) * Theta[2+1] / (2 * 2 + 1)
                  + dtaudx * (Theta[2] - Theta[2] / 10.0);
  //std::cout << "Theta ell p1" << std::endl;

  int last_ell = Constants.n_ell_theta - 1;
  for (int ell = 3; ell < last_ell; ell++) {
    dThetadx[ell] = ck_over_Hp * ell * Theta[ell-1] / (2 * ell + 1) 
                    - ck_over_Hp * (ell + 1) * Theta[ell+1] / (2 * ell + 1)
                    + dtaudx * Theta[ell];
  }
  //std::cout << "Theta ell p2" << std::endl;
  dThetadx[last_ell] = ck_over_Hp * Theta[last_ell - 1] 
                            - c * (last_ell + 1) * Theta[last_ell] / (Hp * eta)
                            + dtaudx * Theta[last_ell];
  //std::cout << "Theta ell p3" << std::endl;

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    //fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Phi(x,k)        << " ";
    fp << get_delta_b(x,k)    << " ";
    fp << get_delta_cdm(x,k)  << " ";
    fp << get_v_b(x,k)        << " ";
    fp << get_v_cdm(x, k)     << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    //fp << get_Source_T(x,k)  << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

