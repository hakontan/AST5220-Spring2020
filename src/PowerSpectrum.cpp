#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  Vector k_array;
  Vector log_k_array = log(k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  double end_pts = 5000.0;
  int n_pts = 1000;
  Vector x_array = Utils::linspace(0.0, end_pts, n_pts);
  Vector j_ell_arr(n_pts);
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for(int j = 0; j < n_pts; j++) {
      j_ell_arr[j] = Utils::j_ell(ell, x_array[j]);

    j_ell_splines[i].create(x_array, j_ell_arr);
    // Make the j_ell_splines[i] spline
    }
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
  
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  //Vector2D result(ells.size(), Vector(k_array.size()));

  double dx;
  double eta;
  double eta_nm1;
  double f_n;
  double f_nm1;
  double Integral = 0.0;
  double eta_0 = cosmo->eta_of_x(0);
  int npts;
  Vector x_array;

  for (int ell = 0; ell < ells.size(); ell ++) {
    for(int ik = 0; ik < k_array.size(); ik++) {
      // Number of points for the x-grid (Callin (2006))
      npts = (int) std::round(10 * ik * eta_0 / (2 * M_PI));
      x_array  = Utils::linspace(Constants.x_start, Constants.x_end, npts);
      for(int j = 1; j < npts; j++) {
        // Assigning timestep for given k (Callin (2006))
        dx = 2 * M_PI * cosmo->Hp_of_x(x_array[j]) / (10.0 * ik);
        // Assigning variables for trapezoidal rule
        eta     = cosmo->eta_of_x(x_array[j]);
        eta_nm1 = cosmo->eta_of_x(x_array[j-1]);
        f_n     = source_function(x_array[j], ik) * j_ell_splines[ell](ik*(eta - eta_0));
        f_nm1   = source_function(x_array[j-1], ik) * j_ell_splines[ell](ik*(eta_nm1 - eta_0));
        // Trapezoidal rule
        Integral += 0.5 * (f_nm1 + f_n) * dx;
      }
      result[ell, ik] = Integral;
      Integral = 0.0;
    }
  }
  

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    // ...
    // ...
    // ...
    
    // Store the result for Source_ell(k) in results[ell][ik]
  

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // ...
    // ...
    // ...
    // ...

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================

  // ...
  // ...
  // ...
  // ...

  Vector result;

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================

  // ...
  // ...
  // ...

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

