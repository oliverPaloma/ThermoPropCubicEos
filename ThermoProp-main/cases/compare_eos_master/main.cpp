#include "../../ThermoPropCubicEos.hpp"

int main ()
{

  // Critical parameters of methane
  std::vector<real> pc{4.e6};      // Critical pressure [Pa]
  std::vector<real> tc{190.6};    // Critical temperature [K]
  std::vector<real> omega{0.008}; // Acentric factor
  std::vector<std::vector<real>> kij(1, std::vector<real>(1)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  auto EoSModel = CubicEOSModel::PengRobinson;
  CubicEOSProps props;
  std::vector<real> z{1.0};  

  real p1 = 3.e7;    // Pressure [Pa]
  real t1 = 180.0;  // Temperature [K]
  compute(props,tc,pc,omega,t1,p1,z,EoSModel,kij);

  auto z1 = p1 * props.V / R / t1;
  auto phi1 = std::exp((double) props.ln_phi[0]); 

  std::cout << std::setprecision(16) << "z1 = " << z1 << std::endl;
  std::cout << "phi1 = " << phi1 << std::endl;

  /*
  using these conditions eos master gives:
  z1 = 0.8990446822387802
  phi1 = 0.16774802530911384

  these code gives (jun 17 2024)
  z1 = 0.8990768566592848
  phi1 = 0.1678533644764063
  */

  return 0;
}



