#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto T = 323.15;
  auto P = 30000000.;
 // auto databasePath = "/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";
  auto databasePath = "/home/paloma-j-oliveira/Documentos/ThermoPropCubicEos/database/test.yml";
  auto components = "CO2 H2O";
  std::vector<double> Tc,Pc,omega;
  std::vector<std::vector<double>> kij(2, std::vector<double>(2)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  auto EoSModel = CubicEOSModel::PengRobinson;
  std::vector<double> z{0.5, 0.5};

  read_database(Tc,Pc,omega,databasePath,components);

  CubicEOSProps props;

  compute(props,Tc,Pc,omega,T,P,z,EoSModel,kij);

  auto a = 111;

  return 0;
}



