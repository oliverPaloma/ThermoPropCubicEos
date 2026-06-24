#include "../../ThermoPropCubicEos.hpp"

int main ()
{
  auto T = 323.15;
  auto P = 30000000.;
  auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
  auto components = "CO2 H2O";
  std::vector<real> Tc,Pc,omega;
  std::vector<std::vector<real>> kij(2, std::vector<real>(2)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  auto EoSModel = CubicEOSModel::PengRobinson;
  std::vector<real> z{0.5, 0.5};

  read_database(Tc,Pc,omega,databasePath,components);

  CubicEOSProps props;

  compute(props,Tc,Pc,omega,T,P,z,EoSModel,kij);
  //props.lnphixk[i][k] = 10;


  auto a = 111;

  return 0;
}
