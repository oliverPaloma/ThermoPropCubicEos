#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto T = 200.;
  auto P = 99961.9;

  auto databasePath = "/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml";  // auto databasePath = "/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";
 
  auto components = "CO2 C1"; //Utilizando esses....
  std::vector<double> Tc,Pc,omega;
  std::vector<std::vector<double>> kij(2, std::vector<double>(2)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes

  auto EoSModel = CubicEOSModel::PengRobinson;

  //std::vector<double> z{0.5, 0.5}; //exemplo original

  std::vector<double> z{0.2, 0.8}; // Fração molar para mistura
  read_database(Tc,Pc,omega,databasePath,components);
  CubicEOSProps props;

  compute(props,Tc,Pc,omega,T,P,z,EoSModel,kij); //essa aqui

  //props.V;

  return 0;
}