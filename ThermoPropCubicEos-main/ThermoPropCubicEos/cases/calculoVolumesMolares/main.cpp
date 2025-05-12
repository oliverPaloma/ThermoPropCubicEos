#include "../../ThermoPropCubicEos.hpp"

 
int main ()
{
  auto T = 350.;
  auto P = 3862600.;

  //auto databasePath = "/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml";  // auto databasePath = "/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";
  auto databasePath = "/home/paloma/Documentos/ThermoPropCubicEos/database/test.yml";

  auto components = "CO2"; 
  std::vector<double> Tc,Pc,omega;
  std::vector<std::vector<double>> kij(1, std::vector<double>(1)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes

  auto EoSModel = CubicEOSModel::PengRobinson;

  //std::vector<double> z{0.5, 0.5}; //exemplo original

  std::vector<double> z{1.}; // Fração molar para mistura
  read_database(Tc,Pc,omega,databasePath,components);
  CubicEOSProps props;
  double value;
  compute(props,Tc,Pc,omega,T,P,z,EoSModel,kij); //essa aqui

 
    std::cout << std::fixed << std::setprecision(16) << props.V << std::endl;



  return 0;
}