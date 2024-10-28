#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto T = 280.; //isoterma temperatura constante
  auto databasePath = "/home/palomajo/Downloads/ThermoPropCubicEos/database/test.yml"; //leitura antiga

  auto components = "CO2"; //componentes escolhidos

  std::vector<double> Tc,Pc,omega;

  auto EoSModel = CubicEOSModel::PengRobinson; //pr escolhido 

  std::vector<double> z{1.}; //fração molar do componente na mistura 
  
  read_database(Tc,Pc,omega,databasePath,components); //chama
  
  auto Vi = 0.05/1000.; // m³/mol
  auto Vf = 0.45/1000.;
  auto npoints = 50;



  calculateIsoterma(EoSModel,Tc,Pc,omega,T,Vi,Vf,npoints); //CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)


  return 0;
}



