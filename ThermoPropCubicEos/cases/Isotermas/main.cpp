#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto T = 280.; //isoterma temperatura constante
 //auto P = 30000000.; obter
  
 // auto databasePath = "/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";

  auto databasePath = "/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml";
  //auto components = "CO2"; //apenas um componete
  auto components = "CO2 H2O"; //apenas um componete
  std::vector<double> Tc,Pc,omega;
 

  //auto EoSModel = CubicEOSModel::PengRobinson; //PR escolhido
  auto EoSModel = CubicEOSModel::VanDerWaals; // vdW escolhido
  //auto EoSModel = CubicEOSModel::SoaveRedlichKwong;  //srk escolhido

  //std::vector<double> z{1.}; //fração molar
  std::vector<double> z{0.5,0.5}; //fração molar para mistura

  read_database(Tc,Pc,omega,databasePath,components);
  

  auto Vi = 0.05/1000.; // m³/mol
  auto Vf = 0.45/1000.;
  auto npoints = 100;

  calculateIsoterma(EoSModel,Tc,Pc,omega,T,Vi,Vf,npoints,z); //CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)


  return 0;
}


