#include "../../ThermoPropCubicEos.hpp"


int main ()
{
   auto T = 350.; //isoterma temperatura constante 
   //auto databasePath = "/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml";  //auto databasePath = "/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";
   //comentado para testar um componente apenas: 

auto databasePath = "/home/palomajo/Downloads/ThermoPropCubicEos-main/ThermoPropCubicEos-main/ThermoPropCubicEos/database/test.yml";

   
   //auto components = "CO2 C1";  //depois testar: //auto components = "NC6 NC10";  //auto components = "CO2 H2O"; 

//testar apenas um componente pra ver se bate com a isotermas
  auto components = "CO2";






  

   std::vector<double> Tc,Pc,omega;
  
   //auto EoSModel = CubicEOSModel::PengRobinson; 
   auto EoSModel = CubicEOSModel::VanDerWaals; 
   //auto EoSModel = CubicEOSModel::SoaveRedlichKwong;  

   //std::vector<double> z{0.2,0.8}; //fração molar
   std::vector<double> z{1.}; //fração molar

   read_database(Tc,Pc,omega,databasePath,components);
   auto ncomp=z.size();

   auto Vi = 0.05/1000.; // m³/mol
   auto Vf = 0.45/1000.;
   auto npoints = 100;

   calculateIsotermaMisture(EoSModel,Tc,Pc,omega,T,Vi,Vf,npoints,z,ncomp); //CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)
   return 0;
}

