#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  //mistura 1
  //auto T = 350.;
  //auto T = 400.;
  //auto T = 450.; 
  
  //mistura 2
  //auto T = 270.;
  //auto T = 320.;
  //auto T = 380.;  /t / \ 

  //auto T = 310.;
  //auto T = 304.;
  auto T = 280.;

  auto databasePath = "/home/paloma/Documentos/ThermoPropCubicEos/database/test.yml"; //"/home/palomajo/Documentos/ThermoPropCubicEos/database/test.yml";"/home/fellipe/Dropbox/Pessoais/Downloads/ThermoPropCubicEos/database/test.yml";
  //auto components = "CO2 H2O"; //apenas um componete //testar com ('NC6,NC10') # Hexane/decane 
  auto components = "CO2"; 
  std::vector<double> z{1.}; //fração molar

  std::vector<double> Tc,Pc,omega;
  //MIX1 (Hamidreza et al., 2013):
  //auto components = "C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 ";
  //std::vector<double> z{0.413506, 0.040300, 0.215300, 0.053900, 0.054300, 0.051500, 0.051900, 0.103900, 0.003470, 0.002680, 0.002070, 0.001590, 0.001230, 0.000950, 0.000730, 0.000566, 0.000437, 0.001671};
  //auto T = 200.; //antes do ponto crítico
  //auto T = 398.5; //no ponto crítico //Temperatura crítica: entre 398,03 K e 399,08 K auto T = 398.5; 
  //auto T = 450.; //depois do ponto crítico
                 //Pressão crítica: entre 10,97 MPa e 10,99 MPa
                  

 //MIX2 (Pedersen et al., 2014):
 //auto components = "N2 CO2 C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20";
 //std::vector<double> z{0.001200, 0.024900, 0.764300, 0.074600, 0.031200, 0.005900, 0.012100, 0.005000, 0.005900, 0.007900, 0.009500, 0.010800, 0.007800, 0.005920, 0.004670, 0.003450, 0.003750, 0.003040, 0.002370, 0.002080, 0.002200, 0.001690, 0.001400, 0.008330};
 //auto T = 200.; //antes do ponto crítico
 //auto T = T = 322.5; //no ponto crítico
 //auto T = 400.; //depois do ponto crítico
                  //Temperatura crítica: entre 321,68 K e 325,70 K auto T = 322.5; 
                  //Pressão crítica: entre 36,78 MPa e 37,12 MPa .

  //auto components = "NC6 NC10";
  //std::vector<double> Tc,Pc,omega;
 
  auto EoSModel = CubicEOSModel::VanDerWaals; // PengRobinson   SoaveRedlichKwong  VanDerWaals
  
  //std::vector<double> z{1.}; //fração molar
  
  read_database(Tc,Pc,omega,databasePath,components);
  
  auto ncomp=z.size();

  auto Vi = 0.05; // m³/mol
  auto Vf = 0.45;
  auto npoints = 100;

 
  //eu tinha rodado sem a regra de mistura, agr coloquei aqui e vou criar a função pra calcular o Vi e vf ideal.
  //calculateIsoterma(EoSModel,Tc,Pc,omega,T,Vi,Vf,npoints); //CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)

  //CalcularVolumeIdeal(EoSModel,Tc,Pc, omega, z, ncomp, Vi, Vf);

  std::cout << "Vi=" << Vi << "   Vf="<< Vf << std::endl; 


  calculateIsotermaMisture(EoSModel,Tc,Pc,omega,T,Vi,Vf,npoints,z,ncomp);

  return 0;
}
