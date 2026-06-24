#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto CalcCritPoint = true;
 auto databasePath = "../../database/test.yml"; //"../../database/SNG1.yml"
 // auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
 
  //MIX2 (Pedersen et al., 2014):
    auto components = "N2 CO2 C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20";
    std::vector<real> z{0.001200, 0.024900, 0.764300, 0.074600, 0.031200, 0.005900, 0.012100, 0.005000, 0.005900, 0.007900, 0.009500, 0.010800, 0.007800, 0.005920, 0.004670, 0.003450, 0.003750, 0.003040, 0.002370, 0.002080, 0.002200, 0.001690, 0.001400, 0.008330};



  auto ncomp = z.size();
  param parameters;
  auto& [props, Tc, Pc, omega, model, kij] = parameters; 
  kij =  std::vector<std::vector<real>>(ncomp, std::vector<real>(ncomp)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  model = CubicEOSModel::PengRobinson;


  read_database(Tc,Pc,omega,databasePath,components);
  if (Tc.size() != ncomp || Pc.size() != ncomp || omega.size() != ncomp)
  {
    std::cout << "Failed to load component data from `" << databasePath << "`." << std::endl;
    return 1;
  }

  // Creating output data file  
  String outputFileName = "output.out";
  std::ofstream output_file (outputFileName);
  output_file << std::scientific << std::setprecision(10);
  if (!output_file.is_open()) 
  {
    std::cout << "Error opening output file: `" << outputFileName << "`." << std::endl; 
    return 0;
  }      
  // End Creating output data file 

  // Creating options for computations 
  IsocurveOpts opts; 
  EquilibriumType EquiType;

  auto LnPhiDerivType = LnPhiDerivativeType::Analytical;
  //auto LnPhiDerivType = LnPhiDerivativeType::Automatic;
  //auto LnPhiDerivType = LnPhiDerivativeType::Numerical;


      std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
 

  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.01,0.03,true,1.e6,true,true);   //std::cout << "Pontos calculados,iniciando bolha;" << std::endl;
  real P0,T0,beta;
  int idx; 
  EquiType = Bubble;
  idx = ncomp; // index for temperature among independent variables, index starts at 0

  //T0 = T0_bubble; // K
  //P0 = P0_bubble; // Pa

  T0 = 243; // K 189.7  245.7
  //260,620583410713
  P0 = 4.04e6; // Pa  1.4e7;
                                          real Temp0_Bubble = T0;  real Pres0_Bubble = P0;
 // computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts); // Begin Compute the bubble - point curve
     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
     std::chrono::duration<double> tempo = end - start;
    

     std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now(); // Begin Compute the dew - point curve
  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.1,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //Original

  // 9° parametro por padrão 0.1, alterado para 0.5 ////08
  createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.5,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //pular o ponto crítico
  
  EquiType = Dew;
  idx = ncomp+1; // index for pressure
  T0 = 243.; // K começa em 370 nos gráfcos... 
  P0 = 101325.; // 101325 Pa =  1 atm

  real Temp0_Dew = T0;  real Pres0_Dew = P0;

  computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts);  // End Compute the dew - point curve 
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> tempo1 = end1 - start1;

 
 std::cout << "Tempo de execução Bubble: " << tempo.count() << " segundos" << ",  T0 ="<< Temp0_Bubble << "  P0 ="<< Pres0_Bubble << std::endl;
  std::cout << "\n Coeficiente usado:     \n Tempo de execução Dew : " << tempo1.count() << " segundos" <<",  \n T0 ="<< Temp0_Dew << "\n  P0 =" << Pres0_Dew << std::endl; //[" << AtoString(LnPhiDerivType) << "]

  switch (LnPhiDerivType){
  case Analytical:std::cout << "Dados do tipo Analytical salvos em: output.out" << std::endl;break;
  case Automatic:  std::cout << "Dados salvos em: output_automatica.out" << std::endl;break;
  case Numerical:  std::cout << "Dados salvos em: output_numerica.out" << std::endl;break;default:
  break;}

  return 0;
}



