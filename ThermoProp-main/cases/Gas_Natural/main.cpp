#include "../../ThermoPropCubicEos.hpp"

int main ()
{
  auto CalcCritPoint = true;
  auto databasePath = "../../database/test.yml"; //"../../database/SNG1.yml"
 // auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
 
  auto components = "C1 C2 IC5";
  std::vector<real> z{0.94085, 0.04468, 0.01447};
   
  auto ncomp = z.size();
  param parameters;
  auto& [props, Tc, Pc, omega, model, kij] = parameters; 
  kij =  std::vector<std::vector<real>>(ncomp, std::vector<real>(ncomp)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  model = CubicEOSModel::PengRobinson;


  read_database(Tc,Pc,omega,databasePath,components);
  if (Tc.size() != ncomp || Pc.size() != ncomp || omega.size() != ncomp)
  { std::cout << "Failed to load component data from `" << databasePath << "`." << std::endl;
    return 1;}

  
  String outputFileName = "output.out"; // Creating output data file  
  std::ofstream output_file (outputFileName);
  output_file << std::scientific << std::setprecision(10);
  if (!output_file.is_open()) 
  { std::cout << "Error opening output file: `" << outputFileName << "`." << std::endl; 
    return 0;}      
  // End Creating output data file 

  // Creating options for computations 
  IsocurveOpts opts; 
  EquilibriumType EquiType;

  auto LnPhiDerivType = LnPhiDerivativeType::Analytical;
  //auto LnPhiDerivType = LnPhiDerivativeType::Automatic;
  //auto LnPhiDerivType = LnPhiDerivativeType::Numerical;

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();//inicio de contagem de tempo
 
  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.01,0.03,true,1.e6,true,true);   //std::cout << "Pontos calculados,iniciando bolha;" << std::endl;
  real P0,T0,beta, T0_Bubble,P0_Bubble, T0_Dew,P0_Dew;
  int idx; 
  EquiType = Bubble;
  idx = ncomp; // index for temperature among independent variables, index starts at 0
  T0 = 180.7; 
  P0 = 101325.; 
  T0_Bubble = T0; // K
  P0_Bubble = P0; // Pa

  //computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts); // Begin Compute the bubble - point curve
     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
     std::chrono::duration<double> tempo = end - start;

     std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now(); // Begin Compute the dew - point curve
  
  
 // createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.1,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //Original
  
  //9° parametro por padrão 0.1, alterado para 0.5 
  createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.89,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //pular o ponto crítico
  EquiType = Dew;
  idx = ncomp+1; // index for pressure 
  T0 = 230. ; 
  P0 = 101325.; 
  T0_Dew = T0;
  P0_Dew = P0;

  computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts);  // End Compute the dew - point curve 

  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> tempo1 = end1 - start1;

 
 std::cout << "Tempo de execução Bubble: " << tempo.count() << " segundos" << ",  T0 ="<< T0_Bubble << "  P0 ="<< P0_Bubble << std::endl;
  std::cout << "Tempo de execução Dew : " << tempo1.count() << " segundos" <<",  T0 ="<< T0_Dew << "  P0 =" << P0_Dew << std::endl; //[" << AtoString(LnPhiDerivType) << "]

  //std::cout << "Gerado por:" << String.LnPhiDerivType << std::endl; não funciona assim


 switch (LnPhiDerivType){
  case Analytical:std::cout << "Dados do tipo Analytical salvos em: output.out" << std::endl;break;
  case Automatic:  std::cout << "Dados salvos em: output_automatica.out" << std::endl;break;
  case Numerical:  std::cout << "Dados salvos em: output_numerica.out" << std::endl;break;default:
  break;}
 

  return 0;
}



