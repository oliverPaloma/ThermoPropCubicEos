#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto CalcCritPoint = true;
  auto databasePath = "../../database/test.yml"; //"../../database/SNG1.yml"
  //auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
 
 //OIL 2
    auto components = "O2 N2 CO2 C1 C2 C3 IC4 NC4 IC5 NC5 C6 McycloC5 Benzene CycloC6 McycloC6 Toluene C2Benzene mpXylene oXylene C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30plus C31 C32 C33 C34 C35 C36plus ";//SumC7plus
    std::vector<real> z{0.00005, 0.00054, 0.01241, 0.60587, 0.08553, 0.05512, 0.00104, 0.00135, 0.00254, 0.00217, 0.00539, 0.00236, 0.00016, 0.00194, 0.00410, 0.00061, 0.00052, 0.00081, 0.00059, 0.00748, 0.01136, 0.01209, 0.01321, 0.01153, 0.01071, 0.01157, 0.00964, 0.00973, 0.00771, 0.00739, 0.00733, 0.00651, 0.00578, 0.00511, 0.00481, 0.00457, 0.00422, 0.00397, 0.00438, 0.00367, 0.00389, 0.00368, 0.00411, 0.00333, 0.00343, 0.00314, 0.00278, 0.00257, 0.0272};//, 0.21691
    


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
  T0 = 1.9488e2; 
  P0 = 2.0265e5; 
  real Temp0_Bubble = T0;  real Pres0_Bubble = P0;

    

  std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now(); // Begin Compute the dew - point curve
  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.1,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //Original

  // 9° parametro por padrão 0.1, alterado para 0.5 
  createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.89,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //pular o ponto crítico
  EquiType = Dew;
  idx = ncomp+1; 
  T0 =  630.; 
  P0 =  101325; 
  real Temp0_Dew = T0;  real Pres0_Dew = P0;

  computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts);  // End Compute the dew - point curve 
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> tempo1 = end1 - start1;

 
 // std::cout << "Tempo de execução Bubble: " << tempo.count() << " segundos" << ",  T0 ="<< Temp0_Bubble << "  P0 ="<< Pres0_Bubble << std::endl;
  std::cout << "Tempo de execução Dew : " << tempo1.count() << " segundos" <<",  T0 ="<< Temp0_Dew << "  P0 =" << Pres0_Dew << "\n Coeficiente usado: \n  "<< std::endl; //[" << AtoString(LnPhiDerivType) << "]

 switch (LnPhiDerivType){
  case Analytical:std::cout << "Dados do tipo Analytical salvos em: output.out" << std::endl;break;
  case Automatic:  std::cout << "Dados salvos em: output_automatica.out" << std::endl;break;
  case Numerical:  std::cout << "Dados salvos em: output_numerica.out" << std::endl;break;default:
  break;}
  
  
  
  return 0;
}



