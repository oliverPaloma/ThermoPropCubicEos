#include "../../ThermoPropCubicEos.hpp"


int main ()
{
  auto CalcCritPoint = true;
 auto databasePath = "../../database/test.yml"; //"../../database/SNG1.yml"
 // auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
 
 //OIL 1
    auto components = "CO2 C1 C2 C3 IC4 NC4 IC5 NC5 C6 McycloC5 Benzene CycloC6 McycloC6 Toluene C2Benzene mpXylene oXylene C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25 C26 C27 C28 C29 C30plus SumC7plus"; 
    std::vector<real> z{0.10172, 0.51249, 0.08113, 0.05354, 0.00010, 0.00003, 0.00008, 0.00007, 0.00129, 0.00084, 0.00005, 0.00088, 0.00255, 0.00022, 0.00054, 0.00057, 0.00032, 0.00437, 0.00940, 0.01206, 0.01425, 0.01314, 0.01267, 0.01360, 0.01227, 0.01221, 0.00978, 0.00940, 0.00946, 0.00840, 0.00720, 0.00673, 0.00620, 0.00561, 0.00547, 0.00522, 0.00471, 0.00462, 0.00447, 0.00444, 0.04792, 0.24360};


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
  idx = ncomp; // index for temperature among independent variables, index starts at 

  T0 = 150.; 
  P0 = 4.04e6; 
                                          real Temp0_Bubble = T0;  real Pres0_Bubble = P0;
 // computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts); // Begin Compute the bubble - point curve
     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
     std::chrono::duration<double> tempo = end - start;
    

     std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now(); // Begin Compute the dew - point curve
  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.1,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //Original

  // 9° parametro por padrão 0.1, alterado para 0.5 
  createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.89,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //pular o ponto crítico
  EquiType = Dew;
  idx = ncomp+1; 
  T0 = 1350.;  
  P0 = 101325.; 
  real Temp_Dew = T0;  real Pres_Dew = P0;                                                                                                     
  
  
  computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts);  // End Compute the dew - point curve 
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> tempo1 = end1 - start1;

 
 //std::cout << "Tempo de execução Bubble: " << tempo.count() << " segundos" << ",  T0 ="<< Temp0_Bubble << "  P0 ="<< Pres0_Bubble << std::endl;
   std::cout << "Tempo de execução Dew : " << tempo1.count() << " segundos" <<",  T0 ="<< Temp_Dew << "  P0 =" << Pres_Dew << std::endl; //[" << AtoString(LnPhiDerivType) << "]

  switch (LnPhiDerivType){
  case Analytical:std::cout << "Dados do tipo Analytical salvos em: output.out" << std::endl;break;
  case Automatic:  std::cout << "Dados salvos em: output_automatica.out" << std::endl;break;
  case Numerical:  std::cout << "Dados salvos em: output_numerica.out" << std::endl;break;default:
  break;}




  return 0;
}



