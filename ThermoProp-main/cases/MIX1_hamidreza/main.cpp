#include "../../ThermoPropCubicEos.hpp"


int main ()
{
 auto CalcCritPoint = true;
 auto databasePath = "../../database/test.yml"; //"../../database/SNG1.yml"
  
 //MIX1 (Hamidreza et al., 2013):
    auto components = "C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16";
    std::vector<real> z{0.413506, 0.040300, 0.215300, 0.053900, 0.054300, 0.051500, 0.051900, 0.103900, 0.003470, 0.002680, 0.002070, 0.001590, 0.001230, 0.000950, 0.000730, 0.000566, 0.000437, 0.001671};
    std::vector<real> MM{16.043, 30.070, 44.097, 58.124, 58.124, 72.151, 72.151, 86.178, 100.205, 114.232, 128.259, 142.286, 156.313, 170.340, 184.367, 198.394, 212.421, 226.448};
    std::vector<int> flag{0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
 
  auto ncomp = z.size();
  param parameters;
  auto& [props, Tc, Pc, omega, model, kij] = parameters; 
  kij =  std::vector<std::vector<real>>(ncomp, std::vector<real>(ncomp)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  model = CubicEOSModel::PengRobinson;


  double zsum = 0.0;
  for (auto zi : z)
  zsum += val(zi);
  std::cout << "Soma(z) = " << zsum << std::endl;

  read_database(Tc,Pc,omega,databasePath,components);
  if (Tc.size() != ncomp || Pc.size() != ncomp || omega.size() != ncomp){
    std::cout << "Failed to load component data from `" << databasePath << "`." << std::endl;
    return 1;}

  
  LumpMethod lumpMethod = LumpMethod::Pedersen; 
  //LumpMethod lumpMethod = LumpMethod::Whitson; //usa nC_first e nC_last para definir o número de pseudocomponentes, não é necessário definir npseudos
  //LumpMethod lumpMethod = LumpMethod::Danesh; 

  CorrMethod corrMethod = CorrMethod::WhitsonRiaziDaubert;
  //CorrMethod corrMethod = CorrMethod::InternalAlkaneMM;
  
  int npseudos = 4; //colocar como input da função? 
  int nC_first = 11;
  int nC_last = 18;
  bool useLumping = true;
  //bool useLumping = false; 
  if (useLumping){
  applyLumping(lumpMethod,corrMethod, components, z, MM, flag, Tc, Pc, omega, npseudos, nC_first, nC_last);
  ncomp = z.size();
  kij =  std::vector<std::vector<real>>(ncomp, std::vector<real>(ncomp));  }
 
  zsum = 0.0;
  for (auto zi : z)
  zsum += val(zi);
  std::cout << "Soma(z) = " << zsum << std::endl;


  String outputFileName = "output.out"; // Creating output data file 
  std::ofstream output_file (outputFileName);
  output_file << std::scientific << std::setprecision(10);
  if (!output_file.is_open()){
    std::cout << "Error opening output file: `" << outputFileName << "`." << std::endl; 
    return 0;}// End Creating output data file      
  

  
  IsocurveOpts opts;// Creating options for computations 
  EquilibriumType EquiType;
  auto LnPhiDerivType = LnPhiDerivativeType::Analytical;
  //auto LnPhiDerivType = LnPhiDerivativeType::Automatic;
  //auto LnPhiDerivType = LnPhiDerivativeType::Numerical;

  std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.02,0.03,true,1.e6,true,true);   //std::cout << "Pontos calculados,iniciando bolha;" << std::endl;
  real P0,T0,beta;
  int idx; 
  EquiType = Bubble;
  idx = ncomp; // index for temperature among independent variables, index starts at 0
  T0 = 460.;              
  P0 = 101325.*6; //Pa 12069876.
  real Temp0_Bubble = T0;  real Pres0_Bubble = P0;
  
  //computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts); // Begin Compute the bubble - point curve
     std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
     std::chrono::duration<double> tempo = end - start;
    

     std::chrono::steady_clock::time_point start1 = std::chrono::steady_clock::now(); // Begin Compute the dew - point curve
  //createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.2,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //Original

  // 9° parametro por padrão 0.1, alterado para 0.5 
  createIsocurveOpts(opts,true,500,tolEqui,LnPhiDerivType,3,0.2,1.e-6,0.6,0.03,true,1.e6,true,true); //0.1, diferença    //std::cout << "Pontos calculados,iniciando dew" << std::endl; //pular o ponto crítico
  
  EquiType = Dew;
  idx = ncomp+1; // index for pressure
  T0 = 460.; // K começa em 370 nos gráfcos...243 350antes
  P0 = 101325.; // 101325 Pa =  1 atm
             real Temp0_Dew = T0;  real Pres0_Dew = P0;

  computeIsocurve(output_file,T0,P0,beta,EquiType,z,idx,parameters,opts);  // End Compute the dew - point curve 
  std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
  std::chrono::duration<double> tempo1 = end1 - start1;

 
  std::cout << "Tempo de execução Bubble: " << tempo.count() << " segundos" << ",  T0 ="<< Temp0_Bubble << "  P0 ="<< Pres0_Bubble << std::endl;
  std::cout << "Tempo de execução Dew : " << tempo1.count() << " segundos" <<",  T0 ="<< Temp0_Dew << "  P0 =" << Pres0_Dew << std::endl; //[" << AtoString(LnPhiDerivType) << "]
  
  switch (LnPhiDerivType){
  case Analytical: std::cout << "Dados do tipo Analytical salvos em: output.out" << std::endl;break;
  case Automatic:  std::cout << "Dados salvos em: output_automatica.out" << std::endl;break;
  case Numerical:  std::cout << "Dados salvos em: output_numerica.out" << std::endl;break;default:
  break;}

  if(useLumping){
    switch (lumpMethod){
      case LumpMethod::Pedersen: std::cout << "Lumping method used: Pedersen" << std::endl;
           std::cout << "npseudos = " << npseudos << std::endl;break;
      case LumpMethod::Whitson: std::cout << "Lumping method used: Whitson" << std::endl;break;
      case LumpMethod::Danesh:  std::cout << "Lumping method used: Danesh" << std::endl;break;default:
      break;}

    switch (corrMethod){
      case CorrMethod::WhitsonRiaziDaubert:std::cout << "Correlation method used: WhitsonRiaziDaubert" << std::endl;break;
      case CorrMethod::InternalAlkaneMM:  std::cout << "Correlation method used: InternalAlkaneMM" << std::endl;break;
      break;}
  }



  return 0;


}


