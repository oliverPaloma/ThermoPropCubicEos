#include "../../../ThermoPropCubicEos.hpp"

int main (){
       auto T = 350.;
       auto P = 100000.;

       std::vector<double> Tcr,Pcr,omega;
      // auto components = "N2 CO2 C1";
      // std::vector<double> x{0.2,0.5,0.3};

//======================Misturas testadas=========================//
       //auto components = "CO2"; 
       //auto components = "C1";  
       //std::vector<double> x{1.};

       //auto components = "CO2 C1";
       //std::vector<double> x{0.2, 0.8};

       //auto components = "C1 C3 NC4";  //gás natural típico
       //std::vector<double> x{0.85, 0.10, 0.05};

       //auto components = "C1 C2 C3 CO2"; // (gás rico com contaminante ácido)
       //std::vector<double> x{0.70, 0.10, 0.10, 0.10};

       //auto components = "C1 C6 C10"; //(mistura leve + médio + pesado)
       //std::vector<double> x{0.60, 0.25, 0.15};
       
       //auto components = "CO2 H2O"; // CO2 + H2O (mistura ácido + água)
       //std::vector<double> x{0.95, 0.05};

       
       auto components = "C1 C7 C14 C20"; // C1 + C7 + C14 + C20 (representando petróleo pesado sintético)
       std::vector<double> x{0.40, 0.25, 0.20, 0.15};

       //auto components = "C1 C2 C3 N2"; // Gás natural típico com nitrogênio
       //std::vector<double> x{0.88, 0.06, 0.04, 0.02};

      
       //auto components = "CO2 C6 C10";  // Frações altas de CO2 com hidrocarbonetos pesados (injeção miscível)
       //std::vector<double> x{0.50, 0.30, 0.20};

       // Mistura ampla de petróleo leve
       //auto components = "C1 NC4 C6 C8";
       //std::vector<double> x{0.50, 0.20, 0.20, 0.10};

       // Mistura pesada típica de C7+ com traços de gás
       //auto components = "C1 C7 C12 C18";
       //std::vector<double> x{0.10, 0.40, 0.30, 0.20};

       // Condensado leve
       //auto components = "C1 C2 IC5 C7";
       //std::vector<double> x{0.55, 0.15, 0.15, 0.15};

       // Mistura com água + hidrocarboneto mais pesado
       //auto components = "H2O C6 C10";
       //std::vector<double> x{0.10, 0.45, 0.45};

       // Mistura rica em etano (gás rico)
       //auto components = "C1 C2 C3 C6";
       //std::vector<double> x{0.50, 0.25, 0.15, 0.10};

       //MIX1 (Hamidreza et al., 2013):
       //auto components = "C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 ";
       //std::vector<double> x{0.413506, 0.040300, 0.215300, 0.053900, 0.054300, 0.051500, 0.051900, 0.103900, 0.003470, 0.002680, 0.002070, 0.001590, 0.001230, 0.000950, 0.000730, 0.000566, 0.000437, 0.001671};             

       //MIX2 (Pedersen et al., 2014):
      // auto components = "N2 CO2 C1 C2 C3 IC4 NC4 IC5 NC5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20";
      // std::vector<double> x{0.001200, 0.024900, 0.764300, 0.074600, 0.031200, 0.005900, 0.012100, 0.005000, 0.005900, 0.007900, 0.009500, 0.010800, 0.007800, 0.005920, 0.004670, 0.003450, 0.003750, 0.003040, 0.002370, 0.002080, 0.002200, 0.001690, 0.001400, 0.008330};


    int nspecies = x.size();
    std::vector<std::vector<double>> BIP(nspecies, std::vector<double>(nspecies));
    
    auto EoSModel = CubicEOSModel::PengRobinson; // PengRobinson   SoaveRedlichKwong  VanDerWaals

    //auto databasePath = "/home/paloma/Área de trabalho/ThermoPropCubicEos/database/test.yml";
    auto databasePath = "/home/paloma/Downloads/ThermoPropCubicEos/database/test.yml";
    read_database(Tcr,Pcr,omega,databasePath,components); //auto x = 1.0; //agora o x vai ser os coeficientes dos compostos.

    //auto epsilon = pow(std::numeric_limits<double>::epsilon(),0.5);
      
auto epsilon_P = 1e-3; //para pressão fotos -3 para todos até c1 c7 c14 c20
   auto epsilon_T = 1e-6; //para temperatura fots -6b para todos até c1 c7 c14 c20.

    auto dt = T*epsilon_T;
    auto dp = P*epsilon_P;
    

    CubicEOSProps props, props1, props2;
    compute(props, Tcr, Pcr, omega, T, P , x, EoSModel, BIP); 
    compute(props1, Tcr, Pcr, omega, T+dt, P, x, EoSModel, BIP); 
    compute(props2, Tcr, Pcr, omega, T,  P+dp, x, EoSModel, BIP); 

    //i < nspecies*nspecies  matricial guardada num vetor
    //i < nspecies vetor
    //i < 1 escalar;


    //auto tamanho = nspecies*nspecies; // matriz
    auto tamanho = nspecies; // vetor 
    //auto tamanho = 1; // escalar

    for (int i = 0; i < tamanho ; ++i) {
         //std::cout << std::fixed << std::setprecision(15);
         auto dN_T = (props1.ln_phi[i]-props.ln_phi[i])/dt;
         auto dN_P = (props2.ln_phi[i]-props.ln_phi[i])/dp;

         //MIX2 (Pedersen et al., 2014)
         //std::cout << "Erro entre Num e anal Temperatura ["<< i  << "]para ["<< components<<"]:   " << abs (dN_T/props.dA_ln_phi_T[i] - 1) << "\n" << std::endl;
         //std::cout << "Erro entre Num e anal Pressão ["<< i  << "]para ["<< components<<"]:   " << abs (dN_P/props.dA_ln_phi_P[i] - 1) << "\n" << std::endl;
         std::cout << "Erro entre Num e anal Temperatura ["<< i  << "]para [MIX2 a T: "<<T<<"]:   " << abs (dN_T/props.dA_ln_phi_T[i] - 1) << "\n" << std::endl;
         std::cout << "Erro entre Num e anal Pressão ["<< i  << "]para [MIX2 a P: "<<P<<"]:   " << abs (dN_P/props.dA_ln_phi_P[i] - 1) << "\n" << std::endl;
         }

return 0; 



}



