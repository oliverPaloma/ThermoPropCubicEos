#include "../../../ThermoPropCubicEos.hpp"

int main (){
    auto T = 350.;
    auto P = 250.;

    std::vector<double> Tcr,Pcr,omega;
    auto components = "N2 CO2 C1";
    std::vector<double> x{0.25, 0.25, 0.50};//fração molar 
           //auto components = "C1";
           //std::vector<double> x{1.};//fração molar 


    int nspecies = x.size();
    std::vector<std::vector<double>> BIP(nspecies, std::vector<double>(nspecies));
    
    auto EoSModel = CubicEOSModel::PengRobinson; // PengRobinson   SoaveRedlichKwong  VanDerWaals

    //auto databasePath = "/home/paloma/Área de trabalho/ThermoPropCubicEos/database/test.yml";
    auto databasePath = "/home/paloma/Downloads/ThermoPropCubicEos/database/test.yml";
    read_database(Tcr,Pcr,omega,databasePath,components); //auto x = 1.0; //agora o x vai ser os coeficientes dos compostos.

   // auto epsilon = pow(std::numeric_limits<double>::epsilon(),0.5);  
    auto epsilon = 1e-5;
    auto dt = T*epsilon;
    auto dp = P*epsilon;

    CubicEOSProps props, props1, props2;
    compute(props, Tcr, Pcr, omega, T,  P , x, EoSModel , BIP); 
    compute(props1, Tcr, Pcr, omega, T+ dt, P  , x, EoSModel , BIP); 
    compute(props2, Tcr, Pcr, omega, T ,  P+dp , x, EoSModel , BIP); 

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
         std::cout << "Erro entre Num e anal Temperatura ["<< i  << "]:   " << abs (dN_T/props.dA_ln_phi_T[i] - 1) << "\n" << std::endl;
         std::cout << "Erro entre Num e anal Pressão ["<< i  << "]:   " << abs (dN_P/props.dA_ln_phi_P[i] - 1) << "\n" << std::endl;
         }

return 0; 
}



