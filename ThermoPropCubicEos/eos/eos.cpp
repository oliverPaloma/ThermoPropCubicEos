#include "eos.hpp"   //2 funções q calcula isotermas e pressao em funcao do volume e temperatura,// cabeçalho do .hpp
                     //colocar por ponteiro os vetores!!!!!! para ele nãp criar um vetor novo a cada vez que entra na função.

auto CalcularVolumeIdeal(CubicEOSModel EoSModel, std::vector<double>Tc, std::vector<double>Pc, std::vector<double>omega, std::vector<double>z,int ncomp, double &Vi, double &Vf)-> void{
    auto OMEGA = computeOmega(EoSModel);
    static std::vector<double> b;
    static int ncomp0;
    auto  b_mistura=0.;

       if (ncomp0 != ncomp) {
            b.resize(ncomp);
            ncomp0 = ncomp;}

       for (auto i = 0; i < ncomp; ++i) { 
            b[i] = OMEGA * (R * Tc[i]) / Pc[i];}

       for (auto i = 0; i < ncomp; ++i){     //b_mistura usando regra de mistura linear
            b_mistura += z[i] * b[i];}

    Vi = 10e-10/b_mistura; //Esse volume molar representa o sistema como um gás ideal
    Vf = 0.74/b_mistura;
}

auto calculateIsotermaComp(CubicEOSModel EoSModel, std::vector<double>Tc, std::vector<double>Pc, std::vector<double>omega , double T, std::vector<double>V, std::vector<double>z,int ncomp) -> void {
        std::string filename = "Arquivos/pressao_T" + std::to_string(static_cast<int>(T)) + ".txt"; // colocar equação utilizada 
        std::ofstream outfile(filename);
        
    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return;}

        outfile << "V(m³/mol)\tP(Pa)\n";
        double P;

    for (const auto& volume : V){ 
        calculatePressureComp(EoSModel, Tc, Pc, omega, T, volume , P, z, ncomp);
        //outfile << volume << "\t" << P << "\n";
        outfile << std::fixed << std::setprecision(10) << volume << "\t" << P << "\n";}
        outfile.close();
        std::cout << "Dados armazenados no arquivo: " << filename << std::endl;
}

auto calculatePressureComp(CubicEOSModel EoSModel, std::vector<double>Tc, std::vector<double>Pc, std::vector<double>omega, double T, double V, double &P,  std::vector<double>z, int ncomp) -> void {
        auto sigma = computesigma(EoSModel);
        auto epsilon = computeepsilon(EoSModel);
        auto psi = computePsi(EoSModel);
        auto OMEGA = computeOmega(EoSModel);
        static std::vector<double> b, a, Tr, alphaTr;
        static int ncomp0;
    
    if (ncomp0 != ncomp) {
        a.resize(ncomp);
        b.resize(ncomp);
        Tr.resize(ncomp);
        alphaTr.resize(ncomp);
        ncomp0 = ncomp;}

    for (int i = 0; i < ncomp0; ++i) { 
            Tr[i] = T / Tc[i]; 
            switch (EoSModel) {
                case CubicEOSModel::VanDerWaals: //van der waals vdW case C return 0.0;
                    alphaTr[i] = 1.; 
                    break;
                case CubicEOSModel::SoaveRedlichKwong: //soave-redlich-kwong SRK    
                    alphaTr[i] = pow(1. + (0.480 + 1.574 * omega[i]- 0.176 * omega[i] * omega[i]) * (1. - sqrt(Tr[i])), 2.);   //for
                    break;
                case CubicEOSModel::PengRobinson: //peng-robinson PR  
                    alphaTr[i] = pow(1. + (0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] * omega[i]) * (1. - sqrt(Tr[i])), 2.);   //for
                    break;
                default:
                std::cout << "Opção inválida." << std::endl;
                return;}
     }

    for (auto i = 0; i < ncomp; ++i) { //Tc.size()
            b[i] = OMEGA * (R * Tc[i]) / Pc[i];  //for
            a[i] = psi * (alphaTr[i] * R * R * Tc[i] * Tc[i]) / Pc[i];}

        auto a_mistura=0., b_mistura=0.;

    for (auto i = 0; i < ncomp; ++i) {
        for (auto j = 0; j < ncomp; ++j) {     // Calcula a_mistura usando regra de mistura com k_ij = 0
            a_mistura += z[i] * z[j] * sqrt(a[i] * a[j]);}}

    for (auto i = 0; i < ncomp; ++i){     //b_mistura usando regra de mistura linear
             b_mistura += z[i] * b[i];}

    P = (R * T) / (V - b_mistura) - (a_mistura / ((V + epsilon * b_mistura) * (V + sigma * b_mistura))); 
}

auto calculateIsotermaMisture(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double &Vi, double &Vf, int npoints, std::vector<double>z,int ncomp)->void{
    std::string filename = "Arquivos/pressao_T" + std::to_string(static_cast<int>(T)) + ".txt"; //std :: string EoSModel
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return;}

    outfile << "V(m³/mol)\tP(Pa)\n";  
    auto V = Vi;
    auto inc = (Vf - Vi) / ((double)npoints - 1.0);
    double P;

    for(auto i = 0; i < npoints; i++) {
        calculatePressureMisture(EoSModel,Tc,Pc,omega,T,V,P,z,ncomp); 
        outfile << V << "\t" << P << "\n";
        V += inc;}

    outfile.close();
    std::cout << "Dados armazenados no arquivo: " << filename << std::endl; 
} 

auto calculatePressureMisture(CubicEOSModel EoSModel,std::vector<double> Tc,std::vector<double>Pc,std::vector<double> omega,double T, double V, double &P, std::vector<double>z, int ncomp)-> void { 
       auto sigma = computesigma(EoSModel);  // Função para calcular a equação de estado com base na seleção do usuário  Eq 3.42
       auto epsilon = computeepsilon(EoSModel); 
       auto psi = computePsi(EoSModel); 
       auto OMEGA = computeOmega(EoSModel); 
       static std::vector<double>b,a, Tr, alphaTr;
       
       static int ncomp0;
       
    if(ncomp0 != ncomp){
        a.resize(ncomp);
        b.resize(ncomp);
        Tr.resize(ncomp);
        alphaTr.resize(ncomp);
        ncomp0=ncomp;}
       
    for (int i = 0; i < ncomp0; ++i){ 
       Tr[i] = T / Tc[i]; //for
            switch (EoSModel){
                case CubicEOSModel::VanDerWaals: //van der waals vdW case C return 0.0;
                    alphaTr[i] = 1.; 
                    break;
                case CubicEOSModel::SoaveRedlichKwong: //soave-redlich-kwong SRK    
                    alphaTr[i] = pow(1. + (0.480 + 1.574 * omega[i]- 0.176 * omega[i] * omega[i]) * (1. - sqrt(Tr[i])), 2.);   //for
                    break;
                case CubicEOSModel::PengRobinson: //peng-robinson PR  
                    alphaTr[i] = pow(1. + (0.37464 + 1.54226 * omega[i] - 0.26992 * omega[i] * omega[i]) * (1. - sqrt(Tr[i])), 2.);   //for
                    break;
                default:std::cout << "Opção inválida." << std::endl;
                return;}
     }


    for (auto i = 0; i < ncomp; ++i) { //Tc.size()
        b[i] = OMEGA * (R * Tc[i]) / Pc[i];  //for
        a[i] = psi * (alphaTr[i] * R * R * Tc[i] * Tc[i]) / Pc[i];        }
        
    auto a_mistura=0., b_mistura=0.;

    for (auto i = 0; i < ncomp; ++i){
        for (auto j = 0; j < ncomp; ++j){     // Calcula a_mistura usando regra de mistura com k_ij = 0
            a_mistura += z[i] * z[j] * sqrt(a[i] * a[j]);}
    }

    for (auto i = 0; i < ncomp; ++i) {     //b_mistura usando regra de mistura linear
        b_mistura += z[i] * b[i];}
        
    P = (R * T) / (V - b_mistura) - (a_mistura / ((V + epsilon * b_mistura) * (V + sigma * b_mistura))); //for
    P= P* 9.86923e-6; //conversão de Pa para atm
}


auto calculateIsoterma(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)->void {
    std::string filename = "Arquivos/pressao_T" + std::to_string(static_cast<int>(T)) + ".txt"; 
    std::ofstream outfile(filename);
    
    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return;}

    outfile << "V(m³/mol)\tP(Pa)\n";  

    auto V = Vi;
    auto inc = (Vf - Vi) / ((double)npoints - 1.0);
    double P;

    for(auto i = 0; i < npoints; i++){
        calculatePressure(EoSModel,Tc,Pc,omega,T,V,P); 
        outfile << V << "\t" << P << "\n";
        //outfile << std::fixed << std::setprecision(10) << V << "\t" << P << "\n";
        V += inc;}

    outfile.close();
    std::cout << "Dados armazenados no arquivo: " << filename << std::endl; 
}

auto calculatePressure(CubicEOSModel EoSModel,std::vector<double> Tc,std::vector<double> Pc,std::vector<double> omega,double T, double V, double &P)-> void {
        auto sigma = computesigma(EoSModel); 
        auto epsilon = computeepsilon(EoSModel); 
        auto psi = computePsi(EoSModel); 
        auto OMEGA = computeOmega(EoSModel); 

        auto alphaTr=0;

        auto Tr = T/Tc[0]; 
            switch (EoSModel){
            case CubicEOSModel::VanDerWaals: 
                alphaTr = 1.; 
                break;
            case CubicEOSModel::SoaveRedlichKwong: 
                alphaTr = pow(1. + (0.480 + 1.574 * omega[0]- 0.176 * omega[0] * omega[0]) * (1. - sqrt(Tr)), 2.);   
                break;
            case CubicEOSModel::PengRobinson: 
                alphaTr = pow(1. + (0.37464 + 1.54226 * omega[0] - 0.26992 * omega[0] * omega[0]) * (1. - sqrt(Tr)), 2.);   
                break;
            default:
                std::cout << "Opção inválida." << std::endl;
                return;}

        auto b = OMEGA * (R * Tc[0]) / Pc[0]; 
        auto a = psi * (alphaTr * R * R * Tc[0] * Tc[0]) / Pc[0]; 
        
        P = (R * T) / (V - b) - (a / ((V + epsilon * b) * (V + sigma * b)));}


//===========================================REFERENCE============================================
// The implementation of this module is based on the presentation of the textbook:
//     Introduction to Chemical Engineering Thermodynamics, 8th Edition, 2017,
//     J.M. Smith, H. Van Ness, M. Abbott, M. Swihart
// More specifically, it is based on the information of the following chapters:
//      3.6 CUBIC EQUATIONS OF STATE, page 95
//     13.6 RESIDUAL PROPERTIES BY CUBIC EQUATIONS OF STATE, page 87
// For more details and derivation of some formulas, check the document `notes/CubicEOS.lyx`.




auto alpha(CubicEOSModel type) -> Fn<AlphaResult(double, double, double)>{   ///Especificação dos Parâmetros das Equações de Estado pg 72 /// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT for a given EOS.
    
    auto alphaVDW = [](double Tr, double TrT, double omega) -> AlphaResult{ // The alpha function for van der Waals EOS (see Table 3.1 of Smith et al. 2017)
        const double alpha = 1.0;
        const double alphaT = 0.0;
        const double alphaTT = 0.0;
        return { alpha, alphaT, alphaTT };};

    
    auto alphaRK = [](double Tr, double TrT, double omega) -> AlphaResult{ // The alpha function for Redlich-Kwong EOS
        const double alpha = 1.0/sqrt(Tr);
        const double alphaTr = -0.5/Tr * alpha;
        const double alphaTrTr = -0.5/Tr * (alphaTr - alpha/Tr);
        const double alphaT = alphaTr*TrT;
        const double alphaTT = alphaTrTr*TrT*TrT;
        return { alpha, alphaT, alphaTT };};

    
    auto alphaSRK = [](double Tr, double TrT, double omega) -> AlphaResult{ // The alpha function for Soave-Redlich-Kwong EOS
        const double m = 0.480 + 1.574*omega - 0.176*omega*omega;
        const double sqrtTr = sqrt(Tr);
        const double aux = 1.0 + m*(1.0 - sqrtTr);
        const double auxTr = -0.5*m/sqrtTr;
        const double auxTrTr = 0.25*m/(Tr*sqrtTr);
        const double alpha = aux*aux;
        const double alphaTr = 2.0*aux*auxTr;
        const double alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        const double alphaT = alphaTr * TrT;
        const double alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };};

   
    auto alphaPR = [](double Tr, double TrT, double omega) -> AlphaResult{   // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.  // The alpha function for Peng-Robinson (1978) EOS
        // Extension of the PPR78 model (predictive 1978, Peng–Robinson EOS
        // with temperature dependent kij calculated through a group
        // contribution method) to systems containing aromatic compounds.
        // Fluid Phase Equilibria, 237(1-2), pp.193–211.

        const double m = omega <= 0.491 ?
            0.374640 + 1.54226*omega - 0.269920*omega*omega :
            0.379642 + 1.48503*omega - 0.164423*omega*omega + 0.016666*omega*omega*omega;
        const double sqrtTr = sqrt(Tr);
        const double aux = 1.0 + m*(1.0 - sqrtTr);
        const double auxTr = -0.5*m/sqrtTr;
        const double auxTrTr = 0.25*m/(Tr*sqrtTr);
        const double alpha = aux*aux;
        const double alphaTr = 2.0*aux*auxTr;
        const double alphaTrTr = 2.0*(auxTr*auxTr + aux*auxTrTr);
        const double alphaT = alphaTr * TrT;
        const double alphaTT = alphaTrTr * TrT*TrT;
        return { alpha, alphaT, alphaTT };};

        switch(type){
            case CubicEOSModel::VanDerWaals: return alphaVDW;
            case CubicEOSModel::RedlichKwong: return alphaRK;
            case CubicEOSModel::SoaveRedlichKwong: return alphaSRK;
            case CubicEOSModel::PengRobinson: return alphaPR;
            default: return alphaPR;}
}


auto computesigma(CubicEOSModel type) -> double {
        switch(type){
            case CubicEOSModel::VanDerWaals: return 0.0;
            case CubicEOSModel::RedlichKwong: return 1.0;
            case CubicEOSModel::SoaveRedlichKwong: return 1.0;
            case CubicEOSModel::PengRobinson: return 1.0 + 1.4142135623730951;
            default: return 1.0 + 1.4142135623730951;}
}

auto computeepsilon(CubicEOSModel type) -> double{
        switch(type){
            case CubicEOSModel::VanDerWaals: return 0.0;
            case CubicEOSModel::RedlichKwong: return 0.0;
            case CubicEOSModel::SoaveRedlichKwong: return 0.0;
            case CubicEOSModel::PengRobinson: return 1.0 - 1.4142135623730951;
            default: return 1.0 - 1.4142135623730951;}
}

auto computeOmega(CubicEOSModel type) -> double{
        switch(type){
            case CubicEOSModel::VanDerWaals: return 1.0/8.0; 
            case CubicEOSModel::RedlichKwong: return 0.08664;
            case CubicEOSModel::SoaveRedlichKwong: return 0.08664;
            case CubicEOSModel::PengRobinson: return 0.0777960739;
            default: return 0.0777960739;}
}

auto computePsi(CubicEOSModel type) -> double{
        switch(type){
            case CubicEOSModel::VanDerWaals: return 27.0/64.0;
            case CubicEOSModel::RedlichKwong: return 0.42748;
            case CubicEOSModel::SoaveRedlichKwong: return 0.42748;
            case CubicEOSModel::PengRobinson: return 0.457235529;
            default: return 0.457235529;}
}

/// Compute the local minimum of pressure along an isotherm of a cubic equation of state.
/// @param a The @eq{a_\mathrm{mix}} variable in the equation of state  
/// @param b The @eq{b_\mathrm{mix}} variable in the equation of state
/// @param s The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @return double

auto computeLocalMinimumPressudoubleongIsotherm(double a, double b, double e, double s, double T) -> double{
    const auto RT = R*T;                                                
    auto V = b;                                                                     
    auto Pprev = 0.0;                                                    
    const auto maxiters = 100;                                          
    const auto tolerance = 1e-6;                                          
    auto i = 0;                                          

    for(; i < maxiters; ++i){
        const auto t  = (V + e*b)*(V + s*b);
        const auto tV = 2*V + b*(e + s);

        const auto w   = 1/b * sqrt(RT/(a*tV));
        const auto wV  = -w/tV;
        const auto wVV = -3*wV/tV;

        const auto q   = 1 + t*w - V/b;
        const auto qV  = tV*w + t*wV - 1/b;
        const auto qVV = t*wVV;

        const auto f = q*qV;
        const auto J = qV*qV + q*qVV;

        const auto dV = -f/J;
    
        V += dV;

        const auto P = RT/(V - b) - a/((V + e*b)*(V + s*b));

        if(abs(P - Pprev) < abs(P) * tolerance)
            return abs(q) < abs(qV) ? P : NaN;
        Pprev = P;
    }

    assert(("Could not compute the minimum pressure along an isotherm of a cubic equation of state.",
      (i == maxiters))); 
    return NaN;
}

/// Compute the residual Gibbs energy of the fluid for a given compressibility factor.
/// @param Z The compressibility factor
/// @param beta The @eq{\beta=Pb/(RT)} variable in the cubic equation of state
/// @param q The @eq{q=a/(bRT)} variable in the cubic equation of state
/// @param epsilon The @eq{\epsilon} parameter in the cubic equation of state
/// @param sigma The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)

auto computeResidualGibbsEnergy(double Z, double beta, double q, double epsilon, double sigma, double T) -> double{
    auto I = 0.0;

    if(epsilon != sigma) // CASE I:  Eq. (13.72) of Smith et al. (2017)
        I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma\beta}{Z+\epsilon\beta}\right) }
    else // CASE II: Eq. (13.74) of Smith et al. (2017)
        I = beta/(Z + epsilon*beta); // @eq{ I=\frac{\beta}{Z+\epsilon\beta} }

    const auto Gres = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)
    return Gres;
}

/// Determine the state of matter of the fluid when three double roots are available (either liquid or gas).
/// @param Zmin The compressibility factor with minimum value
/// @param Zmax The compressibility factor with maximum value
/// @param beta The @eq{\beta=Pb/(RT)} variable in the cubic equation of state
/// @param q The @eq{q=a/(bRT)} variable in the cubic equation of state
/// @param epsilon The @eq{\epsilon} parameter in the cubic equation of state
/// @param sigma The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @return StateOfMatter The state of matter of the fluid, by comparing the residual Gibbs energy of the two states.

auto determinePhysicalStateThreedoubleRoots(double Zmin, double Zmax, double beta, double q, double epsilon, double sigma, double T) -> StateOfMatter{
    const auto Gresmin = computeResidualGibbsEnergy(Zmin, beta, q, epsilon, sigma, T);
    const auto Gresmax = computeResidualGibbsEnergy(Zmax, beta, q, epsilon, sigma, T);
    return Gresmin < Gresmax ? StateOfMatter::liquid : StateOfMatter::gas; 
}

/// Determine the state of matter of the fluid when only one double root is available (either supercritical or low pressure gas).
/// @param a The @eq{a_\mathrm{mix}} variable in the equation of state
/// @param b The @eq{b_\mathrm{mix}} variable in the equation of state
/// @param e The @eq{\epsilon} parameter in the cubic equation of state
/// @param s The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @param P The pressure (in Pa)
/// @return StateOfMatter The state of matter of the fluid, by comparing the residual Gibbs energy of the two states.
auto determinePhysicalStateOnedoubleRoot(double a, double b, double e, double s, 
  double T, double P) -> StateOfMatter
{
    const auto Pmin = computeLocalMinimumPressudoubleongIsotherm(a, b, e, s, T);
    return (Pmin != Pmin) ? StateOfMatter::superCritical : (P < Pmin) ? StateOfMatter::gas : StateOfMatter::liquid;
}

auto compute(CubicEOSProps& props, std::vector<double> &Tcr, std::vector<double> &Pcr, std::vector<double> &omega, double T, double P, std::vector<double> &x, CubicEOSModel &model, std::vector<std::vector<double>> &BIP) -> void {
    
    auto nspecies = x.size();// The number of species in the phase.
    /// The function that calculates the interaction parameters kij and its temperature derivatives.
    static std::vector<double> a(nspecies);     ///< Auxiliary array
    static std::vector<double> aT(nspecies);    ///< Auxiliary array
    static std::vector<double> aTT(nspecies);   ///< Auxiliary array
    static std::vector<double> b(nspecies);     ///< Auxiliary array
    static std::vector<double> abar(nspecies);  ///< Auxiliary array
    static std::vector<double> abarT(nspecies); ///< Auxiliary array
    static std::vector<double> bbar(nspecies);  ///< Auxiliary array

    // Check if the mole fractions are zero or non-initialized
    if(nspecies == 0 || *std::max_element(x.begin(),x.end()) <= 0.0)
        return;

    // Auxiliary variables
    const auto Psi = computePsi(model);
    const auto Omega = computeOmega(model);
    const auto epsilon = computeepsilon(model);
    const auto sigma = computesigma(model);
    const auto alphafn = alpha(model);

    // Calculate the parameters `a` and `b` of the cubic equation of state for each species
    for(auto k = 0; k < nspecies; ++k){
        const double factor = Psi*R*R*(Tcr[k]*Tcr[k])/Pcr[k]; // factor in Eq. (3.45) multiplying alpha
        const auto TrT = 1.0/Tcr[k];
        const auto Tr = T * TrT;
        const auto [alpha, alphaT, alphaTT] = alphafn(Tr, TrT, omega[k]);
        a[k]   = factor*alpha; // see Eq. (3.45)
        aT[k]  = factor*alphaT;
        aTT[k] = factor*alphaTT;
        b[k]   = Omega*R*Tcr[k]/Pcr[k]; // Eq. (3.44)
    }

    //teste das derivações a, aT, aTT (ok)
        //props.ln_phi = a;
        //props.dA_ln_phi_T = aT;
        //return;

        //props.ln_phi = aT;
        //props.dA_ln_phi_T = aTT;
        //return;


    // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
    double amix = {};
    double amixT = {};
    double amixTT = {};
    double amixP = {};

    std::fill(abar.begin(), abar.end(), 0.0);
    std::fill(abarT.begin(), abarT.end(), 0.0);

    for(auto i = 0; i < nspecies; ++i){
        for(auto j = 0; j < nspecies; ++j){
            if (i < BIP.size() && j < BIP[i].size()) {
            const double r = 1.0 - BIP[i][j];
            }else{
            std::cerr << "Índices fora dos limites: i = " << i << ", j = " << j << std::endl;}

            const double r   = 1.0 - BIP[i][j];//const double r   = 0;
            const double rT  = 0.0;
            const double rTT = 0.0;
            //const double rP  = 0.0;//add 12/05/25

            const double s   = sqrt(a[i]*a[j]); // Eq. (13.93)
            const double sT  = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
            const double sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;
            //const double sP  = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aP[j]);//add 12/05/25

            const double aij   = r*s;
            const double aijT  = rT*s + r*sT;
            const double aijTT = rTT*s + 2.0*rT*sT + r*sTT;
            //const double aijP  = rP*s + r*sP;//add 12/05/25 não existe 

            amix   += x[i] * x[j] * aij; // Eq. (13.92) of Smith et al. (2017)
            amixT  += x[i] * x[j] * aijT;
            amixTT += x[i] * x[j] * aijTT;
            //amixP  += x[i] * x[j] * aijP;  //add 12/05/25 não existe 

            abar[i]  += 2 * x[j] * aij;  // see Eq. (13.94)
            abarT[i] += 2 * x[j] * aijT;
        }
    }
    
    //teste das derivações
        //props.ln_phi.resize(1);
        //props.dA_ln_phi_P.resize(1);
        //props.ln_phi[0] = abar[0]; //amix com amixT, amixT com amixTT, amix com amixP
        //props.dA_ln_phi_T[0] = abarT[0];
        //props.dA_ln_phi_P[0] = amixP;
        //abar com abarT (ok)
        //props.ln_phi = abar;
        //props.dA_ln_phi_T = abarT;
        //return;
    
    for(auto i = 0; i < nspecies; ++i){ // Finalize the calculation of `abar` and `abarT`
        abar[i] -= amix;
        abarT[i] -= amixT;}

    // Calculate the parameters bba[i] and bmix of the cubic equation of state
    //     bbar[i] = Omega*R*Tc[i]/Pc[i] as shown in Eq. (3.44)
    //     bmix = sum(x[i] * bbar[i])

    double bmix = {};

    for(auto i = 0; i < nspecies; ++i){
        bbar[i] = Omega*R*Tcr[i]/Pcr[i]; // see Eq. (13.95) and unnumbered equation before Eq. (13.99)
        bmix += x[i] * bbar[i];}  // Eq. (13.91) of Smith et al. (2017)
    
    //teste das derivações bbar com bmix (erro)
        //props.ln_phi.resize(1);
        //props.dA_ln_phi_T.resize(1);
        //props.ln_phi = bbar; 
        //props.dA_ln_phi_T.resize(1);
        //props.dA_ln_phi_T[0] = bmix;
        //return;

    const auto bmixT = 0.0; // no temperature dependence! // Calculate the temperature and pressure derivatives of bmix
    const auto bmixP = 0.0; // no pressure dependence!
    //const auto bmixV = 0.0; // add 14/05/25
    
    const double beta = P*bmix/(R*T); // Eq. (3.46) // Calculate the auxiliary parameter beta and its partial derivatives betaT (at const P) and betaP (at const T)
    const double betaT = -beta/T; 
    const double betaP =  beta/P; 

    //teste das derivações beta com betaT       
        //props.ln_phi.resize(1);
        //props.ln_phi[0] = beta; 
        //props.dA_ln_phi_T.resize(1);
        //props.dA_ln_phi_T[0] = betaT;
        //return;

    //teste das derivações beta com betaP       
        //props.ln_phi.resize(1);
        //props.ln_phi[0] = beta; 
        //props.dA_ln_phi_P.resize(1);
        //props.dA_ln_phi_P[0] = betaP;
        //return;

    // Compute the auxiliary variable q and its partial derivatives qT, qTT (at const P) and qP (at const T)
    const double q = amix/(bmix*R*T); // Eq. (3.47)
    const double qT = q*(amixT/amix - 1.0/T); // === amixT/(bmix*R*T) - amix/(bmix*R*T*T)
    const double qTT = qT*qT/q + q*(amixTT/amix - amixT*amixT/(amix*amix) + 1.0/(T*T)); // === qT*(amixT/amix - 1.0/T) + q*(amixTT/amix - amixT*amixT/amix/amix + 1.0/T/T)
    const double qP = 0.0; // from Eq. (3.47), (dq/dP)_T := 0
    const double qV = 0.0; 

    // Convert Eq. (3.48) into a cubic polynomial Z^3 + AZ^2 + BZ + C = 0, and compute the coefficients A, B, C of the cubic equation of state
    const double A = (epsilon + sigma - 1)*beta - 1;
    const double B = (epsilon*sigma - epsilon - sigma)*beta*beta - (epsilon + sigma - q)*beta;
    const double C = -epsilon*sigma*beta*beta*beta - (epsilon*sigma + q)*beta*beta;

    // Calculate AT := (dA/dT)_P, BT := (dB/dT)_P and CT := (dC/dT)_P (i.e., partial derivatives of A, B, C with respect to T at constant P)
    const double AT = (epsilon + sigma - 1)*betaT;
    const double BT = (epsilon*sigma - epsilon - sigma)*(2*beta*betaT) + qT*beta - (epsilon + sigma - q)*betaT;
    const double CT = -epsilon*sigma*(3*beta*beta*betaT) - qT*beta*beta - (epsilon*sigma + q)*(2*beta*betaT);

    // Calculate AP := (dA/dP)_T, BP := (dB/dP)_T and CP := (dC/dP)_T (i.e., partial derivatives of A, B, C with respect to P at constant T)
    const double AP = (epsilon + sigma - 1)*betaP;
    const double BP = (epsilon*sigma - epsilon - sigma)*(2*beta*betaP) + qP*beta - (epsilon + sigma - q)*betaP;
    const double CP = -epsilon*sigma*(3*beta*beta*betaP) - qP*beta*beta - (epsilon*sigma + q)*(2*beta*betaP);

    auto roots = realRoots(cardano(A, B, C));  // Calculate cubic roots using cardano's method

    assert(roots.size() == 1 || roots.size() == 3); // Ensure there are either 1 or 3 double roots!

    double Z = {}; // Determine the physical state of the fluid phase for given TPx conditions and its compressibility factor

    if(roots.size() == 3){
        const auto Zmax = std::max({roots[0], roots[1], roots[2]});
        const auto Zmin = std::min({roots[0], roots[1], roots[2]});
        props.som = determinePhysicalStateThreedoubleRoots(Zmin, Zmax, beta, q, epsilon, sigma, T);
        Z = (props.som == StateOfMatter::gas) ? Zmax : Zmin;
    }else{
        props.som = determinePhysicalStateOnedoubleRoot(amix, bmix, epsilon, sigma, T, P);
        Z = roots[0];}

    // Calculate ZT := (dZ/dT)_P and ZP := (dZ/dP)_T
    const double ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_T = 3*ZZ*ZT + AT*ZZ + 2*A*Z*ZT + BT*Z + B*ZT + CT = 0 => (3*ZZ + 2*A*Z + B)*ZT = -(AT*ZZ + BT*Z + CT)
    const double ZP = -(AP*Z*Z + BP*Z + CP)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_P = 3*ZZ*ZP + AP*ZZ + 2*A*Z*ZP + BP*Z + B*ZP + CP = 0 => (3*ZZ + 2*A*Z + B)*ZP = -(AP*ZZ + BP*Z + CP)

    double I = {};   //=========================================================================
    double IT = {};  // Calculate the integration factor I, IT := (dI/dT)_P and IP := (dI/dP)_T
    double IP = {};  //=========================================================================

    if(epsilon != sigma){ // CASE I:  Eq. (13.72) of Smith et al. (2017)    
        I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma\beta}{Z+\epsilon\beta}\right) }
        IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{T}+\sigma\beta_{T}}{Z+\sigma\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
        IP = ((ZP + sigma*betaP)/(Z + sigma*beta) - (ZP + epsilon*betaP)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{P}+\sigma\beta_{P}}{Z+\sigma\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
    }else{ // CASE II: Eq. (13.74) of Smith et al. (2017)
        I = beta/(Z + epsilon*beta); // @eq{ I=\frac{\beta}{Z+\epsilon\beta} }
        IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta)); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=I\left(\frac{\beta_{T}}{\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
        IP = I*(betaP/beta - (ZP + epsilon*betaP)/(Z + epsilon*beta)); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=I\left(\frac{\beta_{P}}{\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
        }

    const double V0  =  R*T/P; //===================================================//
    const double V0T =  V0/T;  // Calculate the ideal volume properties of the phase//
    const double V0P = -V0/P;  //===================================================//

    const auto& V  = props.V  = Z*V0;          //===========================================================//
    const auto& VT = props.VT = ZT*V0 + Z*V0T; // Calculate the corrected volumetric properties of the phase//
    const auto& VP = props.VP = ZP*V0 + Z*V0P; //===========================================================//

    //===============================================//
    // Calculate the residual properties of the phase//
    //===============================================//
    const auto& Gres  = props.Gres  = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)
    const auto& Hres  = props.Hres  = R*T*(Z - 1 + T*qT*I); // equation after Eq. (13.74), but using T*qT instead of Tr*qTr, which is equivalent
    const auto& Cpres = props.Cpres = Hres/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT); // from Eq. (2.19), Cp(res) := (dH(res)/dT)P === R*(Z - 1 + T*qT*I) + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT) = H_res/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT)

    //=====================================================//
    // Calculate the fugacity coefficients for each species//
    //=====================================================//
    props.ln_phi.resize(nspecies);
    props.dA_ln_phi_T.resize(nspecies);
    props.dA_ln_phi_P.resize(nspecies);
   
    for(auto k = 0; k < nspecies; ++k){
        const double betak = P*bbar[k]/(R*T);
        const double betakT = -P * bbar[k]/(R*T*T) ; //==============Derivações============
        const double betakP = bbar[k]/(R*T);
        const double betakV = 0 ;

        const double qk    = (1 + abar[k]/amix - bbar[k]/bmix)*q;
        const double qkT    = ((-abar[k] / (amix * amix)) * amixT + (bbar[k] / (bmix * bmix)) * bmixT) * q + (1 + abar[k] / amix - bbar[k] / bmix) * qT; //==============Derivações============
        //erro amixP no lugar de 0. tinha amixP
        const double qkP    = ((-abar[k] / (amix * amix)) * 0. + (bbar[k] / (bmix * bmix)) * bmixP) * q + (1 + abar[k] / amix - bbar[k] / bmix) * qP; //==============Derivações============
        //const double qkV    = ((-abar[k] / (amix * amix)) * amixV + (bbar[k] / (bmix * bmix)) * bmixV) * q + (1 + abar[k] / amix - bbar[k] / bmix) * qV ;

        const double Ak  = (epsilon + sigma - 1.0)*betak - 1.0;
        const double AkT = (epsilon + sigma - 1.0) * betakT; //add 12/05/25
        const double AkP = (epsilon + sigma - 1.0) * betakP; //add 12/05/25

        const double Bk    = ((epsilon*sigma - epsilon - sigma)*(2*betak - beta) + qk - q)*beta - (epsilon + sigma - q)*betak;
        const double BkT = ((C * (2*betak - beta) + qk - q) * betaT + (C * (2*betakT - betaT) + qkT - qT) * beta) - (-qT * betak + (epsilon + sigma - q) * betakT);
        const double BkP = ((C * (2*betak - beta) + qk - q) * betaP + (C * (2*betakP - betaP) + qkP - qP) * beta) - (-qP * betak + (epsilon + sigma - q) * betakP);

        const double C = epsilon*sigma - epsilon - sigma;
        const double Ck  = (epsilon*sigma*(2*beta + 1) + 2*q - qk)*beta*beta - (2*(epsilon*sigma + q) + 3*epsilon*sigma*beta)*beta*betak;
        const double CkT = AT * beta * beta + 2.0 * A * beta * betaT - (BT * beta * betak + B * betaT * betak + B * beta * betakT);
        const double CkP = AP * beta * beta + 2.0 * A * beta * betaP - (BP * beta * betak + B * betaP * betak + B * beta * betakP);

        const double Zk    = -(Ak*Z*Z + (B + Bk)*Z + 2*C + Ck)/(3*Z*Z + 2*A*Z + B);
        const double ZkT = -(AkT*Z*Z + (B + BkT)*Z + 2*C + CkT) / (3*Z*Z + 2*A*Z + B); //add 12/05/25
        const double ZkP = -(AkP*Z*Z + (B + BkP)*Z + 2*C + CkP) / (3*Z*Z + 2*A*Z + B); //add 12/05/25

        const double Ik = (epsilon != sigma) ?
            I + ((Zk + sigma*betak)/(Z + sigma*beta) - (Zk + epsilon*betak)/(Z + epsilon*beta))/(sigma - epsilon) : //true
            I * (1 + betak/beta - (Zk + epsilon*betak)/(Z + epsilon*beta)); //false 

        const double IkT = 0.0;
        const double IkP = 0.0;

        props.ln_phi[k] = Zk - (Zk - betak)/(Z - beta) - log(Z - beta) + q*I - qk*I - q*Ik;

        //props.ln_phi_T_perturbada[k] = (T + dx);
        //props.ln_phi_P_perturbada[k] = (P + dx);

        props.dA_ln_phi_T[k] = ZkT - ((Z - beta)*(ZkT - betakT) - (Zk - betak)*(ZT - betaT)) / pow(Z - beta, 2) - (ZT - betaT) / (Z - beta) + qT * I + q * IT - qkT * I - qk * IT - qT * Ik - q * IkT;   // Derivações em 05/05/25 //std::vector<double> ln_phiV;// Derivações em 05/05/25 //std::vector<double> ln_phiV;
        props.dA_ln_phi_P[k] = ZkP - ((Z - beta)*(ZkP - betakP) - (Zk - betak)*(ZP - betaP)) / pow(Z - beta, 2) - (ZP - betaP) / (Z - beta) + qP * I + q * IP - qkP * I - qk * IP - qP * Ik - q * IkP;

        //props.dN_ln_phi_T[k] = (props.ln_phi_T_perturbada[k] - props.ln_phi[k]) / dx;
        //props.dN_ln_phi_P[k] = (props.ln_phi_P_perturbada[k] - props.ln_phi[k]) / dx;

        //std::cout << std::fixed << std::setprecision(15);
        //std::cout << "Erro derivadaT: " << abs (props.dN_ln_phi_T[k]/props.dA_ln_phi_T[k] - 1) << "\n" << std::endl;
        //std::cout << "Erro derivadaP: " << abs (props.dN_ln_phi_P[k]/props.dA_ln_phi_P[k] - 1) << "\n" << std::endl;
        //props.ln_phiV[k] = ZkV - ((Z - beta)(ZkV - betakV) - (Zk - betak)(ZV - betaV)) / pow(Z - beta, 2) - (ZV - betaV) / (Z - beta) + qV * I + q * IV - qkV * I - qk * IV - qV * Ik - q * IkV; 
    }
}

