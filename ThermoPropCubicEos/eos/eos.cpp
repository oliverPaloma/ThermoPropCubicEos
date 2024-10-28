#include "eos.hpp"

// Função para calcular a equação de estado com base na seleção do usuário

auto calculateIsoterma(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)->void{
        // Constantes para as equações de estado
 std::string filename = "Arquivos/pressao_T" + std::to_string(static_cast<int>(T)) + ".txt"; //std :: string EoSModel
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo!" << std::endl;
        return;
    }
 outfile << "V(m³/mol)\tP(Pa)\n";  
  auto V = Vi;
  auto inc = (Vf - Vi) / ((double)npoints - 1.0);
  double P;

  for(auto i = 0; i < npoints; i++) 
  {
    calculatePressure(EoSModel,Tc,Pc,omega,T,V,P); 
    outfile << V << "\t" << P << "\n";
    V += inc;
  }
    outfile.close();
    std::cout << "Dados armazenados no arquivo: " << filename << std::endl; 
} 

// Função para calcular a equação de estado com base na seleção do usuário  Eq 3.42
auto calculatePressure(CubicEOSModel EoSModel,std::vector<double> Tc,std::vector<double>Pc,std::vector<double> omega,double T, double V, double &P)-> void { //, double &P
 
   
    



    
    
       auto sigma = computesigma(EoSModel); 
       auto epsilon = computeepsilon(EoSModel); 
       auto psi = computePsi(EoSModel); 
       auto OMEGA = computeOmega(EoSModel); 

       auto Tr = T / Tc[0];
       double alphaTr;
        switch (EoSModel) {
            
            case CubicEOSModel::VanDerWaals: //van der waals vdW case C return 0.0;
                alphaTr = 1.; 
                break;

            case CubicEOSModel::SoaveRedlichKwong: //soave-redlich-kwong SRK    
                alphaTr = pow(1. + (0.480 + 1.574 * omega[0]- 0.176 * omega[0] * omega[0]) * (1. - sqrt(Tr)), 2.);  
                break;

            case CubicEOSModel::PengRobinson: //peng-robinson PR  
                alphaTr = pow(1. + (0.37464 + 1.54226 * omega[0] - 0.26992 * omega[0] * omega[0]) * (1. - sqrt(Tr)), 2.);  
                break;

            default:
                std::cout << "Opção inválida." << std::endl;
                return;
        }
        

       auto b = OMEGA * (R * Tc[0]) / Pc[0]; 
       auto alphaT = psi * (alphaTr * R * R * Tc[0] * Tc[0]) / Pc[0]; 
        P = (R * T) / (V - b) - (alphaT / ((V + epsilon * b) * (V + sigma * b)));
        
        
       
    

    
    
    }

/*

               P = (R * T) / (V - b) - (alphaT / ((V + epsilon * b) * (V + sigma * b)));
                //const auto P = RT/(V - b) - a/((V + e*b)*(V + s*b));

                std::cout << "Pressão ("<< toString(eosModel) <<"): " << P << " Pa" << std::endl;
*/
//van der waals vdW     
//psi = 27.0 / 64.0  
//sigma = 0;  
//epsilon = 0;    
//alphaTr = 1;  
// b = (R * Tc) / (8 * Pc);   
//alphaT = psi * (R * R * Tc * Tc) / Pc;  

                
//soave-redlich-kwong SRK
//psi = 0.42748;
//sigma = 1;
//epsilon = 0;
//Tr = T / Tc;
//alphaTr = pow(1 + (0.480 + 1.574 * omega - 0.176 * omega * omega) * (1 - sqrt(Tr)), 2);  
//b = 0.08664 * (R * Tc) / Pc;  
//alphaT = psi * (alphaTr * R * R * Tc * Tc) / Pc; 

               //double sigma = computesigma(eosModel); 
               //double epsilon = computeepsilon(eosModel); 
               //double psi = computePsi(eosModel); 

                 //double omeg = computeOmega(eosModel); //usado para calcular Z
               
              

    
//peng-robinson PR
//psi = 0.45724;
//sigma = 1 + sqrt(2);
//epsilon = 1 - sqrt(2);
//Tr = T / Tc;
//alphaTr = pow(1 + (0.37464 + 1.54226 * omega - 0.26992 * omega * omega) * (1 - sqrt(Tr)), 2);  
//b = 0.07780 * (R * Tc) / Pc; 
//alphaT = psi * (alphaTr * R * R * Tc * Tc) / Pc;  

/*
switch (EoSModel) {

        case CubicEOSModel::VanDerWaals: return 0.0;
        case CubicEOSModel::RedlichKwong: return 1.0;
        case CubicEOSModel::SoaveRedlichKwong: return 1.0;
        case CubicEOSModel::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
        case 0: // Van der Waals
            beta = omega * (Pc / Tc); //pr=p/pc
            q = (27.0 / 64.0) * (Tc / omega);
            P = 1. + beta - q * beta;
            std::cout << "Z (vdW): " << Z << std::endl;
            break;

        case 1: // Soave-Redlich-Kwong
            beta = 0.08664 * (Pc / Tc);
            q = 0.42748 * (Tc / omega);
            P = 1. + beta - q * beta;
            std::cout << "Z (SRK): " << Z << std::endl;
            break;

        case 2: // Peng-Robinson
            beta = 0.07780 * (Pc / Tc);
            q = 0.45724 * (Tc / omega);
            P = 1. + beta - q * beta;
            std::cout << "Z (PR): " << Z << std::endl;
            break;

        default:
            std::cout << "Opção inválida." << std::endl;
            break;
    }
*/


//=================================================================================================
// == REFERENCE ==
//=================================================================================================
// The implementation of this module is based on the presentation of the textbook:
//
//     Introduction to Chemical Engineering Thermodynamics, 8th Edition, 2017,
//     J.M. Smith, H. Van Ness, M. Abbott, M. Swihart
//
// More specifically, it is based on the information of the following chapters:
//
//      3.6 CUBIC EQUATIONS OF STATE, page 95
//     13.6 RESIDUAL PROPERTIES BY CUBIC EQUATIONS OF STATE, page 87
//
//-------------------------------------------------------------------------------------------------
// For more details and derivation of some formulas, check the document `notes/CubicEOS.lyx`.
//=================================================================================================

/// A high-order function that return an `alpha` function that calculates alpha, alphaT and alphaTT for a given EOS.
auto alpha(CubicEOSModel type) -> Fn<AlphaResult(double, double, double)>
{
    // The alpha function for van der Waals EOS (see Table 3.1 of Smith et al. 2017)
    auto alphaVDW = [](double Tr, double TrT, double omega) -> AlphaResult
    {
        const double alpha = 1.0;
        const double alphaT = 0.0;
        const double alphaTT = 0.0;
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Redlich-Kwong EOS
    auto alphaRK = [](double Tr, double TrT, double omega) -> AlphaResult
    {
        const double alpha = 1.0/sqrt(Tr);
        const double alphaTr = -0.5/Tr * alpha;
        const double alphaTrTr = -0.5/Tr * (alphaTr - alpha/Tr);
        const double alphaT = alphaTr*TrT;
        const double alphaTT = alphaTrTr*TrT*TrT;
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Soave-Redlich-Kwong EOS
    auto alphaSRK = [](double Tr, double TrT, double omega) -> AlphaResult
    {
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
        return { alpha, alphaT, alphaTT };
    };

    // The alpha function for Peng-Robinson (1978) EOS
    auto alphaPR = [](double Tr, double TrT, double omega) -> AlphaResult
    { 


        // Jaubert, J.-N., Vitu, S., Mutelet, F. and Corriou, J.-P., 2005.
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
        return { alpha, alphaT, alphaTT };
    };

    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return alphaVDW;
        case CubicEOSModel::RedlichKwong: return alphaRK;
        case CubicEOSModel::SoaveRedlichKwong: return alphaSRK;
        case CubicEOSModel::PengRobinson: return alphaPR;
        default: return alphaPR;
    }
}

auto computesigma(CubicEOSModel type) -> double //computesigma
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 0.0;
        case CubicEOSModel::RedlichKwong: return 1.0;
        case CubicEOSModel::SoaveRedlichKwong: return 1.0;
        case CubicEOSModel::PengRobinson: return 1.0 + 1.4142135623730951;
        default: return 1.0 + 1.4142135623730951;
    }
}

auto computeepsilon(CubicEOSModel type) -> double
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 0.0;
        case CubicEOSModel::RedlichKwong: return 0.0;
        case CubicEOSModel::SoaveRedlichKwong: return 0.0;
        case CubicEOSModel::PengRobinson: return 1.0 - 1.4142135623730951;
        default: return 1.0 - 1.4142135623730951;
    }
}

auto computeOmega(CubicEOSModel type) -> double 
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 1.0/8.0; 
        case CubicEOSModel::RedlichKwong: return 0.08664;
        case CubicEOSModel::SoaveRedlichKwong: return 0.08664;
        case CubicEOSModel::PengRobinson: return 0.0777960739;
        default: return 0.0777960739;
    }
}

auto computePsi(CubicEOSModel type) -> double
{
    switch(type)
    {
        case CubicEOSModel::VanDerWaals: return 27.0/64.0;
        case CubicEOSModel::RedlichKwong: return 0.42748;
        case CubicEOSModel::SoaveRedlichKwong: return 0.42748;
        case CubicEOSModel::PengRobinson: return 0.457235529;
        default: return 0.457235529;
    }
}

/// Compute the local minimum of pressure along an isotherm of a cubic equation of state.
/// @param a The @eq{a_\mathrm{mix}} variable in the equation of state
/// @param b The @eq{b_\mathrm{mix}} variable in the equation of state
/// @param e The @eq{\epsilon} parameter in the cubic equation of state
/// @param s The @eq{\sigma} parameter in the cubic equation of state
/// @param T The temperature (in K)
/// @return double
auto computeLocalMinimumPressudoubleongIsotherm(double a, double b, double e, double s, 
  double T) -> double
{
    const auto RT = R*T;

    auto V = b;
    auto Pprev = 0.0;

    const auto maxiters = 100;
    const auto tolerance = 1e-6;

    auto i = 0;
    for(; i < maxiters; ++i)
    {
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
auto computeResidualGibbsEnergy(double Z, double beta, double q, double epsilon, 
  double sigma, double T) -> double
{
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
auto determinePhysicalStateThreedoubleRoots(double Zmin, double Zmax, double beta, 
  double q, double epsilon, double sigma, double T) -> StateOfMatter
{
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


auto compute(CubicEOSProps& props, std::vector<double> &Tcr, 
  std::vector<double> &Pcr, std::vector<double> &omega, double T, double P, 
  std::vector<double> &x, CubicEOSModel &model, 
  std::vector<std::vector<double>> &BIP) -> void
{

    /// The number of species in the phase.
    auto nspecies = x.size();

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
    for(auto k = 0; k < nspecies; ++k)
    {
        const double factor = Psi*R*R*(Tcr[k]*Tcr[k])/Pcr[k]; // factor in Eq. (3.45) multiplying alpha
        const auto TrT = 1.0/Tcr[k];
        const auto Tr = T * TrT;
        const auto [alpha, alphaT, alphaTT] = alphafn(Tr, TrT, omega[k]);
        a[k]   = factor*alpha; // see Eq. (3.45)
        aT[k]  = factor*alphaT;
        aTT[k] = factor*alphaTT;
        b[k]   = Omega*R*Tcr[k]/Pcr[k]; // Eq. (3.44)
    }

    // Calculate the parameter `amix` of the phase and the partial molar parameters `abar` of each species
    double amix = {};
    double amixT = {};
    double amixTT = {};
    std::fill(abar.begin(), abar.end(), 0.0);
    std::fill(abarT.begin(), abarT.end(), 0.0);
    for(auto i = 0; i < nspecies; ++i)
    {
        for(auto j = 0; j < nspecies; ++j)
        {
            const double r   = 1.0 - BIP[i][j];
            const double rT  = 0.0;
            const double rTT = 0.0;

            const double s   = sqrt(a[i]*a[j]); // Eq. (13.93)
            const double sT  = 0.5*s/(a[i]*a[j]) * (aT[i]*a[j] + a[i]*aT[j]);
            const double sTT = 0.5*s/(a[i]*a[j]) * (aTT[i]*a[j] + 2*aT[i]*aT[j] + a[i]*aTT[j]) - sT*sT/s;

            const double aij   = r*s;
            const double aijT  = rT*s + r*sT;
            const double aijTT = rTT*s + 2.0*rT*sT + r*sTT;

            amix   += x[i] * x[j] * aij; // Eq. (13.92) of Smith et al. (2017)
            amixT  += x[i] * x[j] * aijT;
            amixTT += x[i] * x[j] * aijTT;

            abar[i]  += 2 * x[j] * aij;  // see Eq. (13.94)
            abarT[i] += 2 * x[j] * aijT;
        }
    }

    // Finalize the calculation of `abar` and `abarT`
    for(auto i = 0; i < nspecies; ++i)
    {
        abar[i] -= amix;
        abarT[i] -= amixT;
    }

    // Calculate the parameters bba[i] and bmix of the cubic equation of state
    //     bbar[i] = Omega*R*Tc[i]/Pc[i] as shown in Eq. (3.44)
    //     bmix = sum(x[i] * bbar[i])
    double bmix = {};
    for(auto i = 0; i < nspecies; ++i)
    {
        bbar[i] = Omega*R*Tcr[i]/Pcr[i]; // see Eq. (13.95) and unnumbered equation before Eq. (13.99)
        bmix += x[i] * bbar[i];  // Eq. (13.91) of Smith et al. (2017)
    }

    // Calculate the temperature and pressure derivatives of bmix
    const auto bmixT = 0.0; // no temperature dependence!
    const auto bmixP = 0.0; // no pressure dependence!

    // Calculate the auxiliary parameter beta and its partial derivatives betaT (at const P) and betaP (at const T)
    const double beta = P*bmix/(R*T); // Eq. (3.46)
    const double betaT = -beta/T; // note bmixT = 0
    const double betaP =  beta/P; // note bmixP = 0

    // Compute the auxiliary variable q and its partial derivatives qT, qTT (at const P) and qP (at const T)
    const double q = amix/(bmix*R*T); // Eq. (3.47)
    const double qT = q*(amixT/amix - 1.0/T); // === amixT/(bmix*R*T) - amix/(bmix*R*T*T)
    const double qTT = qT*qT/q + q*(amixTT/amix - amixT*amixT/(amix*amix) + 1.0/(T*T)); // === qT*(amixT/amix - 1.0/T) + q*(amixTT/amix - amixT*amixT/amix/amix + 1.0/T/T)
    const double qP = 0.0; // from Eq. (3.47), (dq/dP)_T := 0

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

    // Calculate cubic roots using cardano's method
    auto roots = realRoots(cardano(A, B, C));

    // Ensure there are either 1 or 3 double roots!
    assert(roots.size() == 1 || roots.size() == 3);

    // Determine the physical state of the fluid phase for given TPx conditions and its compressibility factor
    double Z = {};

    if(roots.size() == 3)
    {
        const auto Zmax = std::max({roots[0], roots[1], roots[2]});
        const auto Zmin = std::min({roots[0], roots[1], roots[2]});
        props.som = determinePhysicalStateThreedoubleRoots(Zmin, Zmax, beta, q, epsilon, sigma, T);
        Z = (props.som == StateOfMatter::gas) ? Zmax : Zmin;
    }
    else
    {
        props.som = determinePhysicalStateOnedoubleRoot(amix, bmix, epsilon, sigma, T, P);
        Z = roots[0];
    }

    // Calculate ZT := (dZ/dT)_P and ZP := (dZ/dP)_T
    const double ZT = -(AT*Z*Z + BT*Z + CT)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_T = 3*ZZ*ZT + AT*ZZ + 2*A*Z*ZT + BT*Z + B*ZT + CT = 0 => (3*ZZ + 2*A*Z + B)*ZT = -(AT*ZZ + BT*Z + CT)
    const double ZP = -(AP*Z*Z + BP*Z + CP)/(3*Z*Z + 2*A*Z + B); // === (ZZZ + A*ZZ + B*Z + C)_P = 3*ZZ*ZP + AP*ZZ + 2*A*Z*ZP + BP*Z + B*ZP + CP = 0 => (3*ZZ + 2*A*Z + B)*ZP = -(AP*ZZ + BP*Z + CP)

    //=========================================================================================
    // Calculate the integration factor I, IT := (dI/dT)_P and IP := (dI/dP)_T
    //=========================================================================================
    double I = {};
    double IT = {};
    double IP = {};

    if(epsilon != sigma) // CASE I:  Eq. (13.72) of Smith et al. (2017)
    {
        I = log((Z + sigma*beta)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I=\frac{1}{\sigma-\epsilon}\ln\left(\frac{Z+\sigma\beta}{Z+\epsilon\beta}\right) }
        IT = ((ZT + sigma*betaT)/(Z + sigma*beta) - (ZT + epsilon*betaT)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{T}+\sigma\beta_{T}}{Z+\sigma\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
        IP = ((ZP + sigma*betaP)/(Z + sigma*beta) - (ZP + epsilon*betaP)/(Z + epsilon*beta))/(sigma - epsilon); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=\frac{1}{\sigma-\epsilon}\left(\frac{Z_{P}+\sigma\beta_{P}}{Z+\sigma\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
    }
    else // CASE II: Eq. (13.74) of Smith et al. (2017)
    {
        I = beta/(Z + epsilon*beta); // @eq{ I=\frac{\beta}{Z+\epsilon\beta} }
        IT = I*(betaT/beta - (ZT + epsilon*betaT)/(Z + epsilon*beta)); // @eq{ I_{T}\equiv\left(\frac{\partial I}{\partial T}\right)_{P}=I\left(\frac{\beta_{T}}{\beta}-\frac{Z_{T}+\epsilon\beta_{T}}{Z+\epsilon\beta}\right) }
        IP = I*(betaP/beta - (ZP + epsilon*betaP)/(Z + epsilon*beta)); // @eq{ I_{P}\equiv\left(\frac{\partial I}{\partial P}\right)_{T}=I\left(\frac{\beta_{P}}{\beta}-\frac{Z_{P}+\epsilon\beta_{P}}{Z+\epsilon\beta}\right) }
    }

    //=========================================================================================
    // Calculate the ideal volume properties of the phase
    //=========================================================================================
    const double V0  =  R*T/P;
    const double V0T =  V0/T;
    const double V0P = -V0/P;

    //=========================================================================================
    // Calculate the corrected volumetric properties of the phase
    //=========================================================================================
    const auto& V  = props.V  = Z*V0;
    const auto& VT = props.VT = ZT*V0 + Z*V0T;
    const auto& VP = props.VP = ZP*V0 + Z*V0P;

    //=========================================================================================
    // Calculate the residual properties of the phase
    //=========================================================================================
    const auto& Gres  = props.Gres  = R*T*(Z - 1 - log(Z - beta) - q*I); // from Eq. (13.74) of Smith et al. (2017)
    const auto& Hres  = props.Hres  = R*T*(Z - 1 + T*qT*I); // equation after Eq. (13.74), but using T*qT instead of Tr*qTr, which is equivalent
    const auto& Cpres = props.Cpres = Hres/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT); // from Eq. (2.19), Cp(res) := (dH(res)/dT)P === R*(Z - 1 + T*qT*I) + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT) = H_res/T + R*T*(ZT + qT*I + T*qTT*I + T*qT*IT)

    //=========================================================================================
    // Calculate the fugacity coefficients for each species
    //=========================================================================================
    props.ln_phi.resize(nspecies);
    for(auto k = 0; k < nspecies; ++k)
    {
        const double betak = P*bbar[k]/(R*T);
        const double qk    = (1 + abar[k]/amix - bbar[k]/bmix)*q;
        const double Ak    = (epsilon + sigma - 1.0)*betak - 1.0;
        const double Bk    = ((epsilon*sigma - epsilon - sigma)*(2*betak - beta) + qk - q)*beta - (epsilon + sigma - q)*betak;
        const double Ck    = (epsilon*sigma*(2*beta + 1) + 2*q - qk)*beta*beta - (2*(epsilon*sigma + q) + 3*epsilon*sigma*beta)*beta*betak;
        const double Zk    = -(Ak*Z*Z + (B + Bk)*Z + 2*C + Ck)/(3*Z*Z + 2*A*Z + B);

        const double Ik = (epsilon != sigma) ?
            I + ((Zk + sigma*betak)/(Z + sigma*beta) - (Zk + epsilon*betak)/(Z + epsilon*beta))/(sigma - epsilon) :
            I * (1 + betak/beta - (Zk + epsilon*betak)/(Z + epsilon*beta));

        props.ln_phi[k] = Zk - (Zk - betak)/(Z - beta) - log(Z - beta) + q*I - qk*I - q*Ik;
    }

}

//2 funções q calcula isotermas e pressao em funcao do volume e temperatura,// cabeçalho do .hpp
