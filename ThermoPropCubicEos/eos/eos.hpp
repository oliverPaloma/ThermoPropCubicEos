#include "../includes.hpp"
#include "../metodosMatematicos/metodosMatematicos.hpp"

enum class StateOfMatter
{
    superCritical,
    gas,
    liquid
};

/// The properties calculated by evaluating a cubic equation of state.
struct CubicEOSProps
{
    /// The molar volume of the phase (in m3/mol).
    double V = {};
    
    /// The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).
    double VT = {};

    /// The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).
    double VP = {};

    /// The residual molar Gibbs energy of the phase (in J/mol).
    double Gres = {};

    /// The residual molar enthalpy of the phase (in J/mol).
    double Hres = {};

    /// The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    double Cpres = {};

    /// The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    double Cvres = {};

    /// The ln fugacity coefficients of the species in the phase.
    std::vector<double> ln_phi;

    /// The state of matter of the fluid phase
    StateOfMatter som;
};

/// The options for the cubic equation of state models.
enum class CubicEOSModel
{
    VanDerWaals,
    RedlichKwong,
    SoaveRedlichKwong,
    PengRobinson
};

    /// The critical temperatures Tcr of the substances (in K).
    /// The critical pressures Pcr of the substances (in Pa).
    /// The acentric factor of the substances.
auto compute(CubicEOSProps& props, std::vector<double> &Tcr, 
  std::vector<double> &Pcr, std::vector<double> &omega, double T, double P, 
  std::vector<double> &x, CubicEOSModel &model, 
  std::vector<std::vector<double>> &BIP) -> void;

auto alpha(CubicEOSModel type) -> Fn<AlphaResult(double, double, double)>; 

auto computesigma(CubicEOSModel type) -> double;

auto computeepsilon(CubicEOSModel type) -> double;

auto computeOmega(CubicEOSModel type) -> double;

auto computePsi(CubicEOSModel type) -> double;

auto computeLocalMinimumPressudoubleongIsotherm(double a, double b, double e, double s, 
  double T) -> double;

auto computeResidualGibbsEnergy(double Z, double beta, double q, double epsilon, 
  double sigma, double T) -> double;

auto determinePhysicalStateThreedoubleRoots(double Zmin, double Zmax, double beta, 
  double q, double epsilon, double sigma, double T) -> StateOfMatter;  

auto determinePhysicalStateOnedoubleRoot(double a, double b, double e, double s, 
  double T, double P) -> StateOfMatter;  

auto calculateIsoterma(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, double T, double Vi, double Vf, int npoints)->void;
auto calculatePressure(CubicEOSModel EoSModel,std::vector<double>Tc,std::vector<double>Pc,std::vector<double> omega,double T, double V, double &P)-> void; 

auto calculateIsotermaMisture(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double> Pc, std::vector<double> omega, 
                       double T, double Vi, double Vf, int npoints, std::vector<double>z, int ncomp)->void;
auto calculatePressureMisture(CubicEOSModel EoSModel,std::vector<double>Tc,std::vector<double>Pc,std::vector<double> omega,double T, double V, double &P, std::vector<double>z,int ncomp)-> void; 

auto calculateIsotermaComp(CubicEOSModel EoSModel, std::vector<double> Tc, std::vector<double>Pc, std::vector<double> omega,double T, std::vector<double>V, std::vector<double>z,  int ncomp) -> void;

auto calculatePressureComp(CubicEOSModel EoSModel, std::vector<double>Tc,std::vector<double>Pc, std::vector<double>omega, double T, double V, double &P,  std::vector<double>z, int ncomp) -> void;

