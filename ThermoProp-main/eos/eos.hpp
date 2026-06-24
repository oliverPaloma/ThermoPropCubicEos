#pragma once 

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
    real V = {};

    /// The temperature derivative of the molar volume at constant pressure (in m3/(mol*K)).
    real VT = {};

    /// The pressure derivative of the molar volume constant temperature (in m3/(mol*Pa)).
    real VP = {};

    /// The residual molar Gibbs energy of the phase (in J/mol).
    real Gres = {};

    /// The residual molar enthalpy of the phase (in J/mol).
    real Hres = {};

    /// The residual molar heat capacity at constant pressure of the phase (in J/(mol*K)).
    real Cpres = {};

    /// The residual molar heat capacity at constant volume of the phase (in J/(mol*K)).
    real Cvres = {};

    /// The ln fugacity coefficients of the species in the phase.
    std::vector<real> ln_phi;
    
    /// The derivative of ln fugacity coefficients of the species in the phase with respect to temperature.
    std::vector<real> ln_phiT;
    
    /// The derivative of ln fugacity coefficients of the species in the phase with respect to pressure.
    std::vector<real> ln_phiP;   
    
    /// The derivative of ln fugacity coefficients of the species in the phase with respect to pressure.
    std::vector<std::vector<real>> ln_phi_xk;      
    
    //props.ln_phixk[i][k] = jhgiujgihik; 

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

/// Struct for cubic eos parameters
struct param
{
  CubicEOSProps props;
  std::vector<real> Tc;
  std::vector<real> Pc;
  std::vector<real> omega;
  CubicEOSModel model;
  std::vector<std::vector<real>> kij;
};

    /// The critical temperatures Tcr of the substances (in K).
    /// The critical pressures Pcr of the substances (in Pa).
    /// The acentric factor of the substances.
auto compute(CubicEOSProps& props, std::vector<real> &Tcr, 
  std::vector<real> &Pcr, std::vector<real> &omega, real T, real P, 
  std::vector<real> &x, CubicEOSModel &model, 
  std::vector<std::vector<real>> &BIP) -> void;

auto alpha(CubicEOSModel type) -> Fn<AlphaResult(real, real, real)>; 

auto computesigma(CubicEOSModel type) -> real;

auto computeepsilon(CubicEOSModel type) -> real;

auto computeOmega(CubicEOSModel type) -> real;

auto computePsi(CubicEOSModel type) -> real;

auto computeLocalMinimumPressurealongIsotherm(real a, real b, real e, real s, 
  real T) -> real;

auto computeResidualGibbsEnergy(real Z, real beta, real q, real epsilon, 
  real sigma, real T) -> real;

auto determinePhysicalStateThreerealRoots(real Zmin, real Zmax, real beta, 
  real q, real epsilon, real sigma, real T) -> StateOfMatter;  

auto determinePhysicalStateOnerealRoot(real a, real b, real e, real s, 
  real T, real P) -> StateOfMatter;  

auto driver_lnphi(CubicEOSProps& props, std::vector<real> &Tc, 
  std::vector<real> &Pc, std::vector<real> &omega, real T, real P, 
  std::vector<real> &z, CubicEOSModel &EoSModel, 
  std::vector<std::vector<real>> &kij) -> std::vector<real>;

 
