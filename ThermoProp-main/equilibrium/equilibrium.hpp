#pragma once 

#include "../includes.hpp"
#include "../eos/eos.hpp"

const static real tolEqui = 1.e-9;

enum EquilibriumType {Bubble,Dew,FractionalBeta};
enum LnPhiDerivativeType {Analytical, Automatic, Numerical};


//testes 
//const char* toString(LnPhiDerivativeType t)
//{
    //switch (t)
    //{
    ///    case Analytical: return "Analytical";
       //case Automatic:  return "Automatic";
        //case Numerical:  return "Numerical";
       /// default:         return "Unknown";
    //}
///}



struct IsocurveOpts
{
  // Settings for Newton iterations
  int MaxNewtonIter;
  real NewtonTol;

  // Evaluates the Jacobian Matrix Analytically if it is true, otherwise the Jacobian Matrix is calculated Numerically
  bool EvalAnalyticalJacobian;

  // Evaluates the Sensibility Matrix if it is true, otherwise it is not calculated
  bool ChangeSpecfVarBasedOnSensi;

  // Selects how to evaluate lnphi derivatives (analytical, automatic differentiation, or numerical)
  LnPhiDerivativeType LnPhiDerivType;

  // Target number of iterations ( for selecting step size )
  int TargetIter;
  
  // Maximum change in specification
  real MaxDeltaSpec;

  // Maximum percentage variation in specified variable
  real MaxPercVarSpecVariable;  

  // Minimum change in specification
  real MinDeltaSpec;  

  // Tolerance for the detection of the critical point
  real CritPointTol;
  
  // Whether to use polynomial regression to generate initial estimates or not
  bool UseCubicInterpolation;

  // Maximum condition number of system matrix in polynomial regression
  real MaxCondInterpol;

  // Set verbosity
  bool Verbose;  

};

struct FlashOpts //essas são as mesmas que em recharacterize_cubic?
{
  // Settings for Newton iterations
  int MaxNewtonIter;
  real NewtonTol;

  // Evaluates the Jacobian Matrix Analytically if it is true, otherwise the Jacobian Matrix is calculated Numerically
  bool EvalAnalyticalJacobian;

  // Selects how to evaluate lnphi derivatives (analytical, automatic differentiation, or numerical)
  LnPhiDerivativeType LnPhiDerivType;

  // Set verbosity
  bool Verbose;
};

auto createFlashOpts (FlashOpts &opts, bool EvalAnalyticalJacobian, int MaxNewtonIter, real NewtonTol,
  LnPhiDerivativeType LnPhiDerivType, bool Verbose) -> void;

auto stabilityAnalysis(real T0, real &T, real dT, real P, int niter, Vec<real> &z, Vec<real> &w, 
  param &parameters) -> void;

auto flashMichelsen(Vec<real> &lnK, real &T, real &P, Vec<real> &x, Vec<real> &y,
  Vec<real> &z, real &beta, param &parameters, int idx, int ncomp, FlashOpts &opts) -> void;   

// The isothermal two-phase flash algorithm described by (Chap. 10):
// Michelsen, M., Mollerup, J., 2007. Thermodynamic Models: Fundamentals and Computational
// Aspects. Tie-Line Publications.
auto createIsocurveOpts(IsocurveOpts &opts, bool EvalAnalyticalJacobian, int MaxNewtonIter, real NewtonTol,
  LnPhiDerivativeType LnPhiDerivType, int TargetIter, real MaxDeltaSpec, real MinDeltaSpec, real MaxPercVarSpecVariable,
  real CritPointTol, bool UseCubicInterpolation, real MaxCondInterpol, bool ChangeSpecfVarBasedOnSensi, bool Verbose) -> void;

// Resolve Newton's Method for Michelsen Isocurve Algorithm
// Michelsen, M., Mollerup, J., 2007. Thermodynamic Models: Fundamentals and Computational
// Aspects. Tie-Line Publications. 
auto ResolveNewtonMichelsenAlgorithm(Vec<real> &X, Vec<real> &dXdS, int &k, real &mix_err,
  real SS, param &parameters, int nvar, FlashOpts &opts, int ncomp,
  real beta, Vec<real> &z, int idx, int iter_number, bool &convergeDensity, LnPhiDerivativeType LnPhiDerivType) -> void;

// It Computes Isocurves 
// beta = 0 -> bubble line
// beta = 1 -> dew line
//
// Returns T (temp.), P(Press.), K(y/x), x(liq. mol.frac), y(vap. mol.frac) vectors for 
// L/V equilibrium
auto computeIsocurve(std::ofstream &file, real T0, real P0, real &beta,
  EquilibriumType &EquiType,Vec<real> &z, int idx, param &parameters,
  IsocurveOpts &opts) -> void;


// Polynomial approximation to compute the initial guess for the subsequent Newton iterations
// Using Paulo Lage's description in Trello:
// https://trello.com/c/iWe1kyn3/287-ajuste-dos-par%C3%A2metros-da-pcsaft-a-misturas-de-hidrocarbonetos
auto PolynomialInitialEstimate(String &InterpolScheme, Vec<real> &X0 , Vec<real> &X1,
  Vec<real> &X2, real SS0, real SS1, real SS2, Vec<real> &dXdS1, 
  Vec<real> &dXdS0) -> void;  

// Linear approximation to compute the initial guess for the subsequent Newton iterations
// for the isothermal two-phase flash algorithm described by (Chap. 10):
// Michelsen, M., Mollerup, J., 2007. Thermodynamic Models: Fundamentals and Computational
// Aspects. Tie-Line Publications.   
auto LinearInitialStimate(Vec<real> &X, Vec<real> &X0, real dS, Vec<real> &dXdS, int n) 
  -> void;

// It calculates F(X;S), dF/dX, and dF/dS
// for the isothermal two-phase flash algorithm described by (Chap. 10):
// Michelsen, M., Mollerup, J., 2007. Thermodynamic Models: Fundamentals and Computational
// Aspects. Tie-Line Publications. 
auto IsocurveEquations(Vec<real> &FF, Vec<real> &dFdX, Vec<real> &dFdS, Vec<real> &X,
  real beta, Vec<real> &z, real SS, int idx, param &parameters,
  bool calc_dFdX, int ncomp, bool &convergeDensity, LnPhiDerivativeType LnPhiDerivType) -> void;