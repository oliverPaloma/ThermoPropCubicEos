#include "equilibrium.hpp"

auto createFlashOpts (FlashOpts &opts, bool EvalAnalyticalJacobian, int MaxNewtonIter, real NewtonTol,
  LnPhiDerivativeType LnPhiDerivType, bool Verbose) -> void
{
  opts.MaxNewtonIter = MaxNewtonIter;
  opts.EvalAnalyticalJacobian = EvalAnalyticalJacobian;
  opts.LnPhiDerivType = LnPhiDerivType;
  opts.NewtonTol = NewtonTol;
  opts.Verbose = Verbose;
}

auto flashMichelsen(Vec<real> &lnK, real &T, real &P, Vec<real> &x, Vec<real> &y,
  Vec<real> &z, real &beta, param &parameters, int idx, int ncomp, FlashOpts &opts) -> void
{
  // T, P, and lnK enter as initial guesses (Specified variable is the exception). 
  // The specified variable will have indx "idx".

  if (opts.Verbose)   std::cout << std::scientific << std::setprecision(8);

  auto a = lnK.size();
  auto b = z.size();
  auto c = x.size();
  auto d = y.size();
  assert(("lnK, z, y, x vectors should have the same size",((a == b) && (b == c) && (c == d))));
  assert(("The number of components should be the size of lnK, z, y, x vectors",((a == ncomp))));

  std::cout << "\nStarting Flash calculation with Beta = " << beta <<"\nGood initial estimates are required for liquid and vapor molar fractions" << std::endl; 

  // Maximum number of Newton iterations
  auto kmax = opts.MaxNewtonIter;

  // Tolerance for the Newton iterations
  auto abs_error = opts.NewtonTol;

  int i;
  real aux, mix_err;
  auto nvar = ncomp+2;
  static Vec<real> X(nvar), dXdS(nvar), deltaX(nvar), FF(nvar), dFdX(nvar*nvar), dFdS(nvar);
  static Vec<int> idxLU(nvar);
  X.resize(nvar);
  deltaX.resize(nvar);
  FF.resize(nvar);
  dFdX.resize(nvar*nvar);
  dFdS.resize(nvar);
  idxLU.resize(nvar);

  for (i = 0; i < ncomp; i++) X[i] = log(y[i]/x[i]); // initial guess for lnK
  X[ncomp] = log(T);
  X[ncomp+1] = log(P);

  auto SS = X[idx];

  // Evaluate function and Jacobian
  static bool convergeDensity;
  IsocurveEquations(FF,dFdX,dFdS,X,beta,z,SS,idx,parameters,true,ncomp,convergeDensity,opts.LnPhiDerivType);
  if (!convergeDensity) 
  {
    std::cout << "Density did not converge during flash calculation" << std::endl;
    return;
  }

  auto Converged = false;
  auto Diverged = false;
  auto k = 0;


  while ((!Converged) && (!Diverged))
  {   
    // Increment iteration counter
    k++;

    // Newton update
    // LU decomposition and Back Substitution to solve the linear system:
    // dFdX * deltaX = -FF
    for (i = 0; i < nvar; i++) deltaX[i] = FF[i];
    ludcmp(dFdX,nvar,idxLU,aux);
    lubksb(dFdX,nvar,idxLU,deltaX); 
    for (i = 0; i < nvar; i++) X[i] -= deltaX[i]; 

    //  Check for convergence
    mix_err = mixed_error(deltaX,X,abs_error);  
    Converged = (mix_err < 1.);
    Diverged = (k >= kmax);

    // Evaluate function and Jacobian Analytically using new X
    if ((!Converged) && (!Diverged)) IsocurveEquations(FF,dFdX,dFdS,X,beta,z,SS,idx,parameters,true,ncomp,convergeDensity,opts.LnPhiDerivType);
  }
  if (opts.Verbose)
  { 
    std::cout << "Variable updated index: " << idx << std::endl; 
    std::cout << "Mixed error in Newton's method (should be < 1): " << mix_err << std::endl;
    std::cout << "Nº of Newton iterations: " << k << std::endl; 
    for (i = 0; i < ncomp; i++) std::cout << "lnK[" << i << "]\t\t";
    std::cout << "T (K)\t\t\tP (MPa)" << std::endl;
    for (i = 0; i < ncomp; i++) std::cout << X[i] << "\t";
    std::cout << exp(X[ncomp]) << "\t\t" << exp(X[ncomp+1])/1.e6 << "\n" << std::endl; 
  }    
  

  for (i = 0; i < ncomp; i++) 
  {
    lnK[i] = X[i];
    aux = exp(lnK[i]);
    x[i] = z[i]/(1.-beta+beta*aux);
    y[i] = aux*x[i]; 
  }
  T = exp(X[ncomp]);
  P = exp(X[ncomp+1]);  
  

  if (k == kmax) std::cout << "Flash calculation did not converge" << std::endl; 
  else std::cout << "Flash calculation converged for T = " << T << " K " << " and P = " << P/1.e6 << " MPa" << std::endl; 
}


auto stabilityAnalysis(real T0, real &T, real dT, real P, int niter, Vec<real> &z, Vec<real> &w, 
  param &parameters) -> void
{
  auto& [props, Tc, Pc, omega, model, kij] = parameters; 

  T = T0;
  real tdp = 1.;
  auto ncomp = z.size();
  static Vec<real> h, lnW, W, w0; 
  h.resize(ncomp);
  lnW.resize(ncomp);
  W.resize(ncomp);
  w0.resize(ncomp);
  auto countEstab = 0;
  real sum, err;
  int i,j;
  auto nsubiter = 50;

  std::fill(w.begin(),w.end(),1./ncomp); 
  std::fill(w0.begin(),w0.end(),1.);
  
  while((tdp > 0.) && (countEstab < niter))
  { 
    compute(props,Tc,Pc,omega,T,P,z,model,kij);  
    for(auto i = 0; i < ncomp; i++) h[i] = log(z[i]) + props.ln_phi[i];

    std::cout << "\nIteration " << countEstab << std::endl;

    // Converging w

    // Succesive substitution
  for(j = 0; j < nsubiter; j++)
    { 
      sum = 0.;
      compute(props,Tc,Pc,omega,T,P,w,model,kij);  
      for(i = 0; i < ncomp; i++)
      {
        lnW[i] = h[i] - props.ln_phi[i]; 
        W[i] = exp(lnW[i]);  
        sum += W[i];    
      }
      err = 0.;
      for(i = 0; i < ncomp; i++) 
      { 
        w[i] = W[i]/sum;
        err += abs((w[i]-w0[i]));
      }
      if (err < tolEqui) break;
      for(i = 0; i < ncomp; i++) w0[i] = w[i];
    }   

   
    // Verifying if converged to trivial solution (usually one phase region)
    auto verTrivSol = false;
    auto logdwz = w;
    for(i = 0; i < ncomp; i++) logdwz[i] = log(w[i]) - log(z[i]);
    if (norm(logdwz,0) < 1.e-3) 
    {
      std::fill(w.begin(),w.end(),1./ncomp);  
      verTrivSol = true;
    }

    if (j == nsubiter) std::cout << "\tMax. number of subiterations (" << nsubiter << ") was reached\n\tNorm of error vector in incipient molar fraction = "<< err <<"\n\tIncipient molar fraction did not converge during stability analysis for T = " << T << " K, Let's update the temperature" << std::endl;
    else std::cout << "\tIteration converged with (" << j+1 << ") subiterations\n\tNorm of error vector in incipient molar fraction = "<< err << std::endl;

    countEstab++;
    T += dT;    
    
    // Jumping Bad solutions
    if (err > 1.e10) {std::cout << "\tJumping bad solution: error in calculating w is = " << err
      << ", next T = " << T << " K" << std::endl; continue;}     
    if (verTrivSol) {std::cout << "\tJumping bad solution: Trivial solution w = z, next T = " << T << " K" << std::endl; continue;}     


    tdp = 0.; 
    for(i = 0; i < ncomp; i++) tdp += W[i]*(log(W[i]) + props.ln_phi[i] - h[i] - 1.); 
    tdp += 1.;    
    
    std::cout << "\tT = " << T-dT << " , tdp = " << tdp << std::endl;
  }

  T -= dT;

  compute(props,Tc,Pc,omega,T,P,w,model,kij);  

  if (countEstab==niter) std::cout << "\nStability analysis reached the maximum number of iterations: " << countEstab <<
    ". The analyzed phase is stable in the provided temperature interval: " << T0 << " - " << T << " K" << std::endl;
  else
  {
    std::cout << "\nA second phase has been found at around T =  " << T << " K,\nMolar volume = " << props.V << " m³/mol \nThe approximated molar fraction: " << std::endl;
    for(i = 0; i < ncomp; i++) std::cout << "w[" << i << "] = " << w[i] << std::endl;
  }
  
}

auto createIsocurveOpts (IsocurveOpts &opts, bool EvalAnalyticalJacobian, int MaxNewtonIter, real NewtonTol,
  LnPhiDerivativeType LnPhiDerivType, int TargetIter, real MaxDeltaSpec, real MinDeltaSpec, real MaxPercVarSpecVariable,
  real CritPointTol, bool UseCubicInterpolation, real MaxCondInterpol, bool ChangeSpecfVarBasedOnSensi, bool Verbose) -> void
{
  opts.MaxNewtonIter = MaxNewtonIter;
  opts.ChangeSpecfVarBasedOnSensi = ChangeSpecfVarBasedOnSensi;
  opts.EvalAnalyticalJacobian = EvalAnalyticalJacobian;
  opts.LnPhiDerivType = LnPhiDerivType;
  opts.NewtonTol = NewtonTol;
  opts.TargetIter = TargetIter;
  opts.MaxDeltaSpec = MaxDeltaSpec;
  opts.MinDeltaSpec = MinDeltaSpec;
  opts.MaxPercVarSpecVariable = MaxPercVarSpecVariable; // if (MaxPercVarSpecVariable > 1.) dS will not be limited by maxX 
  opts.CritPointTol = CritPointTol;
  opts.UseCubicInterpolation = UseCubicInterpolation;
  opts.MaxCondInterpol = MaxCondInterpol;
  opts.Verbose = Verbose;  
}  

auto ResolveNewtonMichelsenAlgorithm(Vec<real> &X, Vec<real> &dXdS, int &k, real &mix_err,
  real SS, param &parameters, int nvar, FlashOpts &opts, int ncomp,
  real beta, Vec<real> &z, int idx, int iter_number, bool &convergeDensity, LnPhiDerivativeType LnPhiDerivType) -> void
{
  int i;
  real aux;
  auto abs_error = opts.NewtonTol;
  auto kmax = opts.MaxNewtonIter;
  static Vec<real> deltaX(nvar), FF(nvar), dFdX(nvar*nvar), dFdS(nvar);
  static Vec<int> idxLU(nvar);
  deltaX.resize(nvar);
  FF.resize(nvar);
  dFdX.resize(nvar*nvar);
  dFdS.resize(nvar);
  idxLU.resize(nvar);

  // Evaluate function and Jacobian
  IsocurveEquations(FF,dFdX,dFdS,X,beta,z,SS,idx,parameters,true,ncomp,convergeDensity,LnPhiDerivType);

  auto Converged = (norm(FF,0) < abs_error);
  auto Diverged = false;
  k = 0;
  
  do
  {   
    // Increment iteration counter
    k++;

    // Newton update
    // LU decomposition and Back Substitution to solve the linear system:
    // dFdX * deltaX = -FF
    for (i = 0; i < nvar; i++) deltaX[i] = FF[i];
    ludcmp(dFdX,nvar,idxLU,aux);
    lubksb(dFdX,nvar,idxLU,deltaX); 
    for (i = 0; i < nvar; i++) X[i] -= deltaX[i]; 

    //  Check for convergence
    mix_err = mixed_error(deltaX,X,abs_error);  
    Converged = (mix_err < 1.);
    Diverged = (k >= kmax);

    // Evaluate function and Jacobian Analytically using new X
    if ((!Converged) && (!Diverged)) IsocurveEquations(FF,dFdX,dFdS,X,beta,z,SS,idx,parameters,true,ncomp,convergeDensity,LnPhiDerivType);
  } while ((!Converged) && (!Diverged));
  

  // Compute sensitivities     
  // LU decomposition and Back Substitution to solve the linear system:
  // dFdX * dXdS = -dFdS
  //IsocurveEquations(FF,dFdX,dFdS,X,beta,z,SS,idx,parameters,true,ncomp,convergeDensity); 
  //ludcmp(dFdX,nvar,idxLU,aux);
  lubksb(dFdX,nvar,idxLU,dFdS);  // dFdS exits as dXdS, dFdX comes from Newton's Method already modified by LU decomposion   
  for (i = 0; i < nvar; i++) dXdS[i] = -dFdS[i]; 

  if (Diverged) 
  {
    std::cout << "\nNewton's method did not Converge for iteration: " 
    << iter_number << std::endl;
    std::cout << "Mixed error =  " << mix_err << std::endl;
  }

}

auto computeIsocurve(std::ofstream &file, real T0, real P0, real &beta, 
  EquilibriumType &EquiType,Vec<real> &z, int idx, param &parameters, 
  IsocurveOpts &opts) -> void
{
  auto& [props, Tc, Pc, omega, model, kij] = parameters;
  
  // If sensibilities are not calculated, the specified variable does not change

  //assert(("Temperature, Pressure, lnK, x, and y vectors should be empty",
  //  (T.empty() && P.empty() && lnK.empty() && x.empty() && y.empty())));

  assert(("z vectors should have non zero size",
    (!z.empty()))); 

  assert(("Initial Temperature and Pressure should be positives",
    ((T0 > 0) && (P0 > 0)))); 

  // Maximum number of Newton iterations
  auto kmax = opts.MaxNewtonIter;

  // Tolerance for the Newton iterations
  auto abs_error = opts.NewtonTol;

  //Target number of iterations
  auto targetIter = opts.TargetIter;

  // Tolerance for the detection of the critical point
  auto crittol = opts.CritPointTol;

  // Maximum change in specification
  auto maxdS = opts.MaxDeltaSpec;  

  // Minimum change in specification
  auto mindS = opts.MinDeltaSpec;    

  // Max % of variation in the specified variable
  // if (maxX > 1.) dS will not be limited by maxX 
  auto maxX = opts.MaxPercVarSpecVariable; 

  // Whether to use polynomial regression or not
  auto UseCubicInterpolation = opts.UseCubicInterpolation; 

  file << std::scientific << std::setprecision(16);

  // Creating log file
  static String logFileName;
  if (EquiType == Bubble) logFileName = "isocurveBubble.log"; 
  else if (EquiType == Dew) logFileName = "isocurveDew.log"; 
  else logFileName = "isocurveFractionalBeta.log"; 
  std::ofstream logFile (logFileName);
  logFile << std::scientific << std::setprecision(5);
  if (!logFile.is_open()) 
  {
    std::cout << "Error opening log file: `" << logFileName << "`" << std::endl;
    return;
  }  

  // Initial Settings
  auto ncomp = z.size();
  int i,j,k;
  
  // normalizing z vector 
  real sumZ = 0.;
  for (auto i = 0; i < ncomp; i++) sumZ += z[i];
  for (auto i = 0; i < ncomp; i++) z[i] /= sumZ;


  int idxNew,MaxIterStability = 1000;
  real dTGetPhase,fac = 1.;
  auto T = T0; // just to create first position, T[0] will be calculated from P[0]
  auto P = P0; 
  static Vec<real> x(ncomp),y(ncomp),lnK(ncomp);
  x.resize(ncomp);
  y.resize(ncomp);
  lnK.resize(ncomp);

  FlashOpts optsFlash;
  optsFlash.EvalAnalyticalJacobian = opts.EvalAnalyticalJacobian;
  optsFlash.LnPhiDerivType = opts.LnPhiDerivType;
  optsFlash.MaxNewtonIter = opts.MaxNewtonIter;
  optsFlash.NewtonTol = opts.NewtonTol;
  optsFlash.Verbose = opts.Verbose;

  switch (EquiType)
  {
    case Bubble:
      beta = 0.; 
      logFile << "\nBubble Line Calculations\nPhase stability calculation to find initial guesses Temperature and Vapor Molar Fraction" << std::endl;                

      // Compute initial guess for equilibrium constants  
      dTGetPhase = 0.5;// 0.1 mais rapido //.005 trocado para .5 em 05/06/26

      // it calculates initial guesses for y and temperature for P0
      stabilityAnalysis(T0,T,dTGetPhase,P,MaxIterStability,z,y,parameters);  //05/06/26
      for (i = 0; i < ncomp; i++) lnK[i] = log((y[i]/z[i]));
      //y[0] = 0.99; y[1] = 0.01; y[2] = 0.0001;
      //for (i = 0; i < ncomp; i++) lnK[i] = log((y[i]/z[i])); 
     
      logFile << "\nStarting Bubble Line (beta = 0) calculations\n" << std::endl;         

      break;

    case Dew:
      beta = 1.; 
      logFile << "\nDew Line Calculations\nPhase stability calculation to find initial guesses Temperature and Liquid Molar Fraction" << std::endl;                

      // Compute initial guess for equilibrium constants 
      dTGetPhase = -1.;

      // it calculates initial guesses for x and temperature for P0
      stabilityAnalysis(T0,T,dTGetPhase,P,MaxIterStability,z,x,parameters);   
      for (i = 0; i < ncomp; i++) lnK[i] = log((z[i]/x[i]));    
      logFile << "\nStarting Dew Line (beta = 1) calculations\n" << std::endl;                
      break;    

    case FractionalBeta:
    {
      logFile << "\nEquilibrium Calculation for Beta = " << beta << std::endl;                
      
      // Compute initial guess for equilibrium constants  
      dTGetPhase = 1.;

      // Initial guesses for y and x 
      real Tdew, Tbub;
      stabilityAnalysis(T0-10.,Tbub,dTGetPhase,P,MaxIterStability,z,y,parameters);    
      stabilityAnalysis(T0+10.,Tdew,-dTGetPhase,P,MaxIterStability,z,x,parameters);   
      
      // Flash at constant beta
      flashMichelsen(lnK,T,P,x,y,z,beta,parameters,idx,ncomp,optsFlash);  
      break;
    }
    default :
      logFile << "Please enter a valid type of equilibrium calculation\nDew, Bubble or FracFractionalBeta";
      return;
  }

  // Initial guess for vector of independent variable
  auto nvar = ncomp+2; // number of independent variables
  static Vec<real> X(nvar), Xplus(nvar), Xminus(nvar), lnK_aux(ncomp), maxAbsdXdS(nvar),
    dXdS(nvar), dXdSminus(nvar);
  X.resize(nvar);
  Xplus.resize(nvar);
  Xminus.resize(nvar);
  lnK_aux.resize(ncomp);
  maxAbsdXdS.resize(nvar);
  dXdS.resize(nvar);
  dXdSminus.resize(nvar);  
  static Vec<int> indx;
  indx.resize(ncomp);
  for (i = 0; i < ncomp; i++) X[i] = lnK[i];
  X[ncomp] = log(T);
  X[ncomp+1] = log(P);

  //Initial specification
  real SS,SSplus,SSminus,aux,mix_err,maxlnK; // SS is S in michelsen algorithm, letter S is already being used
  SS = X[idx];
  static String InterpolScheme;
 
  // Initial step
  auto dS = maxdS; 

  auto count = 0;
  auto CriticalPointReached = false;
  auto MaxEquilIter = 10000;
  static bool convergeDensity, SpecifiedVariableUpdated;
  SpecifiedVariableUpdated = false;
  static real deterQ, C, B; // 2º (B) and 3º (C) terms of helmholtz free energy expansion in taylor series
  while ((!CriticalPointReached) && (count < MaxEquilIter))
  {

    // Display message
    logFile << "Iteration: " << count << std::endl;    
    if (opts.Verbose)
    { 
      logFile << "Index of specified variable before Newton's method: " << idx << std::endl; 
      logFile << "Temperature before Newton's method: " << exp(X[ncomp]) << " K" << std::endl; 
      logFile << "Pressure before Newton's method: " << exp(X[ncomp+1]) << " Pa" << std::endl; 
    }

    // Updating X by Newton's Method
    ResolveNewtonMichelsenAlgorithm(X,dXdS,k,mix_err,SS,parameters,nvar,optsFlash,ncomp,
      beta,z,idx,count,convergeDensity,opts.LnPhiDerivType);

    // Updating T,P,lnK,x,y
    for (i = 0; i < ncomp; i++) 
    {
      aux = exp(X[i]);
      lnK[i] = X[i];
      auto aux2 = z[i]/(1.-beta+beta*aux);
      x[i] = aux2;
      y[i] = aux*aux2; 
    }
    T = exp(X[ncomp]);
    P = exp(X[ncomp+1]);

    if (count == 0)
    {
      file << "\n";
      if (EquiType == Bubble) file << "Bubble line Curve" << std::endl;
      if (EquiType == Dew) file << "Dew line Curve" << std::endl;
      if (EquiType == FractionalBeta) file << "Fractional beta Curve" << std::endl;
      file << "T(K)\t\t\t\t\tP(Pa)\t\t\t\t\tV.Liq(m³/mol)\t\t\tV.Vap(m³/mol)\t\t\t";  
      for (i = 0; i < ncomp; i++) file << "x[" << i << "]" << "\t\t\t\t\t";
      for (i = 0; i < ncomp; i++) file << "y[" << i << "]" << "\t\t\t\t\t"; 
      file << "\n"; 
      compute(props,Tc,Pc,omega,T,P,x,model,kij); 
      auto VLiq = props.V;
      compute(props,Tc,Pc,omega,T,P,y,model,kij); 
      auto VVap = props.V;
      file << T << "\t" << P << "\t" << VLiq << "\t" << VVap;
      for (i = 0; i < ncomp; i++) file << "\t" << x[i];
      for (i = 0; i < ncomp; i++) file << "\t" << y[i];        
      file << "\n";
    }
    else
    {
      compute(props,Tc,Pc,omega,T,P,x,model,kij); 
      auto VLiq = props.V;
      compute(props,Tc,Pc,omega,T,P,y,model,kij); 
      auto VVap = props.V;
      file << T << "\t" << P << "\t" << VLiq << "\t" << VVap;
      for (i = 0; i < ncomp; i++) file << "\t" << x[i];
      for (i = 0; i < ncomp; i++) file << "\t" << y[i];        
      file << "\n";
    }

    // Check if the critical point has been reached
    // Getting lnK values to compute its norm
    for (i = 0; i < ncomp; i++) lnK_aux[i] = X[i];
    maxlnK = norm(lnK_aux,0);
    //CriticalPointReached = (maxlnK < crittol); // dew and bubble as separated curves ///comentado 19 jan 26, regra que para no ponto critio
    CriticalPointReached = (T < 208.15); // complete curve starting from dew  descomentado 19 jan 26
    //CriticalPointReached = (P < 101325.); // complete curve starting from bubble
    

    // Update selection variable based on sensitivities
    SpecifiedVariableUpdated = false;  
    if ((count > 0) && (opts.ChangeSpecfVarBasedOnSensi))
    {
      for (i = 0; i < nvar; i++) maxAbsdXdS[i] = dXdS[i];
      idxNew = idxMaxAbsVec(maxAbsdXdS); 
      if (abs(dXdS[idxNew]) > 1.1*abs(dXdS[idx]))
      { 
        //fac = dXdS[idxNew] < 0. ? -1. : 1. ; 
        SpecifiedVariableUpdated = true;
        idx = idxNew;
        dS *= 0.5;
        SSplus = Xplus[idx];
      }
    }  

    // Update the step
    if (k > targetIter) dS *= 0.5;
    else if (k < targetIter) dS *= 2;

    // Apply upper, lower bound to step size
    if (dS > maxdS) dS = maxdS;
    else if (dS < mindS) dS = mindS;

    // Update Specified variable
    fac = 1.; // dS is positive at the first calculation
    if ((count > 0) && (X[idx] < Xminus[idx])) fac = -1.;
    //if ((count > 0) && (X[idx] < Xminus[idx]) && (idx > (ncomp - 1))) fac = -1.;
    if (maxX <= 1.) dS = abs((dS/X[idx])) < maxX ? dS : abs((maxX * X[idx]));
    SSplus = X[idx]+fac*dS; // specified variable increases or decreases 

    if ((count > 0) && !SpecifiedVariableUpdated && UseCubicInterpolation) 
    {
      // only uses cubic interpolation when specified variable does change
      PolynomialInitialEstimate(InterpolScheme,Xminus,X,Xplus,SSminus,SS,SSplus,dXdS,
        dXdSminus);
    }
    else
    {
      // New initial estimate
      InterpolScheme = "Linear";
      LinearInitialStimate(Xplus,X,dS,dXdS,nvar);
    } 

    // Display message
    if (opts.Verbose)
    { 
      logFile << "Interpolation Scheme: " << InterpolScheme << std::endl; 
      logFile << "Mixed error in Newton's method (should be < 1): " << mix_err << std::endl;
      logFile << "Nº of Newton iterations: " << k << std::endl; 
      logFile << "Delta of specified variable for the next iteration: " << dS << std::endl;
      logFile << "Next Iteration variable updated index: " << idx << std::endl; 
      for (i = 0; i < ncomp; i++) logFile << "lnK[" << i << "]\t\t";
      logFile << "lnT(K)\t\t" << "lnP(Pa)\t\t";
      logFile << std::endl; 
      for (i = 0; i < ncomp; i++) logFile << X[i] << "\t";
      logFile << X[ncomp] << "\t" << X[ncomp+1];
      logFile << std::endl; 
      logFile << "Sensitivity vector: " << std::endl;
      for (i = 0; i < nvar; i++) logFile << "dXdS[" << i << "] = " << dXdS[i] << std::endl; 
      logFile << "\n" << std::endl;        
    }

    // Update variables
    for (i = 0; i < nvar; i++)
    {
      Xminus[i] = X[i];
      X[i] = Xplus[i];
      dXdSminus[i] = dXdS[i];
    }

    SSminus = SS;
    SS = SSplus;

    count++; 
  }

  if (count==MaxEquilIter) logFile << "Max number of line points (" << MaxEquilIter 
    << ") has been reached but critical point was not found";

  logFile.close();  
}

auto PolynomialInitialEstimate(String &InterpolScheme, Vec<real> &X0 , Vec<real> &X1,
  Vec<real> &X2, real SS0, real SS1, real SS2, Vec<real> &dXdS1, 
  Vec<real> &dXdS0) -> void
{
  auto nvar = X0.size();
  InterpolScheme = "Polynomial (Cubic)";
  auto deltaS2S0 = SS2-SS0;
  auto deltaS1S0 = SS1-SS0; 
  for (auto i = 0; i < nvar; i++)
  {
    X2[i] = X0[i] + dXdS0[i]*deltaS2S0 + (3.*(X1[i]-X0[i])-deltaS1S0*(dXdS1[i]+2.*dXdS0[i]))*
      deltaS2S0*deltaS2S0/deltaS1S0/deltaS1S0 + (deltaS1S0*(dXdS1[i]+dXdS0[i])-2.*(X1[i]-X0[i]))*
      deltaS2S0*deltaS2S0*deltaS2S0/deltaS1S0/deltaS1S0/deltaS1S0;   
  }
}

auto LinearInitialStimate(Vec<real> &X, Vec<real> &X0, real dS, Vec<real> &dXdS, int n) 
  -> void
{
  for (auto i = 0; i < n; i++) X[i] = X0[i]+dXdS[i]*dS;
}  

auto IsocurveEquations(Vec<real> &FF, Vec<real> &dFdX, Vec<real> &dFdS, Vec<real> &X,
  real beta, Vec<real> &z, real SS, int idx, param &parameters,
  bool calc_dFdX, int ncomp, bool &convergeDensity, LnPhiDerivativeType LnPhiDerivType) -> void
{
  auto& [props, Tc, Pc, omega, model, kij] = parameters;

  int i,j,k,l;

  using Vector = Eigen::Matrix<real, Eigen::Dynamic, 1>;
  using Matrix = Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>;

  static Vec<real> K(ncomp), x(ncomp), y(ncomp);
  K.resize(ncomp);
  x.resize(ncomp);
  y.resize(ncomp);
  for (i = 0; i < ncomp; i++) K[i] = exp(X[i]);
  real T = exp(X[ncomp]);
  real P = exp(X[ncomp+1]);

  // Mole fractions
  for (i = 0; i < ncomp; i++)
  {
    x[i] = z[i]/(1.-beta+beta*K[i]);
    y[i] = K[i]*x[i];
  }
  
  static Vec<real> Func;

  compute(props,Tc,Pc,omega,T,P,x,model,kij);
  auto lnFug_Liq = props.ln_phi;
  auto lnPhiT_Liq = props.ln_phiT;
  auto lnPhiP_Liq = props.ln_phiP;
  auto lnPhiX_Liq = props.ln_phi_xk;

  Vector dlnFugdT_Liq(ncomp), dlnFugdP_Liq(ncomp);
  Matrix dlnFugdx_Liq(ncomp,ncomp);

  compute(props,Tc,Pc,omega,T,P,y,model,kij);
  auto lnFug_Vap = props.ln_phi;
  auto lnPhiT_Vap = props.ln_phiT;
  auto lnPhiP_Vap = props.ln_phiP;
  auto lnPhiX_Vap = props.ln_phi_xk;

  Vector dlnFugdT_Vap(ncomp), dlnFugdP_Vap(ncomp);
  Matrix dlnFugdx_Vap(ncomp,ncomp);

  if (calc_dFdX)
  {
    switch (LnPhiDerivType)
    {
      case Analytical:
        for (i = 0; i < ncomp; i++)
        {
          dlnFugdT_Liq(i) = lnPhiT_Liq[i];
          dlnFugdP_Liq(i) = lnPhiP_Liq[i];
          dlnFugdT_Vap(i) = lnPhiT_Vap[i];
          dlnFugdP_Vap(i) = lnPhiP_Vap[i];
          for (j = 0; j < ncomp; j++)
          {
            dlnFugdx_Liq(i,j) = lnPhiX_Liq[i][j];
            dlnFugdx_Vap(i,j) = lnPhiX_Vap[i][j];
          }
        }
        break;

      case Automatic:
        dlnFugdT_Liq = jacobian(driver_lnphi, wrt(T), at(props,Tc,Pc,omega,T,P,x,model,kij),Func);
        dlnFugdP_Liq = jacobian(driver_lnphi, wrt(P), at(props,Tc,Pc,omega,T,P,x,model,kij),Func);
        dlnFugdx_Liq = jacobian(driver_lnphi, wrt(x), at(props,Tc,Pc,omega,T,P,x,model,kij),Func);

        dlnFugdT_Vap = jacobian(driver_lnphi, wrt(T), at(props,Tc,Pc,omega,T,P,y,model,kij),Func);
        dlnFugdP_Vap = jacobian(driver_lnphi, wrt(P), at(props,Tc,Pc,omega,T,P,y,model,kij),Func);
        dlnFugdx_Vap = jacobian(driver_lnphi, wrt(y), at(props,Tc,Pc,omega,T,P,y,model,kij),Func);
        break;

      case Numerical:
      {
        auto hT = T*SQRTepsilon;
        auto hP = P*SQRTepsilon;
       // auto hT = SQRTepsilon;
        //auto hP = SQRTepsilon;

        Vec<real> lnphiPlus(ncomp), lnphiMinus(ncomp);

        compute(props,Tc,Pc,omega,T+hT,P,x,model,kij);
        lnphiPlus = props.ln_phi;
        compute(props,Tc,Pc,omega,T-hT,P,x,model,kij);
        lnphiMinus = props.ln_phi;
        for (i = 0; i < ncomp; i++) dlnFugdT_Liq(i) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hT);

        compute(props,Tc,Pc,omega,T,P+hP,x,model,kij);
        lnphiPlus = props.ln_phi;
        compute(props,Tc,Pc,omega,T,P-hP,x,model,kij);
        lnphiMinus = props.ln_phi;
        for (i = 0; i < ncomp; i++) dlnFugdP_Liq(i) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hP);
      for (j = 0; j < ncomp; j++)
        {
          auto hz = x[j]*SQRTepsilon;
          if (hz == 0.) hz = SQRTepsilon;
          Vec<real> x_aux = x;
          x_aux[j] += hz;
          compute(props,Tc,Pc,omega,T,P,x_aux,model,kij);
          lnphiPlus = props.ln_phi;
          x_aux[j] -= 2.*hz;
          compute(props,Tc,Pc,omega,T,P,x_aux,model,kij);
          lnphiMinus = props.ln_phi;
          for (i = 0; i < ncomp; i++) dlnFugdx_Liq(i,j) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hz);
        }

        compute(props,Tc,Pc,omega,T+hT,P,y,model,kij);
        lnphiPlus = props.ln_phi;
        compute(props,Tc,Pc,omega,T-hT,P,y,model,kij);
        lnphiMinus = props.ln_phi;
         for (i = 0; i < ncomp; i++) dlnFugdT_Vap(i) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hT);

        compute(props,Tc,Pc,omega,T,P+hP,y,model,kij);
        lnphiPlus = props.ln_phi;
        compute(props,Tc,Pc,omega,T,P-hP,y,model,kij);
        lnphiMinus = props.ln_phi;
        for (i = 0; i < ncomp; i++) dlnFugdP_Vap(i) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hP);

        for (j = 0; j < ncomp; j++)
        {
          auto hz = y[j]*SQRTepsilon;
          if (hz == 0.) hz = SQRTepsilon;
          Vec<real> y_aux = y;
          y_aux[j] += hz;
          compute(props,Tc,Pc,omega,T,P,y_aux,model,kij);
          lnphiPlus = props.ln_phi;
          y_aux[j] -= 2.*hz;
          compute(props,Tc,Pc,omega,T,P,y_aux,model,kij);
          lnphiMinus = props.ln_phi;
          for (i = 0; i < ncomp; i++) dlnFugdx_Vap(i,j) = (lnphiPlus[i]-lnphiMinus[i])/(2.*hz);
        }

        break;
      }
    }

    // Evaluate isocurve equations


    std::fill(FF.begin(),FF.end(),0.);
    real sum = 0.;
    for (i = 0; i < ncomp; i++)
    {
      FF[i] = X[i] + lnFug_Vap[i] - lnFug_Liq[i];
      sum += y[i] - x[i];
    }
    FF[ncomp] = sum;   // derivatives of this equation wrt T and P are 0
    FF[ncomp+1] = X[idx] - SS;

    static Vec<real> dxdK(ncomp), dydK(ncomp);
    dxdK.resize(ncomp);
    dydK.resize(ncomp);
    std::fill(dFdX.begin(),dFdX.end(),0.);
    
    // Derivatives of mole fractions wrt lnK values
    for (i = 0; i < ncomp; i++)           
    {
      dxdK[i] = -beta*x[i]*x[i]/z[i];
      dydK[i] = x[i]+K[i]*dxdK[i];
    }

    // Evaluate Jacobian
    k = -1;
    for (i = 0; i < ncomp; i++) // i refers to an equation in F
    {
      for (j = 0; j < ncomp; j++) // j refers to component
      {
        auto deltaij = 0.;
        if (i == j) deltaij = 1.;

        // Jacobian wrt lnK
        dFdX[++k] = deltaij + K[j]*(dlnFugdx_Vap(i,j)*dydK[j] - dlnFugdx_Liq(i,j)*dxdK[j]); 
      }

      // Jacobian wrt T      
      dFdX[++k] = T*(dlnFugdT_Vap(i) - dlnFugdT_Liq(i)); 

      // Jacobian wrt P
      dFdX[++k] = P*(dlnFugdP_Vap(i) - dlnFugdP_Liq(i));   
    }

    // Jacobian wrt lnK
    for (i = 0; i < ncomp; i++) dFdX[++k] = K[i]*(dydK[i]-dxdK[i]);
    
    // Jacobian wrt specified variable (lnK , lnT or lnP) 
    dFdX[k+=3+idx] = 1.;
  
    std::fill(dFdS.begin(),dFdS.end(),0.);
    dFdS[ncomp+1] = -1.;

  }

  else
  {
    // Evaluate isocurve equations
    std::fill(FF.begin(),FF.end(),0.);
    real sum = 0.;
    for (i = 0; i < ncomp; i++)
    {
      FF[i] = X[i] + lnFug_Vap[i] - lnFug_Liq[i];
      sum += y[i] - x[i];
    }
    FF[ncomp] = sum;   // derivatives of this equation wrt T and P are 0
    FF[ncomp+1] = X[idx] - SS;
  }

}