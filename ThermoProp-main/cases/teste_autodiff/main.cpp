#include "../../ThermoPropCubicEos.hpp"

auto lnphi(CubicEOSProps& props, std::vector<real> &Tcr, 
  std::vector<real> &Pcr, std::vector<real> &omega, real T, real P, 
  std::vector<real> &z, CubicEOSModel &EoSModel, 
  std::vector<std::vector<real>> &kij) -> std::vector<real>
{
  compute(props,Tcr,Pcr,omega,T,P,z,EoSModel,kij);
  return props.ln_phi;
}  

int main ()
{
  real T = 323.15;
  real P = 30000000.;
  //auto databasePath = "/home/fellipe/Dropbox/Pessoais/Projetos_Pesquisa/modelagem_termodinamica_reservatorio_petroleo/ThermoPropCubicEos/database/test.yml";
auto databasePath = "/home/paloma/Downloads/fellipe/ThermoProp-main/database/test.yml";
  auto components = "CO2 H2O";
  std::vector<real> Tc,Pc,omega;
  std::vector<std::vector<real>> kij(2, std::vector<real>(2)); //  there are three components, this matrix is initialized as a 3x3 matrix with zeroes
  auto EoSModel = CubicEOSModel::PengRobinson;
  CubicEOSProps props;
  std::vector<real> z{0.5, 0.5};
  read_database(Tc,Pc,omega,databasePath,components);
  auto ncomp = z.size();

  std::vector<real> F;
  real h;
  real lnphiPlus, lnphiMinus;

  // testing derivatives of lnphi wrt T
  auto lnphiT_autodiff = jacobian(lnphi, wrt(T), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij),F);
  std::vector<real> lnphiT_numerical(ncomp); 
  h = T*SQRTepsilon;
  for(auto i = 0; i < ncomp; i++)
  {
    compute(props,Tc,Pc,omega,T+h,P,z,EoSModel,kij);
    lnphiPlus = props.ln_phi[i];
    compute(props,Tc,Pc,omega,T-h,P,z,EoSModel,kij);
    lnphiMinus = props.ln_phi[i];   
    lnphiT_numerical[i] = (lnphiPlus - lnphiMinus) / (2.*h);
  }
  std::cout << "Error in Temperature is: \n";
  for(auto i = 0; i < ncomp; i++) std::cout << "ErrorT[" << i << "] = " << lnphiT_numerical[i] / lnphiT_autodiff(i) - 1. << std::endl;
  // End testing derivatives of lnphi wrt T

  // testing derivatives of lnphi wrt P
  auto lnphiP_autodiff = jacobian(lnphi, wrt(P), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij),F);
  std::vector<real> lnphiP_numerical(ncomp); 
  h = P*SQRTepsilon;
  for(auto i = 0; i < ncomp; i++)
  {
    compute(props,Tc,Pc,omega,T,P+h,z,EoSModel,kij);
    lnphiPlus = props.ln_phi[i];
    compute(props,Tc,Pc,omega,T,P-h,z,EoSModel,kij);
    lnphiMinus = props.ln_phi[i];   
    lnphiP_numerical[i] = (lnphiPlus - lnphiMinus) / (2.*h);
  }
  std::cout << "\nError in Pressure is: \n";
  for(auto i = 0; i < ncomp; i++) std::cout << "ErrorP[" << i << "] = " << lnphiP_numerical[i] / lnphiP_autodiff(i) - 1. << std::endl;
  // End testing derivatives of lnphi wrt P

  // testing derivatives of lnphi wrt z
  std::vector<real> lnphiz_numerical(ncomp),z_aux(ncomp); 
  for(auto j = 0; j < ncomp; j++)
  {  
    auto lnphiz_autodiff = jacobian(lnphi, wrt(z[j]), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij),F);
    h = z[j]*SQRTepsilon;
    for(auto i = 0; i < ncomp; i++)
    {
      z_aux = z;
      z_aux[j] += h; 
      compute(props,Tc,Pc,omega,T,P,z_aux,EoSModel,kij);
      lnphiPlus = props.ln_phi[i];
      z_aux = z;
      z_aux[j] -= h;       
      compute(props,Tc,Pc,omega,T,P,z_aux,EoSModel,kij);
      lnphiMinus = props.ln_phi[i];   
      lnphiz_numerical[i] = (lnphiPlus - lnphiMinus) / (2.*h);
    }
    std::cout << "\nError in z[" << j << "] is: \n";
    for(auto i = 0; i < ncomp; i++) std::cout << "Error_z[" << j << "," << i << "] = " << lnphiz_numerical[i] / lnphiz_autodiff(i) - 1. << std::endl;
  }
  // End testing derivatives of lnphi wrt z

  

  return 0;
}



