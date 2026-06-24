#include "../../ThermoPropCubicEos.hpp"

auto lnphi(CubicEOSProps& props, std::vector<real> &Tcr,
  std::vector<real> &Pcr, std::vector<real> &omega, real T, real P,
  std::vector<real> &z, CubicEOSModel &EoSModel,
  std::vector<std::vector<real>> &kij) -> std::vector<real>
{
  compute(props,Tcr,Pcr,omega,T,P,z,EoSModel,kij);
  return props.ln_phi;
}

auto absolute_error(real approx, real exact) -> real
{
  const auto approx_val = val(approx);
  const auto exact_val = val(exact);
  return std::abs(approx_val - exact_val);
}

int main ()
{
  std::vector<real> temperatures{200.0, 300.0, 800.0};
  std::vector<real> pressures{1.0e5, 1.0e6, 1.0e7};

  auto databasePath = "../../database/SNG1.yml";
  auto components = "C1 C2 iC5";
  std::vector<real> Tc,Pc,omega;
  std::vector<real> z{0.94085, 0.04468, 0.01447};
  auto ncomp = z.size();

  std::vector<std::vector<real>> kij(ncomp, std::vector<real>(ncomp));
  auto EoSModel = CubicEOSModel::PengRobinson;
  CubicEOSProps props;

  read_database(Tc,Pc,omega,databasePath,components);
  if (Tc.size() != ncomp || Pc.size() != ncomp || omega.size() != ncomp)
  {
    std::cout << "Failed to load component data from `" << databasePath << "`." << std::endl;
    return 1;
  }

  String outputFileName = "output.out";
  std::ofstream output_file (outputFileName);
  output_file << std::scientific << std::setprecision(10);
  if (!output_file.is_open())
  {
    std::cout << "Error opening output file: `" << outputFileName << "`." << std::endl;
    return 0;
  }

  std::vector<real> F;
  std::vector<real> lnphiT_numerical(ncomp), lnphiP_numerical(ncomp);
  std::vector<real> lnphi_autodiff_T, lnphi_autodiff_P;
  std::vector<std::vector<real>> lnphiz_autodiff(ncomp, std::vector<real>(ncomp));
  std::vector<std::vector<real>> lnphiz_numerical(ncomp, std::vector<real>(ncomp));
  std::vector<std::vector<real>> lnphiz_error_autodiff(ncomp, std::vector<real>(ncomp));
  std::vector<std::vector<real>> lnphiz_error_numerical(ncomp, std::vector<real>(ncomp));
  std::vector<real> z_aux(ncomp);
  real lnphiPlus, lnphiMinus;

  auto print_vector = [&](const String &label, const std::vector<real> &values)
  {
    output_file << "  " << label << " = [";
    for (auto i = 0U; i < values.size(); i++)
    {
      if (i > 0)
        output_file << ", ";
      output_file << values[i];
    }
    output_file << "]" << std::endl;
  };

  auto print_matrix = [&](const String &label, const std::vector<std::vector<real>> &matrix)
  {
    output_file << "  " << label << " =" << std::endl;
    for (const auto &row : matrix)
    {
      output_file << "    [";
      for (auto j = 0U; j < row.size(); j++)
      {
        if (j > 0)
          output_file << ", ";
        output_file << row[j];
      }
      output_file << "]" << std::endl;
    }
  };

  for(auto T : temperatures)
  {
    for(auto P : pressures)
    {
      compute(props,Tc,Pc,omega,T,P,z,EoSModel,kij);

      auto lnphiT_ad = jacobian(lnphi, wrt(T), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij), F);
      lnphi_autodiff_T.assign(lnphiT_ad.data(), lnphiT_ad.data() + ncomp);

      auto lnphiP_ad = jacobian(lnphi, wrt(P), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij), F);
      lnphi_autodiff_P.assign(lnphiP_ad.data(), lnphiP_ad.data() + ncomp);

      // Numerical derivatives with respect to T
      auto hT = T*SQRTepsilon;
      for(auto i = 0; i < ncomp; i++)
      {
        compute(props,Tc,Pc,omega,T+hT,P,z,EoSModel,kij);
        lnphiPlus = props.ln_phi[i];
        compute(props,Tc,Pc,omega,T-hT,P,z,EoSModel,kij);
        lnphiMinus = props.ln_phi[i];
        lnphiT_numerical[i] = (lnphiPlus - lnphiMinus) / (2.*hT);
      }

      // Numerical derivatives with respect to P
      auto hP = P*SQRTepsilon;
      for(auto i = 0; i < ncomp; i++)
      {
        compute(props,Tc,Pc,omega,T,P+hP,z,EoSModel,kij);
        lnphiPlus = props.ln_phi[i];
        compute(props,Tc,Pc,omega,T,P-hP,z,EoSModel,kij);
        lnphiMinus = props.ln_phi[i];
        lnphiP_numerical[i] = (lnphiPlus - lnphiMinus) / (2.*hP);
      }

      // Autodiff and numerical derivatives with respect to composition
      for(auto j = 0; j < ncomp; j++)
      {
        auto lnphiz_ad = jacobian(lnphi, wrt(z[j]), at(props,Tc,Pc,omega,T,P,z,EoSModel,kij), F);
        for(auto i = 0; i < ncomp; i++)
        {
          lnphiz_autodiff[i][j] = lnphiz_ad(i);
        }

        auto hz = z[j]*SQRTepsilon;
        for(auto i = 0; i < ncomp; i++)
        {
          z_aux = z;
          z_aux[j] += hz;
          compute(props,Tc,Pc,omega,T,P,z_aux,EoSModel,kij);
          lnphiPlus = props.ln_phi[i];
          z_aux = z;
          z_aux[j] -= hz;
          compute(props,Tc,Pc,omega,T,P,z_aux,EoSModel,kij);
          lnphiMinus = props.ln_phi[i];
          lnphiz_numerical[i][j] = (lnphiPlus - lnphiMinus) / (2.*hz);
          lnphiz_error_autodiff[i][j] = absolute_error(lnphiz_autodiff[i][j], props.ln_phi_xk[i][j]);
          lnphiz_error_numerical[i][j] = absolute_error(lnphiz_numerical[i][j], props.ln_phi_xk[i][j]);
        }
      }

      output_file << "Caso: T = " << T << " K, P = " << P << " Pa" << std::endl;
      std::vector<real> lnphiT_error_autodiff(ncomp), lnphiT_error_numerical(ncomp);
      std::vector<real> lnphiP_error_autodiff(ncomp), lnphiP_error_numerical(ncomp);

      for (auto i = 0; i < ncomp; i++)
      {
        lnphiT_error_autodiff[i] = absolute_error(lnphi_autodiff_T[i], props.ln_phiT[i]);
        lnphiT_error_numerical[i] = absolute_error(lnphiT_numerical[i], props.ln_phiT[i]);
        lnphiP_error_autodiff[i] = absolute_error(lnphi_autodiff_P[i], props.ln_phiP[i]);
        lnphiP_error_numerical[i] = absolute_error(lnphiP_numerical[i], props.ln_phiP[i]);
      }

      print_vector("ln_phi", props.ln_phi);

      print_vector("dlnphi/dT analitico", props.ln_phiT);
      print_vector("dlnphi/dT autodiff", lnphi_autodiff_T);
      print_vector("dlnphi/dT numerico", lnphiT_numerical);
      print_vector("erro absoluto dlnphi/dT autodiff", lnphiT_error_autodiff);
      print_vector("erro absoluto dlnphi/dT numerico", lnphiT_error_numerical);

      print_vector("dlnphi/dP analitico", props.ln_phiP);
      print_vector("dlnphi/dP autodiff", lnphi_autodiff_P);
      print_vector("dlnphi/dP numerico", lnphiP_numerical);
      print_vector("erro absoluto dlnphi/dP autodiff", lnphiP_error_autodiff);
      print_vector("erro absoluto dlnphi/dP numerico", lnphiP_error_numerical);

      print_matrix("dlnphi/dx analitico", props.ln_phi_xk);
      print_matrix("dlnphi/dx autodiff", lnphiz_autodiff);
      print_matrix("dlnphi/dx numerico", lnphiz_numerical);
      print_matrix("erro absoluto dlnphi/dx autodiff", lnphiz_error_autodiff);
      print_matrix("erro absoluto dlnphi/dx numerico", lnphiz_error_numerical);
      output_file << std::endl;
    }
  }

  std::cout << "Resultados escritos em `" << outputFileName << "`." << std::endl;

  return 0;
}
