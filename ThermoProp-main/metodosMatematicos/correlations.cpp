#include "correlations.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

// ============================================================
// Internal fitted correlations for heavy-fraction characterization.
// These are internal empirical fits inspired by SCN/Whitson-like trends.
// They are NOT literal implementations of the classical Lee-Kesler
// or Riazi-Daubert published correlations.
//
// The Whitson/Riazi-Daubert route uses:
//   - Watson K factor (Whitson 1983 Eq. 15 and 17)
//   - SG/Tb closure via K
//   - Riazi-Daubert for Tc/Pc (Whitson 1983 Table 1, SI units)
//   - Edmister for omega (Whitson 1983 Eq. 21 ref)
//
// The AutoSmoothFromMM policy is a heuristic internal to this software.
// It is NOT a strict Whitson-1983 implementation.
// ============================================================

// ============================================================
// Watson K from MM and SG using Whitson Eq. 17:
// K = 4.5579 * MM^0.15178 * SG^-0.84573
// MM must be in g/mol.
// ============================================================
auto estimateWatsonKFromMMSG(double MM_gmol, double SG) -> double
{
  return 4.5579 * std::pow(MM_gmol, 0.15178) / std::pow(SG, 0.84573);
}

// ============================================================
// SG from MM and Watson K (inverse of Whitson Eq. 17):
// SG = [ (4.5579 * MM^0.15178) / K ] ^ (1 / 0.84573)
// MM must be in g/mol.
// ============================================================
auto estimateSGFromMMAndWatsonK(double MM_gmol, double K) -> double
{
  double base = (4.5579 * std::pow(MM_gmol, 0.15178)) / K;
  return std::pow(base, 1.0 / 0.84573);
}

// ============================================================
// Tb from SG and Watson K using Whitson Eq. 15 (SI version):
// K = 1.21644 * Tb(K)^(1/3) / SG
// Tb(K) = [ (K * SG) / 1.21644 ]^3
// Returns Tb in Kelvin.
// ============================================================
auto estimateTbFromSGAndWatsonK(double SG, double K) -> double
{
  double base = (K * SG) / 1.21644;
  return std::pow(base, 3.0);
}

// ============================================================
// Heuristic Watson K as a smooth sigmoidal function of MM.
// This is an internal policy, NOT from Whitson (1983).
//
// K(MM) = K_min + (K_max - K_min) / (1 + exp(a * (MM - M0)))
//
// Parameters:
//   K_min = 10.5  (heavy pseudocomponents)
//   K_max = 12.2  (lighter pseudocomponents)
//   a     = 0.015
//   M0    = 200.0
//
// MM must be in g/mol. Result is clamped to [10.0, 13.0].
// ============================================================
auto estimateWatsonKAutoSmooth(double MM_gmol) -> double
{
  constexpr double K_min = 10.5;
  constexpr double K_max = 12.2;
  constexpr double a = 0.015;
  constexpr double M0 = 200.0;

  const double x = a * (MM_gmol - M0);
  double K = K_min + (K_max - K_min) / (1.0 + std::exp(x));

  // Safety clamp
  return std::clamp(K, 10.0, 13.0);
}

// ============================================================
// SG from Tb and MM by equating Whitson Eq. 15 and Eq. 17:
// SG = [ (1.21644 * Tb^(1/3)) / (4.5579 * MM^0.15178) ] ^ (1 / 0.15427)
// ============================================================
auto estimateSGFromTbAndMM(double Tb_K, double MM_gmol) -> double
{
  double num = 1.21644 * std::cbrt(Tb_K);
  double den = 4.5579 * std::pow(MM_gmol, 0.15178);
  return std::pow(num / den, 1.0 / 0.15427);
}

// ============================================================
// Tc from Tb and SG using Riazi-Daubert (Whitson 1983, Table 1, SI units):
// Tc(K) = 19.0623 * Tb^0.58848 * SG^0.3596
// ============================================================
auto estimateTcFromTbSG_RiaziDaubert(double Tb_K, double SG) -> double
{
  return 19.0623 * std::pow(Tb_K, 0.58848) * std::pow(SG, 0.3596);
}

// ============================================================
// Pc from Tb and SG using Riazi-Daubert (Whitson 1983, Table 1, SI units):
// If Tb <= 730 K:
//   Pc(kPa) = 5.53028e9 * Tb^-2.3125 * SG^2.3201
// If Tb > 730 K:
//   Pc(kPa) = 1.71589e14 * Tb^-3.86618 * SG^4.2448
// Result is converted from kPa to Pa.
// ============================================================
auto estimatePcFromTbSG_RiaziDaubert(double Tb_K, double SG) -> double
{
  double Pc_kPa = 0.0;
  if (Tb_K <= 730.0) {
    Pc_kPa = 5.53028e9 * std::pow(Tb_K, -2.3125) * std::pow(SG, 2.3201);
  } else {
    Pc_kPa = 1.71589e14 * std::pow(Tb_K, -3.86618) * std::pow(SG, 4.2448);
  }
  return Pc_kPa * 1000.0; // kPa to Pa
}

// ============================================================
// Omega using Edmister equation (Whitson 1983, Eq 21 ref):
// omega = (3/7) * [ log10(Pc / Patm) / (Tc/Tb - 1) ] - 1
// where Patm = 101325 Pa
// ============================================================
auto estimateOmegaFromTcPcTb_Edmister(double Tc_K, double Pc_Pa, double Tb_K) -> double
{
  const double P_atm = 101325.0; // Pa
  double theta = Tc_K / Tb_K;
  if (theta <= 1.0) theta = 1.01; // Avoid division by zero or negative
  double log_term = std::log10(Pc_Pa / P_atm);
  return (3.0 / 7.0) * (log_term / (theta - 1.0)) - 1.0;
}

// ============================================================
// Deprecated — API compatibility wrappers
// ============================================================

auto estimateTbRiaziDaubert(double MM) -> double
{
  // Internal fit: routes to Whitson-based SG + Tb estimation
  double K = estimateWatsonKAutoSmooth(MM);
  double SG = estimateSGFromMMAndWatsonK(MM, K);
  return estimateTbFromSGAndWatsonK(SG, K);
}

auto estimateSGRiaziDaubert(double MM) -> double
{
  // Internal fit: routes to Whitson-based SG estimation
  double K = estimateWatsonKAutoSmooth(MM);
  return estimateSGFromMMAndWatsonK(MM, K);
}

auto estimateTcRiaziDaubert(double MM, double SG) -> double
{
  double K = estimateWatsonKFromMMSG(MM, SG);
  double Tb = estimateTbFromSGAndWatsonK(SG, K);
  return estimateTcFromTbSG_RiaziDaubert(Tb, SG);
}

auto estimatePcRiaziDaubert(double MM, double SG) -> double
{
  double K = estimateWatsonKFromMMSG(MM, SG);
  double Tb = estimateTbFromSGAndWatsonK(SG, K);
  return estimatePcFromTbSG_RiaziDaubert(Tb, SG);
}

auto estimateOmegaRiaziDaubert(double Tc, double Pc, double Tb) -> double
{
  return estimateOmegaFromTcPcTb_Edmister(Tc, Pc, Tb);
}

// ============================================================
// Validation helpers
// ============================================================
static auto isValidSG(double SG) -> bool
{
  return std::isfinite(SG) && SG > 0.0 && SG <= 2.0;
}

static auto isValidTb(double Tb) -> bool
{
  return std::isfinite(Tb) && Tb > 150.0 && Tb <= 1500.0;
}

static auto resetStatusFields(Component &comp) -> void
{
  comp.provenance = CharacterizationProvenance::Unknown;
  comp.characterization_ok = false;
  comp.characterization_message.clear();
}

static auto failStatus(Component &comp, CharacterizationStatus status, const char *msg)
    -> CharacterizationStatus
{
  comp.characterization_message = msg;
  comp.characterization_ok = false;
  return status;
}

// ============================================================
// Internal: Whitson/Riazi-Daubert characterization route
// ============================================================
static auto characterizeWhitsonRiaziDaubert(Component &comp,
                                            const CharacterizationOptions &opts)
    -> CharacterizationStatus
{
  resetStatusFields(comp);

  const bool has_sg = isValidSG(comp.SG);
  const bool has_tb = isValidTb(comp.Tb);

  double K = std::numeric_limits<double>::quiet_NaN();

  // Case A — both SG and Tb provided
  if (has_sg && has_tb) {
    comp.provenance = CharacterizationProvenance::ProvidedTbAndSG;
  }
  // Case B — SG provided, Tb absent
  else if (has_sg && !has_tb) {
    if (opts.watson_k_policy == WatsonKPolicy::Provided &&
        std::isfinite(opts.provided_watson_k)) {
      K = opts.provided_watson_k;
      comp.provenance = CharacterizationProvenance::UserProvidedWatsonK;
    } else {
      K = estimateWatsonKFromMMSG(comp.MM, comp.SG);
      comp.provenance = CharacterizationProvenance::ProvidedSGEstimatedTb;
    }
    comp.Tb = estimateTbFromSGAndWatsonK(comp.SG, K);
  }
  // Case C — Tb provided, SG absent
  else if (!has_sg && has_tb) {
    comp.SG = estimateSGFromTbAndMM(comp.Tb, comp.MM);
    comp.provenance = CharacterizationProvenance::ProvidedTbEstimatedSG;
  }
  // Case D/E — neither provided
  else {
    switch (opts.watson_k_policy) {
      case WatsonKPolicy::CalibrateFromHeavyFractionSG:
        if (!std::isfinite(opts.heavy_fraction_sg)) {
          return failStatus(comp, CharacterizationStatus::HeavyFractionSGRequired,
                            "Heavy-fraction SG required for Watson K calibration.");
        }
        return failStatus(comp, CharacterizationStatus::NotImplemented,
                          "Heavy-fraction SG calibration not yet implemented.");

      case WatsonKPolicy::Provided:
        if (!std::isfinite(opts.provided_watson_k)) {
          return failStatus(comp, CharacterizationStatus::NumericalFailure,
                            "Provided Watson K is invalid.");
        }
        K = opts.provided_watson_k;
        comp.provenance = CharacterizationProvenance::UserProvidedWatsonK;
        break;

      case WatsonKPolicy::AutoSmoothFromMM:
        // Heuristic Whitson-inspired characterization using internal Watson-K policy
        // when heavy-fraction SG is unavailable.
        K = estimateWatsonKAutoSmooth(comp.MM);
        comp.provenance = CharacterizationProvenance::AutoWatsonKFromMM;
        break;
    }

    comp.SG = estimateSGFromMMAndWatsonK(comp.MM, K);
    comp.Tb = estimateTbFromSGAndWatsonK(comp.SG, K);
  }

  // Compute critical properties
  comp.Tc = estimateTcFromTbSG_RiaziDaubert(comp.Tb, comp.SG);
  comp.Pc = estimatePcFromTbSG_RiaziDaubert(comp.Tb, comp.SG);
  comp.omega = estimateOmegaFromTcPcTb_Edmister(comp.Tc, comp.Pc, comp.Tb);

  // Validate results — no silent fallbacks
  if (!isValidSG(comp.SG))
    return failStatus(comp, CharacterizationStatus::InvalidSG, "Invalid SG after characterization.");
  if (!isValidTb(comp.Tb))
    return failStatus(comp, CharacterizationStatus::InvalidTb, "Invalid Tb after characterization.");
  if (!std::isfinite(comp.Tc) || comp.Tc <= comp.Tb || comp.Tc > 2000.0)
    return failStatus(comp, CharacterizationStatus::InvalidTc, "Invalid Tc after characterization.");
  if (!std::isfinite(comp.Pc) || comp.Pc <= 1.0e4 || comp.Pc > 1.0e8)
    return failStatus(comp, CharacterizationStatus::InvalidPc, "Invalid Pc after characterization.");
  if (!std::isfinite(comp.omega) || comp.omega < -0.5 || comp.omega > 2.0)
    return failStatus(comp, CharacterizationStatus::InvalidOmega, "Invalid omega after characterization.");

  comp.characterization_ok = true;
  comp.characterization_message =
      (comp.provenance == CharacterizationProvenance::AutoWatsonKFromMM)
          ? "Heuristic Whitson-inspired characterization using automatic Watson K from MM."
          : "Characterization completed successfully.";

  return CharacterizationStatus::Ok;
}

// ============================================================
// Internal: Alkane MM characterization route
// Uses direct empirical fits from MM only — no Watson K, no SG, no Tb path.
// Equations from documentation.tex (lines 217-236):
//   Tb = 60.806 * MM^0.4046
//   SG = 0.4671 * MM^0.1040
//   omega = 2*SG - 1.2
//   Tc = piecewise function of MM
//   Pc = max(1, min(5, 3.5 - 0.003*MM)) * 1e6
// ============================================================
static auto estimateTbInternalAlkaneMM(double MM) -> double
{
  return 60.806 * std::pow(MM, 0.4046);
}

static auto estimateSGInternalAlkaneMM(double MM) -> double
{
  return 0.4671 * std::pow(MM, 0.1040);
}

static auto estimateOmegaInternalAlkaneMM(double SG) -> double
{
  return 2.0 * SG - 1.2;
}

static auto estimateTcInternalAlkaneMM(double MM) -> double
{
  if (MM < 70.0) return 400.0;
  if (MM < 150.0) return 400.0 + 2.0 * (MM - 70.0);
  if (MM < 300.0) return 560.0 + 1.2 * (MM - 150.0);
  return 740.0 + 0.3 * (MM - 300.0);
}

static auto estimatePcInternalAlkaneMM(double MM) -> double
{
  double val_MPa = std::max(1.0, std::min(5.0, 3.5 - 0.003 * MM));
  return val_MPa * 1.0e6;
}

static auto characterizeInternalAlkane(Component &comp,
                                       const CharacterizationOptions &/*opts*/)
    -> CharacterizationStatus
{
  resetStatusFields(comp);

  comp.SG = estimateSGInternalAlkaneMM(comp.MM);
  comp.Tb = estimateTbInternalAlkaneMM(comp.MM);
  comp.Tc = estimateTcInternalAlkaneMM(comp.MM);
  comp.Pc = estimatePcInternalAlkaneMM(comp.MM);
  comp.omega = estimateOmegaInternalAlkaneMM(comp.SG);

  comp.provenance = CharacterizationProvenance::InternalAlkaneMM;

  // Validate results
  if (!isValidSG(comp.SG))
    return failStatus(comp, CharacterizationStatus::InvalidSG, "Invalid SG in InternalAlkaneMM.");
  if (!isValidTb(comp.Tb))
    return failStatus(comp, CharacterizationStatus::InvalidTb, "Invalid Tb in InternalAlkaneMM.");
  if (!std::isfinite(comp.Tc) || comp.Tc <= comp.Tb || comp.Tc > 2000.0)
    return failStatus(comp, CharacterizationStatus::InvalidTc, "Invalid Tc in InternalAlkaneMM.");
  if (!std::isfinite(comp.Pc) || comp.Pc <= 1.0e4 || comp.Pc > 1.0e8)
    return failStatus(comp, CharacterizationStatus::InvalidPc, "Invalid Pc in InternalAlkaneMM.");
  if (!std::isfinite(comp.omega) || comp.omega < -0.5 || comp.omega > 2.0)
    return failStatus(comp, CharacterizationStatus::InvalidOmega, "Invalid omega in InternalAlkaneMM.");

  comp.characterization_ok = true;
  comp.characterization_message = "InternalAlkaneMM characterization completed.";

  return CharacterizationStatus::Ok;
}

// ============================================================
// High-level characterization — orchestrator
// ============================================================
auto characterizeComponent(Component &comp, const CharacterizationOptions &opts)
    -> CharacterizationStatus
{
  resetStatusFields(comp);

  // Validate MM
  if (!std::isfinite(comp.MM) || comp.MM <= 0.0)
    return failStatus(comp, CharacterizationStatus::InvalidMM, "Invalid molar mass.");

  // Detect likely unit error: MM < 1.0 suggests kg/mol instead of g/mol
  if (comp.MM < 1.0)
    return failStatus(comp, CharacterizationStatus::UnitSuspectedForMM,
                      "MM appears to be in kg/mol; expected g/mol.");

  switch (opts.method) {
    case CorrMethod::WhitsonRiaziDaubert:
      return characterizeWhitsonRiaziDaubert(comp, opts);

    case CorrMethod::InternalAlkaneMM:
      return characterizeInternalAlkane(comp, opts);
  }

  return failStatus(comp, CharacterizationStatus::NotImplemented, "Unknown CorrMethod.");
}

// ============================================================
// Convenience overload — legacy signature
// ============================================================
auto characterizeComponent(Component &comp, CorrMethod method) -> CharacterizationStatus
{
  CharacterizationOptions opts;
  opts.method = method;
  opts.watson_k_policy = WatsonKPolicy::AutoSmoothFromMM;
  return characterizeComponent(comp, opts);
}
