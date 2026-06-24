#pragma once

#include "cubic_eos_types.hpp"

// ============================================================
// Correlation functions for heavy-fraction characterization.
//
// This file provides two distinct characterization routes:
//
//   1. Whitson / Riazi-Daubert / Edmister route (CorrMethod::WhitsonRiaziDaubert)
//      — Uses Watson K factor to bridge MM → SG → Tb, then
//        Riazi-Daubert for Tc/Pc and Edmister for omega.
//      — When heavy-fraction SG is unavailable, a heuristic
//        internal Watson-K policy (AutoSmoothFromMM) is used.
//        This is NOT a strict Whitson-1983 implementation.
//
//   2. Internal Alkane MM route (CorrMethod::InternalAlkaneMM)
//      — Direct empirical fits from MM only.
//
// Legacy function names (LeeKesler, RiaziDaubert) are deprecated
// and retained only for API compatibility. They do NOT implement
// the classical Lee-Kesler or Riazi-Daubert published correlations.
// ============================================================

// ============================================================
// Watson K factor functions (Whitson 1983 basis)
// ============================================================

/// Estimate Watson K from MM and SG using Whitson Eq. 17.
/// MM must be in g/mol. SG is dimensionless (60F/60F).
auto estimateWatsonKFromMMSG(double MM_gmol, double SG) -> double;

/// Estimate SG from MM and Watson K (inverse of Whitson Eq. 17).
/// MM must be in g/mol.
auto estimateSGFromMMAndWatsonK(double MM_gmol, double K) -> double;

/// Estimate normal boiling point Tb from SG and Watson K.
/// Uses Whitson Eq. 15 (SI version): K = 1.21644 * Tb^(1/3) / SG.
/// Returns Tb in Kelvin.
auto estimateTbFromSGAndWatsonK(double SG, double K) -> double;

/// Heuristic Watson K as a smooth sigmoidal function of MM.
/// This is an internal policy, NOT from Whitson (1983).
/// K(MM) = K_min + (K_max - K_min) / (1 + exp(a * (MM - M0)))
/// with K_min=10.5, K_max=12.2, a=0.015, M0=200.
/// Result is clamped to [10.0, 13.0].
auto estimateWatsonKAutoSmooth(double MM_gmol) -> double;

/// Estimate SG from Tb and MM by equating Whitson Eq. 15 and Eq. 17.
auto estimateSGFromTbAndMM(double Tb_K, double MM_gmol) -> double;

// ============================================================
// Critical property functions (Riazi-Daubert / Edmister)
// ============================================================

/// Estimate Tc from Tb and SG using Riazi-Daubert (Whitson 1983, Table 1).
/// Tc(K) = 19.0623 * Tb^0.58848 * SG^0.3596
auto estimateTcFromTbSG_RiaziDaubert(double Tb_K, double SG) -> double;

/// Estimate Pc from Tb and SG using Riazi-Daubert (Whitson 1983, Table 1).
/// Uses piecewise formula: Tb <= 730 K and Tb > 730 K.
/// Returns Pc in Pa.
auto estimatePcFromTbSG_RiaziDaubert(double Tb_K, double SG) -> double;

/// Estimate omega from Tc, Pc, Tb using Edmister equation.
/// omega = (3/7) * [ log10(Pc / Patm) / (Tc/Tb - 1) ] - 1
auto estimateOmegaFromTcPcTb_Edmister(double Tc_K, double Pc_Pa, double Tb_K) -> double;

// ============================================================
// Deprecated — API compatibility wrappers only.
// These do NOT implement Lee-Kesler or classical Riazi-Daubert.
// They route to the Whitson / Riazi-Daubert / Edmister functions.
// ============================================================

/// Deprecated: historical API name only. Routes to estimateTcFromTbSG_RiaziDaubert.
static inline auto estimateTcLeeKesler(double Tb, double SG) -> double
{ return estimateTcFromTbSG_RiaziDaubert(Tb, SG); }

/// Deprecated: historical API name only. Routes to estimatePcFromTbSG_RiaziDaubert.
static inline auto estimatePcLeeKesler(double Tb, double SG) -> double
{ return estimatePcFromTbSG_RiaziDaubert(Tb, SG); }

/// Deprecated: historical API name only. Routes to estimateOmegaFromTcPcTb_Edmister.
static inline auto estimateOmegaLeeKesler(double Tb, double SG, double Tc) -> double
{
  double Pc = estimatePcFromTbSG_RiaziDaubert(Tb, SG);
  return estimateOmegaFromTcPcTb_Edmister(Tc, Pc, Tb);
}

/// Deprecated: historical API name only. Routes to internal fits.
auto estimateTbRiaziDaubert(double MM) -> double;

/// Deprecated: historical API name only. Routes to internal fits.
auto estimateSGRiaziDaubert(double MM) -> double;

/// Deprecated: historical API name only. Routes to Whitson-based functions.
auto estimateTcRiaziDaubert(double MM, double SG) -> double;

/// Deprecated: historical API name only. Routes to Whitson-based functions.
auto estimatePcRiaziDaubert(double MM, double SG) -> double;

/// Deprecated: historical API name only. Routes to Edmister.
auto estimateOmegaRiaziDaubert(double Tc, double Pc, double Tb) -> double;

// ============================================================
// High-level characterization: fills Tc, Pc, omega for a component.
// CorrMethod selects the characterization route.
// WatsonKPolicy controls how K is determined when SG is unavailable.
// Returns CharacterizationStatus — caller decides on fallback behavior.
// ============================================================
auto characterizeComponent(Component &comp, const CharacterizationOptions &opts)
    -> CharacterizationStatus;

/// Convenience overload: legacy signature, routes to the new API with defaults.
auto characterizeComponent(Component &comp, CorrMethod method) -> CharacterizationStatus;
