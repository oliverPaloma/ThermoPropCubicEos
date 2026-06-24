#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>

// ============================================================
// Type aliases — keep them short and explicit
// ============================================================
template<typename T>
using Vec = std::vector<T>;

using String = std::string;

// ============================================================
// Universal gas constant (J / mol / K)
// ============================================================
static const double R_cubic = 8.314462618;

// ============================================================
// Characterization provenance — how properties were obtained
// (Ciclo 3 — Plano_correcoes_3.md)
// ============================================================
enum class CharacterizationProvenance {
  Unknown = 0,
  HeavyFractionSGCalibrated,
  UserProvidedWatsonK,
  AutoWatsonKFromMM,
  InternalAlkaneMM,
  ProvidedTbAndSG,
  ProvidedSGEstimatedTb,
  ProvidedTbEstimatedSG
};

// ============================================================
// Characterization status — result of characterizeComponent
// (Ciclo 3 — Plano_correcoes_3.md)
// ============================================================
enum class CharacterizationStatus {
  Ok = 0,
  InvalidMM,
  InvalidSG,
  InvalidTb,
  InvalidTc,
  InvalidPc,
  InvalidOmega,
  UnitSuspectedForMM,
  HeavyFractionSGRequired,
  NotImplemented,
  NumericalFailure
};

// ============================================================
// Cubic EOS selection
// ============================================================
enum class CubicEOS { PengRobinson, SRK };

inline String eosName(CubicEOS eos)
{
  return (eos == CubicEOS::PengRobinson) ? "Peng-Robinson" : "SRK";
}

// ============================================================
// Correlation selection
// CorrMethod must have real effect — controls the characterization route.
// (Ciclo 3 — Plano_correcoes_3.md)
// ============================================================
enum class CorrMethod {
  WhitsonRiaziDaubert,
  InternalAlkaneMM
};

inline String corrName(CorrMethod m)
{
  switch (m) {
    case CorrMethod::WhitsonRiaziDaubert: return "Whitson-RiaziDaubert";
    case CorrMethod::InternalAlkaneMM:   return "Internal-Alkane-MM";
  }
  return "Unknown";
}

// ============================================================
// Watson K policy — how Watson K factor is determined
// (Ciclo 3 — Plano_correcoes_3.md)
// ============================================================
enum class WatsonKPolicy {
  CalibrateFromHeavyFractionSG,
  Provided,
  AutoSmoothFromMM
};

// ============================================================
// Characterization options — passed to characterizeComponent
// (Ciclo 3 — Plano_correcoes_3.md)
// ============================================================
struct CharacterizationOptions {
  CorrMethod method = CorrMethod::WhitsonRiaziDaubert;
  WatsonKPolicy watson_k_policy = WatsonKPolicy::AutoSmoothFromMM;
  double provided_watson_k = std::numeric_limits<double>::quiet_NaN();
  double heavy_fraction_sg = std::numeric_limits<double>::quiet_NaN();
  bool strict_validation = true;
};

// ============================================================
// Component description
// ============================================================
struct Component
{
  String name;       // Component name (e.g. "C7", "Pseudo1")
  double z;          // Molar fraction
  double MM;         // Molar mass (g/mol) — MUST be in g/mol, not kg/mol
  int  flag;         // 0 = well-known, 1 = SCN fraction to be lumped
  double Tb;         // Normal boiling point (K) — estimated or from input
  double SG;         // Specific gravity (60°F/60°F) — estimated or from input
  double Tc;         // Critical temperature (K)
  double Pc;         // Critical pressure (Pa)
  double omega;      // Acentric factor
  double a_c;        // Cubic EOS parameter a at Tc (Pa·m⁶/mol²)
  double b;          // Cubic EOS parameter b (m³/mol)

  // Characterization provenance and status (Ciclo 3)
  CharacterizationProvenance provenance = CharacterizationProvenance::Unknown;
  bool characterization_ok = false;
  String characterization_message;
};

// ============================================================
// Mixture after recharacterization
// ============================================================
struct Mixture
{
  Vec<Component> comps;       // All components (light + pseudos)
  Vec<Vec<double>> kij;       // Binary interaction matrix (n x n)
  int n_light;                // Number of well-known (flag=0) components
  int n_pseudo;               // Number of lumped pseudocomponents
};

// ============================================================
// Lumping method selection
// ============================================================
enum class LumpMethod { Whitson, Pedersen, Danesh };

inline String lumpName(LumpMethod m)
{
  switch (m) {
    case LumpMethod::Whitson:  return "Whitson";
    case LumpMethod::Pedersen: return "Pedersen";
    case LumpMethod::Danesh:   return "Danesh";
  }
  return "Unknown";
}

// ============================================================
// Characterization status to string
// ============================================================
inline String statusString(CharacterizationStatus s)
{
  switch (s) {
    case CharacterizationStatus::Ok:                    return "Ok";
    case CharacterizationStatus::InvalidMM:             return "InvalidMM";
    case CharacterizationStatus::InvalidSG:             return "InvalidSG";
    case CharacterizationStatus::InvalidTb:             return "InvalidTb";
    case CharacterizationStatus::InvalidTc:             return "InvalidTc";
    case CharacterizationStatus::InvalidPc:             return "InvalidPc";
    case CharacterizationStatus::InvalidOmega:          return "InvalidOmega";
    case CharacterizationStatus::UnitSuspectedForMM:    return "UnitSuspectedForMM";
    case CharacterizationStatus::HeavyFractionSGRequired: return "HeavyFractionSGRequired";
    case CharacterizationStatus::NotImplemented:        return "NotImplemented";
    case CharacterizationStatus::NumericalFailure:      return "NumericalFailure";
  }
  return "Unknown";
}

inline String provenanceString(CharacterizationProvenance p)
{
  switch (p) {
    case CharacterizationProvenance::Unknown:                    return "Unknown";
    case CharacterizationProvenance::HeavyFractionSGCalibrated:  return "HeavyFractionSGCalibrated";
    case CharacterizationProvenance::UserProvidedWatsonK:        return "UserProvidedWatsonK";
    case CharacterizationProvenance::AutoWatsonKFromMM:          return "AutoWatsonKFromMM";
    case CharacterizationProvenance::InternalAlkaneMM:           return "InternalAlkaneMM";
    case CharacterizationProvenance::ProvidedTbAndSG:            return "ProvidedTbAndSG";
    case CharacterizationProvenance::ProvidedSGEstimatedTb:      return "ProvidedSGEstimatedTb";
    case CharacterizationProvenance::ProvidedTbEstimatedSG:      return "ProvidedTbEstimatedSG";
  }
  return "Unknown";
}
