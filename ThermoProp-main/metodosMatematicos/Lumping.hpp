#pragma once
#include "cubic_eos_types.hpp"
#include "../includes.hpp"
#include "correlations.hpp"
// ============================================================
// Lumping methods: Whitson, Pedersen, Danesh
// These operate on the heavy fractions (flag=1) only.
// Light components (flag=0) are passed through unchanged.
// ============================================================

// ============================================================
// Lumped pseudo-component data
// ============================================================
struct LumpedPseudo
{
  double z;       // Molar fraction of the pseudo
  double MM;      // Molar mass (g/mol)
  double Tb;      // Normal boiling point (K) — molar average
  double SG;      // Specific gravity — molar average
  int    n_members;  // How many original SCN fractions were lumped
  Vec<int> member_indices;  // Indices into the heavy[] vector
};

// ============================================================
// Whitson (1993) — logarithmic cut-point spacing
// Number of pseudos is computed automatically from carbon range.
// ============================================================
auto lumpWhitson(
  const Vec<Component> &heavy,   // Only flag=1 components
  int nC_first,                   // First carbon number (e.g. 7 for C7)
  int nC_last,                    // Last carbon number (e.g. 20 for C20+)
  Vec<LumpedPseudo> &pseudos     // Output
) -> int;                         // Returns npseudos

// ============================================================
// Pedersen — equal MM*z per group
// User specifies the number of pseudocomponents.
// ============================================================
auto lumpPedersen(
  const Vec<Component> &heavy,   // Only flag=1 components
  int npseudos,                   // User-specified number
  Vec<LumpedPseudo> &pseudos     // Output
) -> void;

// ============================================================
// Danesh — equal ln(MM)*z per group
// User specifies the number of pseudocomponents.
// ============================================================
auto lumpDanesh(
  const Vec<Component> &heavy,   // Only flag=1 components
  int npseudos,                   // User-specified number
  Vec<LumpedPseudo> &pseudos     // Output
) -> void;

auto applyLumping(LumpMethod lumpMethod, CorrMethod corrMethod, const String &components, Vec<real> &z,
                  Vec<real> &MM, Vec<int> &flag, Vec<real> &Tc, Vec<real> &Pc, Vec<real> &omega, int npseudos, int nC_first, int nC_last)-> void;
static auto buildLumpedMixture(const Vec<Component> &lights,
                               const Vec<LumpedPseudo> &pseudos,
                               CorrMethod corrMethod,
                               Mixture &mix) -> void;
// ============================================================
