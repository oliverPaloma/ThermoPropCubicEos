#include "Lumping.hpp"




auto applyLumping(LumpMethod lumpMethod, CorrMethod corrMethod, const String &components, Vec<real> &z,
                  Vec<real> &MM, Vec<int> &flag, Vec<real> &Tc, Vec<real> &Pc, Vec<real> &omega, int npseudos, int nC_first, int nC_last)-> void{
    Vec<Component> comps; 

    auto ncomp = z.size();

    std::stringstream ss(components);//
            std::vector<std::string> compNames;//
            std::string name;//
            while (ss >> name)//
            compNames.push_back(name);//
    
    for (int i = 0; i < ncomp; i++){
            Component c{};
            //
            c.name  = compNames[i]; //trocado 11/06/26
            c.z     = val(z[i]);
            c.MM    = val(MM[i]);
            c.flag = flag[i];
            c.Tc    = val(Tc[i]);
            c.Pc    = val(Pc[i]);
            c.omega = val(omega[i]);
            comps.push_back(c);}

    Vec<Component> lights, heavies;

    for (const auto &c : comps){ //separa os componentes leves e pesados com base na flag
        if (c.flag == 0)
        lights.push_back(c);
        else
        heavies.push_back(c);}
 
      std::cout << "lights = " << lights.size() << std::endl;
      std::cout << "heavies = " << heavies.size() << std::endl;
 

  Vec<LumpedPseudo> pseudos;

  switch (lumpMethod){
    case LumpMethod::Pedersen:
        lumpPedersen(heavies, npseudos, pseudos);
        break;
    case LumpMethod::Danesh:
        lumpDanesh(heavies, npseudos, pseudos);
        break;
    case LumpMethod::Whitson:
        lumpWhitson(heavies, nC_first, nC_last, pseudos);
        break;
      }


    Mixture lump_mix;
    buildLumpedMixture(lights, pseudos, corrMethod, lump_mix);

    z.clear();
    Tc.clear();
    Pc.clear();
    omega.clear();
    MM.clear();
    flag.clear();

    for (const auto& c : lump_mix.comps){
    z.push_back(c.z);
    Tc.push_back(c.Tc);
    Pc.push_back(c.Pc);
    omega.push_back(c.omega);
    MM.push_back(c.MM);
    flag.push_back(c.flag);}
    
    z.shrink_to_fit();
    Tc.shrink_to_fit();
    Pc.shrink_to_fit();
    omega.shrink_to_fit();
    MM.shrink_to_fit();
    flag.shrink_to_fit();
}


// ============================================================
// Whitson (1993) — logarithmic cut-point spacing
// npseudos = 1 + 3.3 * log10(nC_last - nC_first)
// Cut points are logarithmically spaced between MM_min and MM_max.
// Components are assigned to the pseudo whose range contains their MM.
// Pseudo properties are molar-average of members.
// ============================================================
auto lumpWhitson(const Vec<Component> &heavy, int nC_first, int nC_last, Vec<LumpedPseudo> &pseudos) -> int{ 
  assert(!heavy.empty() && "Heavy fraction vector is empty");
  assert(nC_last > nC_first && "nC_last must be > nC_first");

  auto npseudos = static_cast<int>(1.0 + 3.3 * std::log10(nC_last - nC_first));
  // Log-spaced cut points
  auto MM_min = heavy.front().MM;
  auto MM_max = heavy.back().MM;
  auto ratio = MM_max / MM_min;

  Vec<double> MMcut(npseudos + 1);
  for (int i = 0; i <= npseudos; i++) {
    MMcut[i] = MM_min * std::pow(ratio, static_cast<double>(i) / npseudos);
  }

  // Accumulate z per pseudo (first pass — just sums)
  // Use semi-open intervals to avoid duplicate assignment at boundaries:
  //   groups j < npseudos-1: MMcut[j] <= MM < MMcut[j+1]
  //   last group:            MMcut[npseudos-1] <= MM <= MMcut[npseudos]
  Vec<double> z_sum(npseudos, 0.0);
  for (const auto &c : heavy) {
    for (int j = 0; j < npseudos; j++) {
      bool in_group = false;
      if (j < npseudos - 1) {
        in_group = (c.MM >= MMcut[j] && c.MM < MMcut[j + 1]);
      } else {
        in_group = (c.MM >= MMcut[j] && c.MM <= MMcut[j + 1]);
      }
      if (in_group) {
        z_sum[j] += c.z;
        break;
      }
    }
  }

  // Build pseudos with molar-averaged properties
  pseudos.resize(npseudos);
  for (int j = 0; j < npseudos; j++) {
    pseudos[j].z = z_sum[j];
    pseudos[j].MM = 0.0;
    pseudos[j].Tb = 0.0;
    pseudos[j].SG = 0.0;
    pseudos[j].n_members = 0;
    pseudos[j].member_indices.clear();

    if (z_sum[j] > 1e-30) {
      for (int ci = 0; ci < static_cast<int>(heavy.size()); ci++) {
        const auto &c = heavy[ci];
        bool in_group = false;
        if (j < npseudos - 1) {
          in_group = (c.MM >= MMcut[j] && c.MM < MMcut[j + 1]);
        } else {
          in_group = (c.MM >= MMcut[j] && c.MM <= MMcut[j + 1]);
        }
        if (in_group) {
          auto phi = c.z / z_sum[j];
          pseudos[j].MM += phi * c.MM;
          pseudos[j].Tb += phi * c.Tb;
          pseudos[j].SG += phi * c.SG;
          pseudos[j].n_members++;
          pseudos[j].member_indices.push_back(ci);
        }
      }
    }
  }

  return npseudos;
}

// ============================================================
// Pedersen — equal MM*z per group
// Total mixture molar mass (MM*z sum) is divided equally among
// npseudos groups. Components are traversed in order (assumed
// sorted by increasing MM) and accumulated until the group
// target is reached. At the boundary, the algorithm decides
// whether to include or exclude the boundary component based
// on which gives a closer match to the target.
// ============================================================

auto lumpPedersen(const Vec<Component> &heavy,int npseudos,Vec<LumpedPseudo> &pseudos) -> void{
  assert(!heavy.empty() && "Heavy fraction vector is empty");
  assert(npseudos > 0 && "npseudos must be > 0");

  auto n = static_cast<int>(heavy.size());

  // Total MM*z
  auto MMtotal = 0.0;
  for (const auto &c : heavy) MMtotal += c.MM * c.z;

  auto MMgroup = MMtotal / static_cast<double>(npseudos);

  pseudos.resize(npseudos);
  auto idx = 0;

  for (int j = 0; j < npseudos; j++) {
     
    auto accum = 0.0;
    auto z_accum = 0.0;
    Vec<int> members;
    int end_idx = idx;  // tracks where this group ends

    for (int i = idx; i < n; i++) {
      accum += heavy[i].MM * heavy[i].z;
      z_accum += heavy[i].z;
      members.push_back(i);

      // Check if we've exceeded the group target
      if (accum > MMgroup && j < npseudos - 1) {
        auto plus = std::abs(accum - MMgroup);
        auto minus = std::abs(accum - MMgroup - heavy[i].MM * heavy[i].z);

        if (minus < plus) {
          // Exclude component i from this group
          members.pop_back();
          z_accum -= heavy[i].z;
          accum -= heavy[i].MM * heavy[i].z;
          end_idx = i;  // component i starts next group
        }
        else {
          end_idx = i + 1;  // component i+1 starts next group
        }
        break;
      }

      // Last group absorbs everything remaining
      if (j == npseudos - 1 && i == n - 1) {
        end_idx = n;
      }
    }

    // Build pseudo from accumulated members
    pseudos[j].member_indices = members;
    pseudos[j].n_members = static_cast<int>(members.size());
    if (z_accum > 1e-30) {
      for (int k : members) {
        auto phi = heavy[k].z / z_accum;
        pseudos[j].MM += phi * heavy[k].MM;
        pseudos[j].Tb += phi * heavy[k].Tb;
        pseudos[j].SG += phi * heavy[k].SG;
      }
    }
    pseudos[j].z = z_accum;
    idx = end_idx;
    
    
//teste
    //std::cout << "\nPseudo " << j << " members:\n";
    //for (auto k : members){
    //std::cout << heavy[k].name << " MM=" << heavy[k].MM << " z=" << heavy[k].z << std::endl;}




  }
}

// ============================================================
// Danesh — equal ln(MM)*z per group
// Identical algorithm to Pedersen, but the grouping criterion
// is ln(MM)*z instead of MM*z. This reduces the dominance of
// very heavy components (C30+) in the grouping decision.
// ============================================================
auto lumpDanesh(const Vec<Component> &heavy,int npseudos,Vec<LumpedPseudo> &pseudos) -> void{
  assert(!heavy.empty() && "Heavy fraction vector is empty");
  assert(npseudos > 0 && "npseudos must be > 0");

  auto n = static_cast<int>(heavy.size());

  // Total ln(MM)*z
  auto lnTotal = 0.0;
  for (const auto &c : heavy) lnTotal += std::log(c.MM) * c.z;

  auto lnGroup = lnTotal / static_cast<double>(npseudos);

  pseudos.resize(npseudos);
  auto idx = 0;

  for (int j = 0; j < npseudos; j++) {
    auto accum = 0.0;
    auto z_accum = 0.0;
    Vec<int> members;
    int end_idx = idx;

    for (int i = idx; i < n; i++) {
      accum += std::log(heavy[i].MM) * heavy[i].z;
      z_accum += heavy[i].z;
      members.push_back(i);

      // Check if we've exceeded the group target
      if (accum > lnGroup && j < npseudos - 1) {
        auto plus = std::abs(accum - lnGroup);
        auto minus = std::abs(accum - lnGroup - std::log(heavy[i].MM) * heavy[i].z);

        if (minus < plus) {
          // Exclude component i from this group
          members.pop_back();
          z_accum -= heavy[i].z;
          accum -= std::log(heavy[i].MM) * heavy[i].z;
          end_idx = i;
        }
        else {
          end_idx = i + 1;
        }
        break;
      }

      // Last group absorbs everything remaining
      if (j == npseudos - 1 && i == n - 1) {
        end_idx = n;
      }
    }

    // Build pseudo from accumulated members
    pseudos[j].member_indices = members;
    pseudos[j].n_members = static_cast<int>(members.size());
    if (z_accum > 1e-30) {
      for (int k : members) {
        auto phi = heavy[k].z / z_accum;
        pseudos[j].MM += phi * heavy[k].MM;
        pseudos[j].Tb += phi * heavy[k].Tb;
        pseudos[j].SG += phi * heavy[k].SG;
      }
    }
    pseudos[j].z = z_accum;
    idx = end_idx;
  }
}




// bulding lumped mixture from lights + pseudos
// ============================================================
static auto buildLumpedMixture(const Vec<Component> &lights,
                               const Vec<LumpedPseudo> &pseudos,
                               CorrMethod corrMethod,
                               Mixture &mix) -> void{
  mix.n_light = static_cast<int>(lights.size());
  mix.n_pseudo = static_cast<int>(pseudos.size());
  mix.comps.reserve(lights.size() + pseudos.size());

  for (const auto &c : lights) mix.comps.push_back(c);

  for (int i = 0; i < mix.n_pseudo; i++) {
    Component pc;//variavel do tipo Component
    pc.name = "Pseudo" + std::to_string(i + 1);
    pc.z = pseudos[i].z;
    pc.MM = pseudos[i].MM;
    pc.Tb = pseudos[i].Tb;
    pc.SG = pseudos[i].SG;
    pc.flag = 1;
    // Pseudocomponents use the same correlation method as the original
    // mixture. InternalAlkaneMM is a known limitation for pseudos because
    // they are mixtures of real SCN fractions, not pure alkanes.
    // See documentation.tex "Limitações Conhecidas e Recomendações" section.
    CharacterizationOptions opts;
    opts.method = corrMethod;
    opts.watson_k_policy = WatsonKPolicy::AutoSmoothFromMM;
    auto st = characterizeComponent(pc, opts);
    //teste
    std::cout<< pc.name<< " MM=" << pc.MM<< " SG=" << pc.SG<< " Tb=" << pc.Tb<< " Tc=" << pc.Tc<< " Pc=" << pc.Pc<< " omega=" << pc.omega<< std::endl;
    
    if (st != CharacterizationStatus::Ok) {
      std::cerr << "Warning: characterization failed for '" << pc.name
                << "': " << pc.characterization_message << "\n";
    }
    // Warn when using InternalAlkaneMM for pseudocomponents
    if (corrMethod == CorrMethod::InternalAlkaneMM && st == CharacterizationStatus::Ok) {
      pc.characterization_message += " [WARNING: InternalAlkaneMM for pseudocomponents — see documentation]";
    }
    mix.comps.push_back(pc);
  }

  auto n = static_cast<int>(mix.comps.size());
  mix.kij.assign(n, Vec<double>(n, 0.0));

  


}
