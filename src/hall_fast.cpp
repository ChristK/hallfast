// =============================================================================
//  hall_fast.cpp
//
//  Rcpp bindings for the fast annual-resolution solver of the Hall et al.
//  (Lancet 2011) adult body-weight model. The model itself (equations, closed
//  forms, integrator and documentation) is in hall_model.h, which has no R
//  dependency and can be reused in a C++ microsimulation core.
//
//  PARALLELISM
//  -----------
//  Individuals are independent. If compiled with OpenMP (-fopenmp), the outer
//  loop runs in parallel; otherwise the pragmas are ignored. No R API calls are
//  made inside the parallel region.
//
//  Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos.
//  Licence: GNU General Public License version 3 (GPL-3); see LICENSE.md.
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include "hall_model.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;


// Core solver (internal). Use the R wrappers hall_bw_annual() / hall_bw_long().
//
// Individual i occupies rows  start[i] + y*stride,  y = 0..len[i]-1  of the
// flat input/output arrays (column-major n x Y matrix: start = i, stride = n;
// long format sorted by id: start = first row of the id, stride = 1).
// dNa, pcarb and k_app may be given per row (kind 2, same layout as dEI), per
// individual (kind 1, length n) or as a single value (kind 0); this avoids
// expanding constants to full matrices for very large populations.
// k_app is the appetite feedback gain (kcal/d per kg); 0 = no feedback (the
// Lancet 2011 model), > 0 = Hall (2024) proportional feedback, in which case
// dEI is the intervention push u = -P and the realised intake change is
// u - k * (BW - BW0).
// [[Rcpp::export(.hall_core)]]
List hall_core(NumericVector bw, NumericVector ht, NumericVector age,
               IntegerVector female, NumericVector PAL, NumericVector EI,
               NumericVector fat, NumericVector pcarb_base,
               NumericVector dEI, NumericVector dNa, int dNa_kind,
               NumericVector pcarb, int pcarb_kind,
               NumericVector k_app, int k_kind,
               IntegerVector start, IntegerVector len, int stride,
               double days, double atol, bool bw_compat,
               bool full_output, int nthreads) {
  const int n = bw.size();
  const R_xlen_t nrow = dEI.size();
  // Glycogen energy term: rho_G * dG/dt (appendix) or dG/dt (bw v1.0.0)
  const double s_glyc = hall::glycogen_scale(bw_compat);

  // Outputs (allocated here, filled through raw pointers inside threads)
  NumericVector BW(nrow, NA_REAL);
  NumericVector FM, LM, ECFo, GLY, ATo, EIR;
  IntegerVector NEV;
  if (full_output) {
    FM = NumericVector(nrow, NA_REAL);  LM  = NumericVector(nrow, NA_REAL);
    ECFo = NumericVector(nrow, NA_REAL); GLY = NumericVector(nrow, NA_REAL);
    ATo = NumericVector(nrow, NA_REAL); NEV = IntegerVector(nrow, NA_INTEGER);
    EIR = NumericVector(nrow, NA_REAL);
  }
  NumericVector EI0(n, NA_REAL), F0(n, NA_REAL), L0(n, NA_REAL);

  // Raw pointers (thread-safe access; no R API inside the parallel loop)
  const double *pbw = bw.begin(), *pht = ht.begin(), *page = age.begin(),
               *pPAL = PAL.begin(), *pEI = EI.begin(), *pfat = fat.begin(),
               *ppb = pcarb_base.begin(), *pdEI = dEI.begin(),
               *pdNa = dNa.begin(), *ppc = pcarb.begin(), *pk = k_app.begin();
  const int *pfem = female.begin(), *pst = start.begin(), *plen = len.begin();
  double *pBW = BW.begin(), *pEI0 = EI0.begin(), *pF0 = F0.begin(), *pL0 = L0.begin();
  double *pFM = full_output ? FM.begin() : nullptr, *pLM = full_output ? LM.begin() : nullptr,
         *pECF = full_output ? ECFo.begin() : nullptr, *pG = full_output ? GLY.begin() : nullptr,
         *pAT = full_output ? ATo.begin() : nullptr, *pEIR = full_output ? EIR.begin() : nullptr;
  int *pNEV = full_output ? NEV.begin() : nullptr;

#ifdef _OPENMP
  if (nthreads < 1) nthreads = 1;
  #pragma omp parallel for schedule(dynamic, 512) num_threads(nthreads)
#endif
  for (int i = 0; i < n; ++i) {
    if (std::isnan(pbw[i]) || std::isnan(pht[i]) || std::isnan(page[i]) ||
        std::isnan(pPAL[i]) || std::isnan(ppb[i]) || pfem[i] == NA_INTEGER) continue;

    const hall::Baseline b = hall::make_baseline(pbw[i], pht[i], page[i], pfem[i],
                                                 pPAL[i], pEI[i], pfat[i], ppb[i], bw_compat);
    pEI0[i] = b.EI0; pF0[i] = b.F0; pL0[i] = b.L0;

    hall::State s = hall::initial_state(b);         // AT(0) = 0 (energy balance)
    double h_hint = days;                            // first-step hint
    for (int y = 0; y < plen[i]; ++y) {
      const R_xlen_t r = (R_xlen_t)pst[i] + (R_xlen_t)y * stride;
      const double vNa = pdNa[dNa_kind == 2 ? r : (dNa_kind == 1 ? i : 0)];
      const double vpc = ppc[pcarb_kind == 2 ? r : (pcarb_kind == 1 ? i : 0)];
      const double vk  = pk[k_kind == 2 ? r : (k_kind == 1 ? i : 0)];
      if (std::isnan(pdEI[r]) || std::isnan(vNa) || std::isnan(vpc) || std::isnan(vk)) break;
      int ne = 0;
      const bool ok = (vk > 0.0)
        ? hall::advance_feedback(b, s, pdEI[r], vNa, vpc, vk, days, atol, s_glyc, h_hint, ne)
        : hall::advance(b, s, pdEI[r], vNa, vpc, days, atol, s_glyc, h_hint, ne);
      if (!ok) break;
      const double F = hall::fat_mass(b, s);
      pBW[r] = hall::body_weight(b, s);
      if (full_output) {
        pFM[r] = F; pLM[r] = s.L; pECF[r] = s.ECF; pG[r] = s.G; pAT[r] = s.AT; pNEV[r] = ne;
        pEIR[r] = s.dEI;
      }
    }
  }

  List out = List::create(Named("BW") = BW, Named("EI0") = EI0,
                          Named("F0") = F0, Named("L0") = L0);
  if (full_output) {
    out["FM"] = FM; out["LM"] = LM; out["ECF"] = ECFo; out["GLY"] = GLY;
    out["AT"] = ATo; out["nevals"] = NEV; out["EIR"] = EIR;
  }
  return out;
}

// Report whether this build has OpenMP (internal).
// [[Rcpp::export(.hall_has_openmp)]]
bool hall_has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
