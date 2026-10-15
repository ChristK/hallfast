// =============================================================================
//  reference_ode.cpp
//
//  Independent high-precision reference solver for the Hall et al. (2011)
//  adult model (web appendix eqs 1-9): the full four-compartment ODE system
//  (adaptive thermogenesis AT, extracellular fluid ECF, glycogen G and lean
//  mass L; fat mass from the Forbes relation), integrated jointly with
//  classical fourth-order Runge-Kutta at a small fixed step.
//
//  Unlike hall_model.h it uses NO closed-form solutions, NO quasi-steady-state
//  reductions and NO change of variables, so agreement between the two checks
//  the algebra of the fast solver. Dietary inputs are piecewise constant on
//  [days*y, days*(y+1)) with no look-ahead at interval boundaries.
//
//  bw_compat = false: appendix (rho_G * dG/dt; K from EE(0) = EI(0)).
//  bw_compat = true : bw v1.0.0 (dG/dt in kg/d; K from RMR x PAL).
//  k_app (kcal/d per kg, n x Y): appetite feedback of Hall (2024); realised
//  intake change x = max(-EI0, dEI - k (BW - BW0)) is evaluated from the
//  CURRENT full state at every RK4 stage.
//
//  Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos.
//  Licence: GNU General Public License version 3 (GPL-3); see LICENSE.md.
// =============================================================================
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Internal; use hall_bw_reference() from R.
// [[Rcpp::export(.hall_reference_core)]]
List hall_reference_core(NumericVector bw, NumericVector ht, NumericVector age,
                         IntegerVector female, NumericVector PAL,
                         NumericVector EI, NumericVector fat, NumericVector pcarb_base,
                         NumericMatrix dEI, NumericMatrix dNa, NumericMatrix pcarb,
                         NumericMatrix k_app, double days, int steps_per_day,
                         bool bw_compat) {
  // Constants of the appendix (kJ converted to kcal with 1/4.184)
  const double roG = 4206.501, Na = 3220, zetaNa = 3000, zetaCI = 4000,
    roF = 9440.727, roL = 1816.444, gammaF = 3.107075, gammaL = 21.98853,
    etaF = 179.2543, etaL = 229.4455, betaTEF = 0.1, betaAT = 0.14, tauAT = 14.0;
  const double C = 10.4 * (roL / roF);
  const double alfa1 = -(1 + etaL / roL) * C, alfa2 = -(1 + etaF / roF);
  const double sG = bw_compat ? 1.0 : roG;
  const int n = bw.size(), Y = dEI.ncol();
  const double h = 1.0 / steps_per_day;
  const long nstep = std::lround(days * steps_per_day);
  NumericMatrix BWo(n, Y), TISo(n, Y);
  std::fill(BWo.begin(), BWo.end(), NA_REAL);
  std::fill(TISo.begin(), TISo.end(), NA_REAL);

  for (int i = 0; i < n; i++) {
    if (NumericVector::is_na(bw[i]) || NumericVector::is_na(ht[i]) ||
        NumericVector::is_na(age[i]) || NumericVector::is_na(PAL[i]) ||
        female[i] == NA_INTEGER) continue;
    const bool fem = female[i] == 1;
    const double rmr = 9.99 * bw[i] + 625 * ht[i] - 4.92 * age[i] + (fem ? -161.0 : 5.0);
    const double ecf0 = fem ? (-4.0 + 5.98 * ht[i] + 0.167 * bw[i])
                            : (0.025 * age[i] + 9.57 * ht[i] + 0.191 * bw[i] - 12.4);
    const double G0 = 0.5;
    double F0 = fat[i];
    if (NumericVector::is_na(F0)) {
      const double lb = std::log(bw[i] / (ht[i] * ht[i]));
      F0 = fem ? bw[i] * (0.14 * age[i] + 39.96 * lb - 102.01) / 100
               : bw[i] * (0.14 * age[i] + 37.31 * lb - 103.94) / 100;
    }
    const double L0 = bw[i] - (ecf0 + F0 + 3.7 * G0);
    const double EI0 = NumericVector::is_na(EI[i]) ? rmr * PAL[i] : EI[i];
    const double delta = ((1 - betaTEF) * PAL[i] - 1) * rmr / bw[i];
    const double K = (bw_compat ? rmr * PAL[i] : EI0) - gammaL * L0 - gammaF * F0 - delta * bw[i];
    const double CIb = pcarb_base[i] * EI0, kG = CIb / (G0 * G0);

    double y[4] = {0.0, ecf0, G0, L0};   // AT, ECF, G, L
    for (int yr = 0; yr < Y; yr++) {
      const double u = dEI(i, yr), dn = dNa(i, yr), pc = pcarb(i, yr), kk = k_app(i, yr);
      if (std::isnan(u) || std::isnan(dn) || std::isnan(pc) || std::isnan(kk)) break;
      auto deriv = [&](const double* s, double* ds) {
        const double AT = s[0], ECF = s[1], G = s[2], L = s[3];
        const double F = F0 * std::exp(roL * (L - L0) / (roF * C));
        const double W = L + F + ECF + 3.7 * G;
        const double d  = std::max(-EI0, u - kk * (W - bw[i]));   // realised dEI
        const double TI = EI0 + d, CI = pc * TI;
        const double dG = (CI - kG * G * G) / roG;
        ds[0] = (betaAT * d - AT) / tauAT;
        ds[1] = (dn - zetaNa * (ECF - ecf0) - zetaCI * (1.0 - CI / CIb)) / Na;
        ds[2] = dG;
        const double R3 = K + delta * W + betaTEF * d + AT - TI + sG * dG;
        ds[3] = (C / roL) * (R3 + gammaL * L + gammaF * F) / (alfa1 + alfa2 * F);
      };
      for (long k = 0; k < nstep; k++) {
        double k1[4], k2[4], k3[4], k4[4], tmp[4];
        deriv(y, k1);
        for (int j = 0; j < 4; j++) tmp[j] = y[j] + 0.5 * h * k1[j];
        deriv(tmp, k2);
        for (int j = 0; j < 4; j++) tmp[j] = y[j] + 0.5 * h * k2[j];
        deriv(tmp, k3);
        for (int j = 0; j < 4; j++) tmp[j] = y[j] + h * k3[j];
        deriv(tmp, k4);
        for (int j = 0; j < 4; j++) y[j] += h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0;
      }
      const double F = F0 * std::exp(roL * (y[3] - L0) / (roF * C));
      BWo(i, yr)  = y[3] + F + y[1] + 3.7 * y[2];
      TISo(i, yr) = y[3] + F;
    }
  }
  return List::create(Named("BW") = BWo, Named("tissue") = TISo);
}
