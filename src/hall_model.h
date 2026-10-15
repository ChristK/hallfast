// =============================================================================
//  hall_model.h
//
//  Fast, annual-resolution solver for the adult dynamic body-weight model of
//  Hall KD, Sacks G, Chandramohan D, et al. Quantification of the effect of
//  energy imbalance on bodyweight. Lancet 2011; 378: 826-37 (web appendix,
//  "Dynamic Mathematical Model of Body Weight Change in Adults").
//  The web appendix is the reference specification. A compatibility mode
//  reproduces the R package 'bw' (adult_weight(), v1.0.0, Camacho-Garcia-
//  Formenti & Zepeda-Tello), which implements the same model with two
//  deviations (see MODEL VARIANTS).
//
//  PURPOSE
//  -------
//  bw::adult_weight() integrates four coupled ODEs (adaptive thermogenesis AT,
//  extracellular fluid ECF, glycogen G and lean mass L) with a fixed one-day
//  Runge-Kutta 4 step, and stores every compartment for every day. For a
//  microsimulation that only needs body weight at the end of each year this is
//  wasteful. This file solves the same equations at annual resolution with
//  about 1/30 of the right-hand-side evaluations and no daily storage.
//
//  WHY IT IS FAST (and still accurate)
//  -----------------------------------
//  Within a year, the change in energy intake (dEI), the change in sodium
//  intake (dNa) and the carbohydrate share (pcarb) are constant. Then:
//
//   1. AT, ECF and G evolve independently of L and F (their equations contain
//      only themselves and the dietary inputs). Each has an EXACT closed-form
//      solution:
//         AT(t)  = AT_inf  + (AT_s  - AT_inf ) exp(-t / tau_AT)            [eq 7]
//         ECF(t) = ECF_inf + (ECF_s - ECF_inf) exp(-t * xi_Na / [Na])      [eq 2]
//         G(t)   = G_inf (G_s + G_inf th) / (G_inf + G_s th),
//                  th = tanh(kappa t), kappa = k_G G_inf / rho_G           [eq 1]
//      (Riccati equation rho_G dG/dt = CI - k_G G^2.)
//
//   2. The partitioning p = C/(C+F), C = 10.4 kg x rho_L/rho_F (eq 3), implies
//      dL/dF = 10.4 kg / F, so F = F0 * exp((L - L0) / 10.4 kg), and fat and
//      lean mass reduce to ONE scalar ODE in L. Substituting eqs 3 and 5 gives
//         dL/dt = (C/rho_L) [EI - K - gF F - gL L - delta BW - TEF - AT
//                            - rho_G dG/dt] / [(1+eL/rho_L) C + (1+eF/rho_F) F]
//      (algebraically identical to eq 9 combined with eq 3).
//
//   3. ECF and glycogen relax within 1-2 days. Their effect on L over the year
//      is added as an exact time-integral ("kick") at the start of the year.
//      The error of that is O(|dL/dt| x 1 day), which is negligible.
//
//   4. The remaining scalar ODE for L, including the exact AT transient
//      (tau_AT = 14 d), is integrated with an adaptive Dormand-Prince 5(4)
//      embedded Runge-Kutta pair with absolute error control on L. An exact
//      change of variables, L = M - b0 tau_AT A_AT exp(-t/tau_AT), moves most
//      of the AT transient out of the integrated variable. Typically 2-4
//      steps (10-26 RHS evaluations) per person-year are enough.
//
//  MODEL VARIANTS (argument bw_compat)
//  -----------------------------------
//  bw_compat = false (default, appendix):
//    * glycogen energy enters the energy balance as rho_G * dG/dt (kcal/d),
//      as in eqs 3 and 9;
//    * K is set from the initial energy-balance condition EE(0) = EI(0),
//      using the supplied baseline EI if given, otherwise RMR x PAL.
//  bw_compat = true (reproduces bw v1.0.0):
//    * bw adds dG/dt in kg/d (missing the rho_G factor), i.e. the glycogen
//      energy term is ~4,200 times too small;
//    * bw always sets K from RMR x PAL, so a user-supplied EI that differs from
//      RMR x PAL creates an energy imbalance (weight drift) at baseline.
//  Both variants share all constants and the remaining equations, which agree
//  with the appendix (constants converted from kJ to kcal with 1/4.184).
//  The appendix does not state how initial ECF is obtained; as bw does, the
//  Silva et al. sex-specific regressions are used.
//  bw's day-365 output also contains a one-day look-ahead artefact at each
//  change in EIchange (its RK4 k4 stage reads the next day's row), and with an
//  EIchange of exactly `days` columns it simulates only days - 1 days. These
//  are numerical/indexing issues in bw, not features of the model.
//
//  STANDALONE C++ USE (no R dependency; C++11)
//  --------------------------------------------
//    #include "hall_model.h"
//    hall::Baseline b = hall::make_baseline(bw, ht, age, female, PAL,
//                                           NAN /*EI*/, NAN /*fat*/, 0.5,
//                                           false /*bw_compat*/);
//    hall::State s   = hall::initial_state(b);
//    double h = 365.0; int nev = 0;
//    for (int y = 0; y < Y; ++y) {
//      hall::advance(b, s, dEI[y], 0.0 /*dNa*/, 0.5 /*pcarb*/, 365.0, 1e-4,
//                    hall::glycogen_scale(false), h, nev);
//      double BW = hall::body_weight(b, s);   // kg at the end of year y
//    }
//  Everything here is thread-safe (no shared mutable state).
//
//  Part of the R package 'hallfast'. Copyright (C) 2026 Chris Kypridemos.
//  Licence: GNU General Public License version 3 (GPL-3); see LICENSE.md.
// =============================================================================

#ifndef HALL_MODEL_H
#define HALL_MODEL_H

#include <cmath>
#include <algorithm>

namespace hall {

// ---------------------------------------------------------------------------
// Model constants (kcal units): appendix values converted from kJ with 1/4.184;
// identical to bw::Adult::getParameters()
// ---------------------------------------------------------------------------
constexpr double roG     = 4206.501;  // glycogen energy density   17.6 MJ/kg
constexpr double Na      = 3220.0;    // [Na] extracellular, mg/L  (3.22 mg/ml)
constexpr double zetaNa  = 3000.0;    // xi_Na, mg/L/d
constexpr double zetaCI  = 4000.0;    // xi_CI, mg/d
constexpr double roF     = 9440.727;  // fat energy density        39.5 MJ/kg
constexpr double roL     = 1816.444;  // lean energy density        7.6 MJ/kg
constexpr double gammaF  = 3.107075;  // RMR per kg fat            13 kJ/kg/d
constexpr double gammaL  = 21.98853;  // RMR per kg lean           92 kJ/kg/d
constexpr double etaF    = 179.2543;  // fat synthesis cost       750 kJ/kg
constexpr double etaL    = 229.4455;  // lean synthesis cost      960 kJ/kg
constexpr double betaTEF = 0.1;       // thermic effect of feeding
constexpr double betaAT  = 0.14;      // adaptive thermogenesis
constexpr double tauAT   = 14.0;      // AT time constant (days)
constexpr double G_base  = 0.5;       // baseline glycogen (kg)
constexpr double hydG    = 3.7;       // 1 kg glycogen + 2.7 kg bound water

// Derived constants (as in bw)
constexpr double Cf      = 10.4 * (roL / roF);            // C in p = C/(C+F)
constexpr double alfa1   = -(1.0 + etaL / roL) * Cf;
constexpr double alfa2   = -(1.0 + etaF / roF);
constexpr double invForb = roL / (roF * Cf);              // = 1/10.4 per kg
constexpr double CroL    = Cf / roL;
constexpr double tauECF  = Na / zetaNa;                   // ~1.07 days

// ---------------------------------------------------------------------------
// Per-individual baseline quantities (fixed for the whole simulation)
// ---------------------------------------------------------------------------
struct Baseline {
  double ht, BW0, EI0, F0, L0, ECF0, delta, K, CIb, kG;
  bool   equilibrium;   // true if baseline is an energy-balance steady state (EE0 == EI0)
};

// Builds the baseline state and constants (appendix; bw::Adult::build()).
//  female    : 1 = female, 0 = male
//  EI_in     : baseline energy intake (kcal/d); NaN -> Mifflin-St Jeor RMR x PAL
//  fat_in    : baseline fat mass (kg); NaN -> Jackson et al. equation (eq 4)
//  bw_compat : if true, K is set from RMR x PAL even when EI_in is supplied
//              (bw behaviour); otherwise K satisfies EE(0) = EI(0) (appendix).
inline Baseline make_baseline(double bw, double ht, double age, int female,
                              double PAL, double EI_in, double fat_in,
                              double pcarb_base, bool bw_compat) {
  Baseline b;
  b.ht = ht;
  b.BW0 = bw;
  // Mifflin-St Jeor resting metabolic rate (kcal/d)
  const double rmr = 9.99 * bw + 625.0 * ht - 4.92 * age + (female ? -161.0 : 5.0);
  // Silva extracellular fluid equations (kg)
  b.ECF0 = female ? (-4.0 + 5.98 * ht + 0.167 * bw)
                  : (0.025 * age + 9.57 * ht + 0.191 * bw - 12.4);
  // Jackson et al. body fat (eq 4), unless supplied
  if (std::isnan(fat_in)) {
    const double lbmi = std::log(bw / (ht * ht));
    b.F0 = female ? bw * (0.14 * age + 39.96 * lbmi - 102.01) / 100.0
                  : bw * (0.14 * age + 37.31 * lbmi - 103.94) / 100.0;
  } else {
    b.F0 = fat_in;
  }
  b.L0    = bw - (b.ECF0 + b.F0 + hydG * G_base);
  b.EI0   = std::isnan(EI_in) ? rmr * PAL : EI_in;
  b.delta = ((1.0 - betaTEF) * PAL - 1.0) * rmr / bw;        // eq 8
  // K from the initial energy balance: EE(0) = K + gF F0 + gL L0 + delta BW = EI(0)
  const double EE0 = bw_compat ? rmr * PAL : b.EI0;
  b.K     = EE0 - gammaL * b.L0 - gammaF * b.F0 - b.delta * bw;
  // Always true for the appendix model; false for bw_compat with EI != RMR x PAL
  b.equilibrium = (EE0 == b.EI0);
  b.CIb   = pcarb_base * b.EI0;
  b.kG    = b.CIb / (G_base * G_base);
  return b;
}

// Time-varying state carried from one interval to the next.
// dEI is the realised change in energy intake from baseline (kcal/d) at the
// end of the interval; without appetite feedback it equals the input dEI.
struct State {
  double L, AT, ECF, G, dEI;
};

// log( cosh(x) + c*sinh(x) ) for x >= 0, c >= 0, without overflow.
inline double log_cosh_c_sinh(double x, double c) {
  return x + std::log(0.5 * (1.0 + c) + 0.5 * (1.0 - c) * std::exp(-2.0 * x));
}

// ---------------------------------------------------------------------------
// Advance one individual by one interval of length T days with constant
// dEI (kcal/d), dNa (mg/d) and pcarb. Returns false if integration failed.
// nevals (output) counts right-hand-side evaluations (diagnostic).
// ---------------------------------------------------------------------------
inline bool advance(const Baseline& b, State& s, double dEI, double dNa,
                    double pcarb, double T, double atol, double s_glyc,
                    double& h_hint, int& nevals) {
  const double TI = b.EI0 + dEI;                 // total intake
  const double CI = std::max(0.0, pcarb * TI);   // carbohydrate intake

  // Exact shortcut: baseline is an energy-balance equilibrium, there is no
  // dietary change and the state is still at baseline (e.g. years before an
  // intervention starts) -> nothing changes.
  if (b.equilibrium && dEI == 0.0 && dNa == 0.0 && CI == b.CIb && s.L == b.L0 &&
      s.AT == 0.0 && s.ECF == b.ECF0 && s.G == G_base) { s.dEI = 0.0; return true; }

  // Steady states of the fast compartments for this year's diet
  const double AT_inf  = betaAT * dEI;                                 // eq 7
  const double ECF_inf = b.ECF0 + (dNa - zetaCI * (1.0 - CI / b.CIb)) / zetaNa; // eq 2
  const double bG      = b.kG / roG;
  const double G_inf   = std::sqrt(CI / b.kG);                         // eq 1
  const double kappa   = bG * G_inf;

  // Transient amplitudes at the start of the year
  const double A_AT  = s.AT - AT_inf;
  const double A_ECF = s.ECF - ECF_inf;
  const double Gs    = s.G;

  // Glycogen at the end of the year (closed form of the Riccati equation)
  double G_T, intG;   // intG = integral_0^T (G(t) - G_inf) dt
  if (G_inf > 0.0) {
    const double x  = kappa * T;
    const double th = std::tanh(x);
    G_T  = G_inf * (Gs + G_inf * th) / (G_inf + Gs * th);
    intG = log_cosh_c_sinh(x, Gs / G_inf) / bG - G_inf * T;
  } else {             // zero carbohydrate intake: dG/dt = -bG G^2
    G_T  = Gs / (1.0 + bG * Gs * T);
    intG = std::log1p(bG * Gs * T) / bG;
  }

  // Part of the energy-balance numerator that does not depend on L or t:
  //   K + delta*(ECF_inf + 3.7 G_inf) + TEF + AT_inf - TI
  const double N0 = b.K + b.delta * (ECF_inf + hydG * G_inf) + betaTEF * dEI
                  + AT_inf - TI;
  const double F0 = b.F0, L0 = b.L0, delta = b.delta;

  // ---- Exact contribution of the ECF and glycogen transients ("kick") -------
  // dL/dt contains b(L) * [delta*(ECF-ECF_inf) + delta*3.7*(G-G_inf) + s*dG/dt]
  // with b(L) = (C/rho_L)/(alfa1 + alfa2 F). Both compartments relax within
  // ~1-2 days, so b(L) is constant over that time to high accuracy.
  const double Fs = F0 * std::exp(invForb * (s.L - L0));
  const double bs = CroL / (alfa1 + alfa2 * Fs);       // b(L) at the start
  s.L += bs * (delta * A_ECF * tauECF * (1.0 - std::exp(-T / tauECF))
               + delta * hydG * intG + s_glyc * (G_T - Gs));

  // ---- Change of variables removing the AT transient (exact) ---------------
  // dL/dt = a(L) + b(L) * A_AT e^{-t/tau}. With L = M - bs*tau*A_AT e^{-t/tau}
  // (bs = b at the start), dM/dt = a(L) + (b(L) - bs) A_AT e^{-t/tau}: the
  // integrator sees only the small residual (b(L) - bs) of the 14-day
  // transient. The RHS below receives atf = A_AT e^{-t/tau}.
  auto rhs = [&](double atf, double M) -> double {
    ++nevals;
    const double L = M - bs * tauAT * atf;
    const double F = F0 * std::exp(invForb * (L - L0));
    const double N = N0 + gammaF * F + gammaL * L + delta * (L + F) + atf;
    return CroL * N / (alfa1 + alfa2 * F) - bs * atf;
  };

  // ---- Adaptive Dormand-Prince 5(4) on the scalar ODE for M ----------------
  double t = 0.0, M = s.L + bs * tauAT * A_AT;
  double Et = A_AT;                               // A_AT e^{-t/tau} at the current t
  // Start with a short step only if the AT transient is large enough to matter
  // (its total effect on L is ~ |b A_AT| tau).
  const bool at_matters = std::fabs(A_AT * bs) * tauAT > 0.1 * atol;
  double h = std::min(T, at_matters ? std::min(h_hint, 7.0) : h_hint);
  double k1 = rhs(Et, M);
  int guard = 0;
  while (t < T) {
    if (++guard > 100000 || h < 1e-10) return false;   // pathological input
    if (t + h > T) h = T - t;
    // e^{-c h/tau} for the stage nodes c = 1/5, 3/10, 4/5, 8/9, 1 from one exp:
    // with q = e^{-h/(90 tau)} they are q^18, q^27, q^72, q^80, q^90.
    const double q = std::exp(-h / (90.0 * tauAT));
    const double q2 = q * q, q4 = q2 * q2, q8 = q4 * q4, q9 = q8 * q, q18 = q9 * q9,
                 q27 = q18 * q9, q36 = q18 * q18, q72 = q36 * q36, q80 = q72 * q8,
                 q90 = q80 * q8 * q2;
    const double k2 = rhs(Et * q18, M + h * (k1 / 5.0));
    const double k3 = rhs(Et * q27, M + h * (3.0*k1/40.0 + 9.0*k2/40.0));
    const double k4 = rhs(Et * q72, M + h * (44.0*k1/45.0 - 56.0*k2/15.0 + 32.0*k3/9.0));
    const double k5 = rhs(Et * q80, M + h * (19372.0*k1/6561.0 - 25360.0*k2/2187.0
                                             + 64448.0*k3/6561.0 - 212.0*k4/729.0));
    const double k6 = rhs(Et * q90, M + h * (9017.0*k1/3168.0 - 355.0*k2/33.0
                                             + 46732.0*k3/5247.0 + 49.0*k4/176.0
                                             - 5103.0*k5/18656.0));
    const double Mn = M + h * (35.0*k1/384.0 + 500.0*k3/1113.0 + 125.0*k4/192.0
                               - 2187.0*k5/6784.0 + 11.0*k6/84.0);
    const double k7 = rhs(Et * q90, Mn);                  // FSAL
    const double err = std::fabs(h * (71.0*k1/57600.0 - 71.0*k3/16695.0 + 71.0*k4/1920.0
                                      - 17253.0*k5/339200.0 + 22.0*k6/525.0 - k7/40.0));
    if (!std::isfinite(err)) return false;
    if (err <= atol) {            // accept
      t += h; M = Mn; k1 = k7; Et *= q90;
    }
    const double fac = (err > 0.0) ? 0.9 * std::pow(atol / err, 0.2) : 5.0;
    h *= std::min(5.0, std::max(0.2, fac));
  }
  h_hint = h;                     // carry a step-size hint into next year
  s.L = M - bs * tauAT * Et;      // back-transform M -> L at t = T

  // ---- Fast compartments at the end of the year (exact) --------------------
  s.AT  = AT_inf  + A_AT  * std::exp(-T / tauAT);
  s.ECF = ECF_inf + A_ECF * std::exp(-T / tauECF);
  s.G   = G_T;
  s.dEI = dEI;
  return std::isfinite(s.L);
}

// ---------------------------------------------------------------------------
// APPETITE FEEDBACK (optional extension; Hall 2024, Obesity 32:1163-8)
// ---------------------------------------------------------------------------
// Hall's later work closes the loop on intake with proportional feedback
// control of appetite (Polidori et al. 2016, Obesity 24:2289-95; Hall,
// Sanghvi & Gobel 2017, Obesity 25:2088-91; Hall 2024):
//     dI(t) = -P(t) - k * dW(t)
// where P(t) is the intervention "push" (kcal/d; positive = intake cut), dW(t)
// the change in body weight from the weight-maintenance baseline, and k the
// feedback gain (kcal/d per kg). Here the input u = -P (the change in intake
// that would occur without any appetite response; same sign convention as
// bw's EIchange), so the realised change in intake is
//     x(t) = u - k * (BW(t) - BW0),       truncated so that EI0 + x >= 0.
// TEF, AT and carbohydrate intake (hence glycogen and ECF) respond to the
// realised x(t), exactly as in the full model with time-varying intake.
//
// Because x(t) now varies within the interval, AT no longer has a closed form
// and is integrated together with L (2-D adaptive Dormand-Prince 5(4)). ECF
// and G (time constants ~1 day) are written as their quasi-steady state at
// the CURRENT intake plus the exact decaying transient triggered by the step
// change at the start of the interval:
//     ECF(t) = ECF*(x(t)) + e0 exp(-t/tau_ECF)
//     G(t)   = G*(x(t))   + [G_R(t) - G*(x0)]
// where G_R is the closed-form Riccati solution at the post-step intake x0.
// With k = 0 this is exact; with k > 0 it neglects the ~1-day lag of ECF and
// G behind their moving steady states (order tau_fast * d(steady state)/dt).
//
// Realised intake in closed form. BW contains ECF*(x) = ECFa + sE x and
// G*(x) = sqrt(pcarb (EI0 + x) / kG). With y = G*(x), x = kG y^2/pcarb - EI0,
// and x - u + k (BW - BW0) = 0 is the quadratic qa y^2 + qb y + c = 0
// (qa, qb > 0); its non-negative root is y = -2c / (qb + sqrt(qb^2 - 4 qa c)).
// c is the residual at zero intake; c >= 0 means intake is truncated at zero.
//
// Glycogen energy. rho_G dG/dt contains (dG*/dx) dx/dt, and dx/dt = -k dBW/dt
// depends on dL/dt; the resulting linear relation is solved exactly.
//
// Two treatments of the 1-2-day ECF/G step transients:
//   * "kick" (default): their effect on L over the interval, including the
//     brief intake response k * (transient water) and its TEF/AT effects, is
//     applied as an exact time-integral at the start of the interval (as in
//     advance()). The integrator then only sees the slow dynamics: ~2x faster.
//   * "resolved": the transients are kept explicitly in the right-hand side
//     and resolved by the integrator. Used automatically when the realised
//     intake after the step is below 25% of baseline (intake near zero, e.g.
//     bariatric surgery), where the kick linearisation is less accurate.
// ---------------------------------------------------------------------------
namespace detail {

// Constants of the intake solve for one interval
struct IntakeModel {
  double EI0, u, k, pcarb, kG, sE, ECFa, qa, qb, BW0;
  // Solve for the realised intake given base = BW - BW0 excluding sE*x and 3.7*G*(x).
  inline void solve(double base, double& x, double& Gst, double& dGst, bool& clamped) const {
    const double c = -EI0 - u + k * (base - sE * EI0);
    clamped = (c >= 0.0);
    if (clamped)          { x = -EI0; Gst = 0.0; }
    else if (pcarb > 0.0) { const double y = -2.0 * c / (qb + std::sqrt(qb * qb - 4.0 * qa * c));
                            Gst = y; x = kG * y * y / pcarb - EI0; }
    else                  { x = (u - k * base) / (1.0 + k * sE); Gst = 0.0; }
    dGst = (Gst > 0.0) ? pcarb / (2.0 * kG * Gst) : 0.0;
  }
};

// Adaptive Dormand-Prince 5(4) for y = (L, AT) on [0, T]. rhs(t, y, dy).
// Error norm: max(|err_L| / atol, |err_AT| / (100 atol)); an AT error of
// 0.01 kcal/d corresponds to well under 1e-4 kg of lean mass.
template <class RHS>
inline bool dp5_2d(RHS& rhs, double* y, double T, double atol, double& h) {
  const double atolAT = 100.0 * atol;
  double t = 0.0, k1[2], k2[2], k3[2], k4[2], k5[2], k6[2], k7[2], yt[2], yn[2];
  rhs(t, y, k1);
  int guard = 0;
  while (t < T) {
    if (++guard > 100000 || h < 1e-10) return false;
    if (t + h > T) h = T - t;
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + h * (k1[j] / 5.0);
    rhs(t + h / 5.0, yt, k2);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + h * (3.0*k1[j]/40.0 + 9.0*k2[j]/40.0);
    rhs(t + 3.0 * h / 10.0, yt, k3);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + h * (44.0*k1[j]/45.0 - 56.0*k2[j]/15.0 + 32.0*k3[j]/9.0);
    rhs(t + 4.0 * h / 5.0, yt, k4);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + h * (19372.0*k1[j]/6561.0 - 25360.0*k2[j]/2187.0
                                                  + 64448.0*k3[j]/6561.0 - 212.0*k4[j]/729.0);
    rhs(t + 8.0 * h / 9.0, yt, k5);
    for (int j = 0; j < 2; ++j) yt[j] = y[j] + h * (9017.0*k1[j]/3168.0 - 355.0*k2[j]/33.0
                                                  + 46732.0*k3[j]/5247.0 + 49.0*k4[j]/176.0
                                                  - 5103.0*k5[j]/18656.0);
    rhs(t + h, yt, k6);
    for (int j = 0; j < 2; ++j) yn[j] = y[j] + h * (35.0*k1[j]/384.0 + 500.0*k3[j]/1113.0
                                                  + 125.0*k4[j]/192.0 - 2187.0*k5[j]/6784.0
                                                  + 11.0*k6[j]/84.0);
    rhs(t + h, yn, k7);                                   // FSAL
    double err = 0.0;
    for (int j = 0; j < 2; ++j) {
      const double ej = std::fabs(h * (71.0*k1[j]/57600.0 - 71.0*k3[j]/16695.0 + 71.0*k4[j]/1920.0
                                       - 17253.0*k5[j]/339200.0 + 22.0*k6[j]/525.0 - k7[j]/40.0));
      err = std::max(err, ej / (j == 0 ? atol : atolAT));
    }
    if (!std::isfinite(err)) return false;
    if (err <= 1.0) { t += h; y[0] = yn[0]; y[1] = yn[1]; k1[0] = k7[0]; k1[1] = k7[1]; }
    const double fac = (err > 0.0) ? 0.9 * std::pow(1.0 / err, 0.2) : 5.0;
    h *= std::min(5.0, std::max(0.2, fac));
  }
  return true;
}

} // namespace detail

inline bool advance_feedback(const Baseline& b, State& s, double u, double dNa,
                             double pcarb, double k, double T, double atol,
                             double s_glyc, double& h_hint, int& nevals) {
  const double EI0 = b.EI0, kG = b.kG, bG = kG / roG;
  const double F0 = b.F0, L0 = b.L0, delta = b.delta, K = b.K;
  detail::IntakeModel im;
  im.EI0 = EI0; im.u = u; im.k = k; im.pcarb = pcarb; im.kG = kG; im.BW0 = b.BW0;
  im.sE   = zetaCI * pcarb / (b.CIb * zetaNa);
  im.ECFa = b.ECF0 + (dNa - zetaCI * (1.0 - pcarb * EI0 / b.CIb)) / zetaNa;
  im.qa   = (pcarb > 0.0) ? (kG / pcarb) * (1.0 + k * im.sE) : 0.0;
  im.qb   = k * hydG;
  auto Gstar = [&](double x) { const double c = pcarb * (EI0 + x); return c > 0.0 ? std::sqrt(c / kG) : 0.0; };

  // Realised intake just after the step (the state is continuous: BW(0+) = BW_s)
  const double Fs  = F0 * std::exp(invForb * (s.L - L0));
  const double BWs = s.L + Fs + s.ECF + hydG * s.G;
  const double x0  = std::max(-EI0, u - k * (BWs - b.BW0));
  // Step transients at intake x0 (exact): e(t) = e0 exp(-t/tau_ECF); g(t) = G_R(t) - G*(x0)
  const double e0    = s.ECF - (im.ECFa + im.sE * x0);
  const double CI0   = std::max(0.0, pcarb * (EI0 + x0));
  const double Ginf0 = std::sqrt(CI0 / kG), kap0 = bG * Ginf0, Gs = s.G, Gst0 = Gstar(x0);
  auto G_R = [&](double t) {
    if (Ginf0 > 0.0) { const double th = std::tanh(kap0 * t);
                       return Ginf0 * (Gs + Ginf0 * th) / (Ginf0 + Gs * th); }
    return Gs / (1.0 + bG * Gs * t);
  };
  const double bL0 = CroL / (alfa1 + alfa2 * Fs);
  const bool resolved = (EI0 + x0) < 0.25 * EI0;        // intake near zero
  double y[2] = {s.L, s.AT};

  if (!resolved) {
    // ---- "kick": integral of the step transients applied at t = 0 ----------
    double G_T, intG;     // Riccati solution at T and integral of (G_R - Ginf0)
    if (Ginf0 > 0.0) {
      const double xx = kap0 * T, th = std::tanh(xx);
      G_T  = Ginf0 * (Gs + Ginf0 * th) / (Ginf0 + Gs * th);
      intG = log_cosh_c_sinh(xx, Gs / Ginf0) / bG - Ginf0 * T;
    } else {
      G_T  = Gs / (1.0 + bG * Gs * T);
      intG = std::log1p(bG * Gs * T) / bG;
    }
    // Integral of the transient water w(t) = e(t) + 3.7 (G_R(t) - Ginf0)
    const double Iw = e0 * tauECF * (1.0 - std::exp(-T / tauECF)) + hydG * intG;
    // Intake responds to transient water: dx_tr = -k w / (1 + k S0), where S0
    // accounts for the quasi-steady ECF*/G* response to dx_tr itself.
    const double dG0   = (Gst0 > 0.0) ? pcarb / (2.0 * kG * Gst0) : 0.0;
    const double scale = 1.0 / (1.0 + k * (im.sE + hydG * dG0));
    // L: delta * BW_tr - (1 - beta_TEF) * dx_tr + glycogen energy of the transient
    y[0] += bL0 * ((delta + (1.0 - betaTEF) * k) * Iw * scale + s_glyc * (G_T - Gs));
    // AT: the ~1-day intake pulse enters AT as an impulse of area beta_AT * dx_tr
    y[1] += -betaAT * k * Iw * scale / tauAT;

    auto rhs = [&](double, const double* yy, double* dy) {
      ++nevals;
      const double L = yy[0];
      const double F = F0 * std::exp(invForb * (L - L0));
      const double base = L + F + im.ECFa - b.BW0;
      double x, Gst, dGst; bool cl;
      im.solve(base, x, Gst, dGst, cl);
      const double BW = base + b.BW0 + im.sE * x + hydG * Gst;
      const double bL = CroL / (alfa1 + alfa2 * F);
      const double Nr = K + gammaF * F + gammaL * L + delta * BW + betaTEF * x + yy[1] - (EI0 + x);
      double xdot = 0.0;
      if (!cl) {
        const double Fp = invForb * F, S = im.sE + hydG * dGst;
        xdot = -k * (1.0 + Fp) * bL * Nr / (1.0 + k * S + k * (1.0 + Fp) * bL * s_glyc * dGst);
      }
      dy[0] = bL * (Nr + s_glyc * dGst * xdot);
      dy[1] = (betaAT * x - yy[1]) / tauAT;
    };
    // Resolve the 14-day AT transient of the step only if it matters
    const double AT_jump = std::fabs(betaAT * x0 - y[1]);
    double h = std::min(T, AT_jump * std::fabs(bL0) * tauAT > 0.1 * atol ? std::min(h_hint, 7.0) : h_hint);
    if (!detail::dp5_2d(rhs, y, T, atol, h)) return false;
    h_hint = h;
    // End-of-interval state
    const double F = F0 * std::exp(invForb * (y[0] - L0));
    double x, Gst, dGst; bool cl;
    im.solve(y[0] + F + im.ECFa - b.BW0, x, Gst, dGst, cl);
    s.L = y[0]; s.AT = y[1]; s.dEI = x;
    s.ECF = im.ECFa + im.sE * x + e0 * std::exp(-T / tauECF);
    s.G   = Gst + (G_T - Gst0);
  } else {
    // ---- "resolved": transients explicit in the right-hand side ------------
    auto rhs = [&](double t, const double* yy, double* dy) {
      ++nevals;
      const double L = yy[0];
      const double F = F0 * std::exp(invForb * (L - L0));
      const double e = e0 * std::exp(-t / tauECF);
      const double GR = G_R(t), g = GR - Gst0, gp = (CI0 - kG * GR * GR) / roG;
      const double base = L + F + im.ECFa + e + hydG * g - b.BW0;
      double x, Gst, dGst; bool cl;
      im.solve(base, x, Gst, dGst, cl);
      const double BW = base + b.BW0 + im.sE * x + hydG * Gst;
      const double bL = CroL / (alfa1 + alfa2 * F);
      const double Nr = K + gammaF * F + gammaL * L + delta * BW + betaTEF * x + yy[1]
                      - (EI0 + x) + s_glyc * gp;
      double xdot = 0.0;
      if (!cl) {
        // x' = -k BW';  BW' = (1 + F/10.4) L' + (sE + 3.7 G*') x' + e' + 3.7 g'
        const double Fp = invForb * F, S = im.sE + hydG * dGst;
        xdot = -k * ((1.0 + Fp) * bL * Nr - e / tauECF + hydG * gp)
               / (1.0 + k * S + k * (1.0 + Fp) * bL * s_glyc * dGst);
      }
      dy[0] = bL * (Nr + s_glyc * dGst * xdot);
      dy[1] = (betaAT * x - yy[1]) / tauAT;
    };
    const bool fast_tr = std::fabs(e0) > 1e-9 || std::fabs(Gs - Gst0) > 1e-9;
    double h = std::min(T, fast_tr ? std::min(h_hint, 0.5) : std::min(h_hint, 7.0));
    if (!detail::dp5_2d(rhs, y, T, atol, h)) return false;
    h_hint = h;
    const double F = F0 * std::exp(invForb * (y[0] - L0));
    const double e = e0 * std::exp(-T / tauECF), g = G_R(T) - Gst0;
    double x, Gst, dGst; bool cl;
    im.solve(y[0] + F + im.ECFa + e + hydG * g - b.BW0, x, Gst, dGst, cl);
    s.L = y[0]; s.AT = y[1]; s.dEI = x;
    s.ECF = im.ECFa + im.sE * x + e;
    s.G   = Gst + g;
  }
  return std::isfinite(s.L) && std::isfinite(s.AT);
}

// ---------------------------------------------------------------------------
// Convenience helpers
// ---------------------------------------------------------------------------
// State at baseline: energy balance, so AT(0) = 0; ECF and G at baseline.
inline State initial_state(const Baseline& b) {
  return State{b.L0, 0.0, b.ECF0, G_base, 0.0};
}
// Fat mass implied by lean mass (Forbes partitioning, see point 2 above).
inline double fat_mass(const Baseline& b, const State& s) {
  return b.F0 * std::exp(invForb * (s.L - b.L0));
}
// Body weight = fat + lean + ECF + glycogen with its bound water.
inline double body_weight(const Baseline& b, const State& s) {
  return s.L + fat_mass(b, s) + s.ECF + hydG * s.G;
}
// Multiplier of dG/dt in the energy balance: rho_G (appendix) or 1 (bw v1.0.0).
inline double glycogen_scale(bool bw_compat) { return bw_compat ? 1.0 : roG; }

} // namespace hall

#endif // HALL_MODEL_H
