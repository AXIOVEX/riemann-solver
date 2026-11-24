# The widening bounds D_k — derivation

This file is the package-relative, sufficient reference for the
tail bounds used in the widened moment balls
W_k = I_k + [0, D_k]. It reproduces the derivation of the
program's round-3 adjudication
(coscientist-adjudication-round3.md, §1.5), where it was
independently re-derived and then executed; that adjudication is
historical background for this release and is **not** required
for replay or verification — everything needed is stated here
and in the header of `e2cert3.py`, whose `widening_bounds()`
evaluates the closed form below.

## Setup

Φ(u) = Σ_{n≥1} φ_n(u) with reduced summand
φ_n(u) = q_n(2q_n e^{2u} − 3) e^{5u/2 − q_n e^{2u}}, q_n = πn².
The historical scripts enclosed only the six-term truncation's
moments I_k; the omitted contribution to a_k on the integration
interval is

δa_k = (1/(2k)!) ∫₀⁸ Σ_{n≥7} φ_n(u) u^{2k} du ≥ 0.

(The omitted sum's part over [8,∞) needs no separate bound: it
rides inside the integration-domain tail ball the moment
construction already adds, and the confirmation enclosures below
covered it separately as well, at the 1e-900 scale.) The goal is
a proved bound δa_k ≤ D_k for every k = 0,…,32.

## Derivation

1. **Term majorant.** 0 < φ_n(u) ≤ 2q_n² e^{(9/2)u − q_n e^{2u}}:
   positivity from 2q_n e^{2u} − 3 ≥ 2π − 3 > 0; the upper bound
   drops the negative −3q_n e^{5u/2} part of the prefactor.
2. **Exponent bound.** e^{2u} ≥ 1 + 2u gives
   (9/2)u − q_n e^{2u} ≤ −q_n − (2q_n − 9/2)u, so
   φ_n(u) ≤ 2q_n² e^{−q_n} e^{−(2q_n − 9/2)u}.
3. **Integration.** The integrand is positive, so integrating
   u^{2k} over [0,∞) ⊇ [0,8] only enlarges:
   δa_k ≤ Σ_{n≥7} 2q_n² e^{−q_n} / (2q_n − 9/2)^{2k+1}.
4. **Ratio bound.** The successive-term ratio of that sum is at
   most ((n+1)/n)⁴ e^{−π(2n+1)} (the denominator ratio
   (2q_n − 9/2)/(2q_{n+1} − 9/2) ≤ 1), hence for n ≥ 7 at most
   r* = (8/7)⁴ e^{−15π} = 5.8387849549…e-21.
5. **Geometric sum.** Summing from the n = 7 term:

   D_k = 2(49π)² e^{−49π} / [(98π − 9/2)^{2k+1} (1 − r*)].

   Evaluated in Arb at 256 bits: D_0 = 2.184234e-65,
   D_1 = 2.373212e-70, D_2 = 2.578540e-75. The full vector
   (D_0,…,D_32) is recorded, as ball strings, both in
   `results_adj3.json` and in the replay payload
   `results_e2cert3.json` (33/33 string-identical; see
   FULL-GRID-REPLAY.md).

## Confirmation against independent enclosures

The bound is not merely proved, it was numerically confirmed in
the adjudication: the omitted contributions δa_k were enclosed
directly in Arb at 512 bits (acb_integral of the n = 7..25 tail
callback over [0,8], plus a covering ball [0, 1e-900] for the
n > 25 and [8,∞) parts), giving δa_0 = 2.148808e-65,
δa_1 = 2.259718e-70, δa_2 = 2.316534e-75, decreasing thereafter.
The inequality δa_k ≤ D_k holds for **every** k = 0,…,32
(tightest at k = 0, where D_0 exceeds δa_0 by only 1.6%). The
δa_k enclosures are recorded in `results_adj3.json`
(`failure_set`).

## Consequence

Since trunc(a_k) ∈ I_k by the original acb_integral enclosure
and 0 ≤ δa_k ≤ D_k by the bound above, the true moment a_k lies
in W_k = I_k + [0, D_k] for every k. Containment is proved, not
asserted, and the Sturm certificates run from W_k are
certificates about the true kernel's coefficients.
