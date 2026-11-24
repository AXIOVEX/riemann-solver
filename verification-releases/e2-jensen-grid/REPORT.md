# E2 Report — Jensen hyperbolicity grid, with E1 rung 2 (Newton) cross-check (2026-11-22)

**Correction (2026-11-23, Co-scientist round-2 adjudication R2-1, SUSTAINED; Spec 005 R5-01).** The polynomials computed in this run are not the standard Jensen polynomials of Ξ. With H(z) := ∫₀^∞ Φ(u)cosh(u√z)du = Σ a_k z^k = Σ γ_k z^k/k!, the standard Jensen family is J_H^{d,n}(x) = Σ C(d,k) γ_{n+k} x^k with γ_k = k!·a_k. The polynomials computed here, Σ C(d,k) a_{n+k} x^k, are the Jensen polynomials of B(z) = Σ a_k z^k/k!, the Borel transform of H. The original Reading below said: "The Jensen polynomials of Ξ are hyperbolic on the whole registered grid d ≤ 8, n ≤ 10" and "one non-hyperbolic J^{d,n} anywhere would refute RH itself". Corrected: the 77/77 numerical hyperbolicity and the Newton (rung-2) results in this report are results about the Borel family. The falsifier direction still holds for this family, but by a different argument: {1/k!} is a multiplier sequence (Pólya–Schur; its generating function Σ z^k/(k!)² = I₀(2√z) lies in the Laguerre–Pólya class), so the Borel transform preserves the LP class, and under RH every Borel-family Jensen polynomial is hyperbolic. The standard ξ Jensen grid was certified separately in the round-2 adjudication's consequence run (231/231; see the addendum in EXTENSION.md and results_e2gamma.json, sha256 438d4e7f…). The identification phrases in the Reading below are corrected in place under this note.

**Venue:** tower-compute (load ~9 on 32 cores at dispatch; tower-external not needed), run dir ~/riemann-runs/e2-20261122T143947, TC-03 venv, mpmath 1.4.1 at **100 dps**. **Payload:** results_e2.json (73,521 bytes) read back from the run directory; the verdict below is quoted from that file, not from the launch wrapper (which died with a broken pipe after launch, as it sometimes does — the detached run completed and wrote its payload at 14:41 ET).

**Classification (Spec 002 R1): NUMERICAL.** Roots are mpmath `polyroots` values at 100-digit floating point. This is **not** a Sturm/interval certification: Spec 002 R5's CERTIFIED label for the grid requires rigorous methods, and that upgrade remains open. Hyperbolicity on any finite grid, certified or not, proves nothing beyond the grid (§1 of the spec).

## Registration (fixed in the script header before computing)

- Grid: d = 2..8, n = 0..10, J^{d,n}(x) = Σ_{k=0}^{d} C(d,k) a_{n+k} x^k, with a_k = ∫₀^∞ Φ(u) u^{2k} du / (2k)! recomputed fresh at 100 dps (quadrature on [0,8], the E1 discipline). The grid was not extended after seeing results.
- **Count note:** the ranges define 7 × 11 = **77** grid points, and all 77 were computed. The payload header's "count: 63" is an arithmetic slip inherited from the batch's task statement; the registered ranges govern, and the payload's per-polynomial records are the 77.
- Newton form implemented: for c_k = C(d,k) a_{n+k}, the binomial-normalized coefficients are b_k = c_k / C(d,k) = a_{n+k}, and the inequalities are b_k² − b_{k−1} b_{k+1} ≥ 0 for 1 ≤ k ≤ d−1; reported per polynomial as the minimum normalized margin (b_k² − b_{k−1}b_{k+1}) / b_k² with its argmin k.
- Falsifier thresholds fixed in advance: relative max |Im(root)| > 10⁻⁴⁰, or normalized Newton margin < −10⁻⁴⁰, would be a falsifier candidate requiring interval certification.

## Controls — PASS, so the batch counts

J^{2,0}, J^{3,0}, J^{3,1}, J^{4,0} all hyperbolic with max |Im| = 0.0, reproducing TC-02 control C-D. The batch is not a NO VERDICT.

## Moment cross-check against E1

a_0..a_18 recomputed at 100 dps agree with the E1 payload (80 dps) to a max relative deviation of **2.2×10⁻⁸⁰** (at k = 11) — exactly the E1 run's precision floor. a_0 = 0.1242801945470785274781934349213494299518; a_18 = 9.446×10⁻⁵⁰.

## Grid result — 77 of 77 hyperbolic, no falsifier candidates

- **Every polynomial on the grid is hyperbolic at 100 dps: max |Im(root)| = 0.0 for all 77** (polyroots returned exactly real roots, as in TC-02's C-D). Worst root residual |J(root)| on the grid: 4.5×10⁻⁹⁵ (J^{8,0}).
- **Roots are widely separated everywhere:** the smallest pairwise root separation on the whole grid is **32.60** (J^{8,0}), then 37.21 (J^{7,0}), 43.37 (J^{6,0}). No near-degenerate root pair anywhere on the grid, so the hyperbolicity readings are not precision-marginal.
- **Payload verdict string:** "ALL HYPERBOLIC ON GRID (NUMERICAL, 100 dps)". Falsifier candidates: none.

## E1 rung 2 — Newton inequalities hold on the whole grid

- All 77 polynomials satisfy every Newton inequality; the smallest normalized margin on the grid is **0.074704 at J^{8,10}**, then 0.078702 (J^{7,10} and J^{8,9}), 0.083171 (J^{6,10}, J^{7,9}). Margins are ~36 orders of magnitude above the falsifier floor.
- **Pattern, and the cross-check closing:** the minimum margin of J^{d,n} is always attained at the last index k = d−1, i.e. at m = n+d−1, and equals 1 − m / ((m+1) r_m) where r_m is E1 rung 1's Turán ratio. Check at the grid's extreme: r_17 = 1.020694 gives 1 − 17/(18 × 1.020694) = 0.074704, the J^{8,10} value exactly. Rung 2 is therefore rung 1 re-expressed on the Jensen grid — the two computations agree, and the Newton margin dataset to carry forward is the single curve 1 − m/((m+1) r_m) against m, decreasing toward 0 as r_m → 1, tracked next by extending the moment range rather than by widening this grid.

## Reading

The Jensen polynomials of B, the Borel transform of H (not of Ξ directly — see the correction note at the head of this report), are hyperbolic on the whole registered grid d ≤ 8, n ≤ 10 — beyond the TC-02 control points and consistent with everything tabulated in the literature — and the Newton shadow of the Laguerre ladder holds with wide margins. This is the expected outcome under RH (one non-hyperbolic J^{d,n} anywhere would refute RH itself, for this family via the 1/k! multiplier argument in the correction note), reported as NUMERICAL evidence at 100 dps for this grid only. The open upgrades, in order: interval/Sturm certification of this same grid to earn the CERTIFIED label, and extension to higher d, n as venue time allows.
