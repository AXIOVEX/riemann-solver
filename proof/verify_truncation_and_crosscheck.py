# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Tristen Kyle Pierson / BitConcepts Research
"""Truncation-error certificate for the [0,1] verifier + cross-validation.

The [0,1] interval verifiers (verify_logconcavity_rigorous.py,
verify_logconcavity_arb.py) evaluate the 5-term partial sum Phi_5, not
Phi.  Normalization throughout: Phi = 4 * sum_n phi_n, as in both
verifiers (so margins and errors are in the same units).  This script certifies, uniformly on [0,1], that the truncation
error in Q = Phi''*Phi - (Phi')^2 is far below the certified margin.

Part 1 (uniform truncation bound).  With delta = Phi - Phi_5
= sum_{n>=6} phi_n,
    Q_Phi - Q_{Phi_5} = Phi_5*delta'' + Phi_5''*delta - 2*Phi_5'*delta'
                        + delta*delta'' - (delta')^2.
Uniform majorants for u in [0,1], n >= 6 (growth factors at u = 1,
decay e^{-pi n^2} at u = 0):
    |g_n|   <= (2 pi^2 n^4 + 3 pi n^2) e^{9/2}
    |g_n'|  <= (9 pi^2 n^4 + (15/2) pi n^2) e^{9/2}
    |g_n''| <= ((81/2) pi^2 n^4 + (75/4) pi n^2) e^{9/2}
    |E_n|   <= e^{-pi n^2},  |E_n'| <= 2 pi n^2 e^2 |E_n|,
    |E_n''| <= (4 pi n^2 e^2 + 4 pi^2 n^4 e^4) |E_n|,
give |phi_n^{(k)}| <= M_k(n); the Delta_k = sum_{n>=6} M_k(n) are
evaluated in mpmath.iv (finite sum n = 6..40 plus a certified geometric
tail).  The sup norms of Phi_5^{(k)} on [0,1] are enclosed by iv
evaluation of the finite sums on a grid of subintervals.  The propagated
total is compared against the certified margin of the rigorous [0,1]
run (read from results/verify_logconcavity_rigorous.json when present).

Part 2 (cross-validation).  Q_{Phi_5}(u) is computed in float mpmath and
in interval arithmetic at 10 points of [0,1]; the enclosures must
contain the float values and be negative.

Exit status: 0 iff Part 1 certifies (total < margin) AND Part 2 passes
(fail-closed).  Output: results/truncation_crosscheck.json.

History: the previous version bounded only delta (not delta', delta''),
compared it directly to the Q margin, never exited nonzero, and printed
its headline bound unconditionally.  The manuscript's former "safety
factor exceeds 2.9 x 10^30" rested on that non-uniform derivation.
See Spec 004 / round-1 adjudication, Finding C-3.2.
"""

import json
import os
import sys

import mpmath as mp
from mpmath import iv

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DPS = 60
N_MAJORANT = 40        # majorant sums n = 6..40, tail bounded beyond
N_GRID = 20            # subintervals for the Phi_5 sup-norm enclosures


# ---------------------------------------------------------------------------
# Closed-form phi_n and derivatives (works with mp or iv via duck typing)
# ---------------------------------------------------------------------------

def phi_n_all(n, u, lib):
    pi = lib.pi
    n2 = n * n
    n4 = n2 * n2
    e5h = lib.exp(lib.mpf(5) / 2 * u)
    e9h = lib.exp(lib.mpf(9) / 2 * u)
    e2u = lib.exp(2 * u)
    e4u = e2u * e2u
    g = 2 * pi ** 2 * n4 * e9h - 3 * pi * n2 * e5h
    gd1 = 9 * pi ** 2 * n4 * e9h - lib.mpf(15) / 2 * pi * n2 * e5h
    gd2 = lib.mpf(81) / 2 * pi ** 2 * n4 * e9h - lib.mpf(75) / 4 * pi * n2 * e5h
    E = lib.exp(-pi * n2 * e2u)
    Ed1 = -2 * pi * n2 * e2u * E
    Ed2 = (-4 * pi * n2 * e2u + 4 * pi ** 2 * n4 * e4u) * E
    return g * E, gd1 * E + g * Ed1, gd2 * E + 2 * gd1 * Ed1 + g * Ed2


def majorant(k, n):
    """iv upper bound M_k(n) >= |phi_n^{(k)}(u)| for all u in [0,1]."""
    pi = iv.pi
    n2 = iv.mpf(n * n)
    n4 = n2 * n2
    e92 = iv.exp(iv.mpf(9) / 2)
    D = iv.exp(-pi * n2)
    G0 = (2 * pi ** 2 * n4 + 3 * pi * n2) * e92
    if k == 0:
        return G0 * D
    G1 = (9 * pi ** 2 * n4 + iv.mpf(15) / 2 * pi * n2) * e92
    if k == 1:
        return (G1 + G0 * 2 * pi * n2 * iv.exp(2)) * D
    G2 = (iv.mpf(81) / 2 * pi ** 2 * n4 + iv.mpf(75) / 4 * pi * n2) * e92
    Efac = 4 * pi * n2 * iv.exp(2) + 4 * pi ** 2 * n4 * iv.exp(4)
    return (G2 + 2 * G1 * 2 * pi * n2 * iv.exp(2) + G0 * Efac) * D


def delta_bounds():
    """Uniform iv bounds Delta_k >= sup_{[0,1]} |delta^{(k)}|, k = 0,1,2,
    with delta = 4 * sum_{n>=6} phi_n (the verifiers' normalization)."""
    deltas = []
    for k in range(3):
        total = iv.mpf(0)
        for n in range(6, N_MAJORANT + 1):
            total += majorant(k, n)
        # Certified ratio bound for n >= N_MAJORANT: the majorants are
        # polynomials in n (degree <= 8) times e^{-pi n^2}, so
        # M(n+1)/M(n) <= ((n+1)/n)^8 * e^{-pi*(2n+1)} < 1/2.
        n = N_MAJORANT
        r = (iv.mpf(n + 1) / n) ** 8 * iv.exp(-iv.pi * (2 * n + 1))
        assert float(r.b) < 0.5, "majorant ratio bound not < 1/2"
        total += majorant(k, n + 1) / (1 - r)
        deltas.append(4 * total)
    return deltas


def phi5_sup_norms():
    """mpf upper bounds on sup_{[0,1]} |Phi_5^{(k)}| via a subinterval grid."""
    sups = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    for j in range(N_GRID):
        u_iv = iv.mpf([j / N_GRID, (j + 1) / N_GRID])
        sums = [iv.mpf(0), iv.mpf(0), iv.mpf(0)]
        for n in range(1, 6):
            p, pd1, pd2 = phi_n_all(n, u_iv, iv)
            sums[0] += p
            sums[1] += pd1
            sums[2] += pd2
        for k in range(3):
            sups[k] = max(sups[k], abs(4 * sums[k]).b)
    return sups


def main():
    mp.mp.dps = DPS
    iv.dps = DPS

    print("=" * 72)
    print("  TRUNCATION CERTIFICATE [0,1] + CROSS-VALIDATION")
    print("  Uniform bounds on delta, delta', delta'' (mpmath.iv, 60 digits)")
    print("=" * 72)

    # ---- Part 1: uniform truncation bound
    D0, D1, D2 = delta_bounds()
    S0, S1, S2 = phi5_sup_norms()
    total = D2 * S0 + S2 * D0 + 2 * S1 * D1 + D0 * D2 + D1 ** 2
    total_hi = float(total.b)

    margin = 3.36e-12  # fallback: the certified [0,1] margin of record
    rig_path = os.path.join(ROOT, "results", "verify_logconcavity_rigorous.json")
    if os.path.exists(rig_path):
        with open(rig_path) as f:
            rig = json.load(f)
        if rig.get("max_Q_upper") is not None:
            margin = abs(float(rig["max_Q_upper"]))

    safety = margin / total_hi if total_hi > 0 else float("inf")
    part1_ok = total_hi < margin

    print(f"  Delta_0 (sup |delta|)   <= {float(D0.b):.6e}")
    print(f"  Delta_1 (sup |delta'|)  <= {float(D1.b):.6e}")
    print(f"  Delta_2 (sup |delta''|) <= {float(D2.b):.6e}")
    print(f"  sup |Phi_5|   <= {float(S0):.6f}")
    print(f"  sup |Phi_5'|  <= {float(S1):.6f}")
    print(f"  sup |Phi_5''| <= {float(S2):.6f}")
    print(f"  Propagated |Q_Phi - Q_Phi5| <= {float(total.b):.6e} on [0,1]")
    print(f"  Certified [0,1] margin: {margin:.4e}")
    print(f"  Safety factor: {safety:.4e}")
    print(f"  Part 1: {'CERTIFIED' if part1_ok else '*** NOT CERTIFIED ***'}")

    # ---- Part 2: float vs interval cross-validation of Q at 10 points
    print()
    print("  Cross-validation (float mpmath vs mpmath.iv):")
    xv_ok = True
    for i in range(10):
        u_f = 0.05 + i * 0.1
        s = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        for n in range(1, 6):
            p, pd1, pd2 = phi_n_all(n, mp.mpf(u_f), mp)
            s[0] += p
            s[1] += pd1
            s[2] += pd2
        s = [4 * x for x in s]
        q_float = s[2] * s[0] - s[1] ** 2
        u_iv = iv.mpf([u_f - 1e-9, u_f + 1e-9])
        si = [iv.mpf(0), iv.mpf(0), iv.mpf(0)]
        for n in range(1, 6):
            p, pd1, pd2 = phi_n_all(n, u_iv, iv)
            si[0] += p
            si[1] += pd1
            si[2] += pd2
        si = [4 * x for x in si]
        q_iv = si[2] * si[0] - si[1] ** 2
        contains = float(q_iv.a) <= float(q_float) <= float(q_iv.b)
        negative = float(q_iv.b) < 0
        ok = contains and negative
        xv_ok = xv_ok and ok
        print(f"    u={u_f:.2f}: Q_float={float(q_float):>14.6e}  "
              f"Q_iv=[{float(q_iv.a):.4e}, {float(q_iv.b):.4e}]  "
              f"{'OK' if ok else '*** MISMATCH ***'}")
    print(f"  Part 2: {'PASSED' if xv_ok else '*** FAILED ***'}")

    passed = bool(part1_ok and xv_ok)
    print()
    print(f"  *** {'ALL CHECKS PASSED' if passed else 'FAILURES — REVIEW ABOVE'} ***")
    print("=" * 72)

    result = {
        "delta_0_bound": float(D0.b),
        "delta_1_bound": float(D1.b),
        "delta_2_bound": float(D2.b),
        "sup_phi5": [float(s) for s in (S0, S1, S2)],
        "truncation_bound": total_hi,
        "q_phi_margin": margin,
        "safety_factor": safety,
        "cross_validation_passed": bool(xv_ok),
        "passed": passed,
    }
    out_path = os.path.join(ROOT, "results", "truncation_crosscheck.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print("  -> results/truncation_crosscheck.json")
    return result


if __name__ == "__main__":
    result = main()
    sys.exit(0 if result["passed"] else 1)
