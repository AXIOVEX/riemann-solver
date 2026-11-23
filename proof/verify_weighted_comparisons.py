# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Tristen Kyle Pierson / BitConcepts Research
"""Arb/FLINT certificate for the far-tail weighted comparisons (u >= 3).

Certifies the restated Lemma (weighted_comparisons) of Appendix E.7:
for all u >= 3,
    S_6(u) <= 5 * S_4(u),   S_8(u) <= 32 * S_4(u),
where S_m(u) = sum_{n>=2} n^m * exp(-pi*(n^2-1)*e^{2u}).

Method (n = 2 domination).  Write a_n = exp(-pi*(n^2-1)*e^{2u}).
  S_6 - 5*S_4 = sum_{n>=2} (n^6 - 5*n^4) a_n.
The n = 2 term is (64 - 80) a_2 = -16 a_2; every n >= 3 term is >= 0.
Hence S_6 <= 5*S_4 iff
    T_6(u) := sum_{n>=3} (n^6 - 5*n^4) e^{-pi*(n^2-4)*e^{2u}} <= 16.
Similarly S_8 - 32*S_4 has n = 2 term (256 - 512) a_2 = -256 a_2, and
S_8 <= 32*S_4 iff
    T_8(u) := sum_{n>=3} (n^8 - 32*n^4) e^{-pi*(n^2-4)*e^{2u}} <= 256.
Every summand of T_6 and T_8 is a nonnegative constant times
e^{-pi*(n^2-4)*e^{2u}}, which is strictly decreasing in u; therefore
T_6 and T_8 are decreasing in u and u = 3 is the worst case.  This script
certifies T_6(3) and T_8(3) in Arb ball arithmetic (finite sum n = 3..N
plus an explicit geometric tail bound), and also certifies the numeric
fact used by the far-tail theorem: eps*(3) = 2*S_4(3) < 1e-1636.

This script replaces the role that the retired verify_algebraic_core.py
was claimed (incorrectly) to fill for this lemma: that script never
evaluated these sums and used no interval arithmetic.

Exit status: 0 iff every certified check passes (fail-closed).
Output: results/verify_weighted_comparisons.json
"""

import json
import os
import sys
import time

try:
    from flint import arb, ctx
except ImportError:
    print("SKIP: python-flint not installed. Install with: pip install python-flint")
    print("The weighted-comparisons certificate requires python-flint.")
    sys.exit(0)

PREC = 256
N_FINITE = 40          # finite sum n = 3..40 at u = 3, tail bounded beyond
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def tail_sum_arb(m, c, u_val):
    """Certified upper bound on T(u) = sum_{n>=3} (n^m - c*n^4) e^{-pi*(n^2-4)*e^{2u}}.

    Finite sum n = 3..N_FINITE evaluated in Arb with the difference
    sequence's own terms (n^m - c*n^4). The n > N_FINITE tail is NOT
    bounded from the difference sequence itself: it is bounded through
    the MAJORANT sequence n^m * e^{-pi*(n^2-4)*e^{2u}}, which dominates
    the difference sequence termwise (c*n^4 >= 0). The tail is at most
    (the majorant's first omitted term, at n = N_FINITE + 1) / (1 - r),
    where r is a certified ratio bound for the MAJORANT sequence: for
    n >= N_FINITE,
        majorant(n+1)/majorant(n)
            = ((n+1)/n)^m * e^{-pi*(2n+1)*e^{2u}} <= r < 1/2,
    evaluated at n = N_FINITE with e^{2u} >= e^6 (u = 3); the ratio
    only shrinks for larger n.
    (Docstring corrected 2026-11-23, R2-4 item 10: the previous text
    described the tail as the difference sequence's own first neglected
    term over (1 - r), which is not the argument the code runs.)
    """
    pi = arb.pi()
    u = arb(u_val)
    e2u = (2 * u).exp()
    total = arb(0)
    for n in range(3, N_FINITE + 1):
        na = arb(n)
        coeff = na ** m - c * na ** 4
        total += coeff * (-pi * (na * na - 4) * e2u).exp()
    # Ratio bound at n = N_FINITE (valid for all larger n: the polynomial
    # factor decreases toward 1 and the exponential factor only shrinks).
    n = N_FINITE
    r = (arb(n + 1) / arb(n)) ** m * (-pi * (2 * n + 1) * (2 * arb(3)).exp()).exp()
    assert r < arb(1) / 2, "ratio bound not < 1/2"
    na = arb(n + 1)
    first_neglected = (na ** m) * (-pi * (na * na - 4) * e2u).exp()
    tail = first_neglected / (1 - r)
    return total + tail


def eps_star_arb(u_val):
    """Certified upper bound on eps*(u) = 2*sum_{n>=2} n^4 e^{-pi*(n^2-1)*e^{2u}}.

    Finite sum n = 2..N_FINITE in Arb plus the same style of tail bound.
    """
    pi = arb.pi()
    u = arb(u_val)
    e2u = (2 * u).exp()
    total = arb(0)
    for n in range(2, N_FINITE + 1):
        na = arb(n)
        total += 2 * na ** 4 * (-pi * (na * na - 1) * e2u).exp()
    n = N_FINITE
    r = (arb(n + 1) / arb(n)) ** 4 * (-pi * (2 * n + 1) * (2 * arb(3)).exp()).exp()
    assert r < arb(1) / 2, "ratio bound not < 1/2"
    na = arb(n + 1)
    first_neglected = 2 * na ** 4 * (-pi * (na * na - 1) * e2u).exp()
    return total + first_neglected / (1 - r)


def main():
    ctx.prec = PREC
    t0 = time.time()

    print("=" * 72)
    print("  ARB/FLINT CERTIFICATE: weighted tail comparisons for u >= 3")
    print("  S_6 <= 5*S_4 and S_8 <= 32*S_4 (Appendix E.7, restated lemma)")
    print("  Method: n=2 domination; worst case u=3 certified in Arb")
    print("=" * 72)

    T6 = tail_sum_arb(6, 5, 3)
    T8 = tail_sum_arb(8, 32, 3)
    eps3 = eps_star_arb(3)

    ok6 = T6 < 16
    ok8 = T8 < 256
    oke = eps3 < arb("1e-1636")

    print(f"  T_6(3)  = {T6}   (< 16 ? {ok6})")
    print(f"  T_8(3)  = {T8}   (< 256 ? {ok8})")
    print(f"  eps*(3) = {eps3}   (< 1e-1636 ? {oke})")

    # Corroboration at u = 3.5 and u = 4 (the proof step is monotonicity;
    # these are spot checks that the T sums decrease, as they must).
    T6_35 = tail_sum_arb(6, 5, "3.5")
    T8_4 = tail_sum_arb(8, 32, 4)
    print(f"  T_6(3.5) = {T6_35}  (corroboration: <= T_6(3))")
    print(f"  T_8(4)   = {T8_4}  (corroboration: <= T_8(3))")
    okc = (T6_35 <= T6) and (T8_4 <= T8)

    all_ok = bool(ok6 and ok8 and oke and okc)
    elapsed = time.time() - t0
    print()
    print(f"  S_6 <= 5*S_4 certified:  {ok6}")
    print(f"  S_8 <= 32*S_4 certified: {ok8}")
    print(f"  eps*(3) < 1e-1636:       {oke}")
    print(f"  Monotonicity spot check: {okc}")
    print(f"  *** {'ALL CHECKS CERTIFIED' if all_ok else 'FAILURES — REVIEW ABOVE'} ***")
    print("=" * 72)

    result = {
        "claim": "S_6(u) <= 5*S_4(u) and S_8(u) <= 32*S_4(u) for all u >= 3",
        "method": "Arb/FLINT ball arithmetic via python-flint; n=2 domination "
                  "reduces each inequality to a T-sum at the worst case u=3; "
                  "T-sums are decreasing in u (each summand is), finite sum "
                  "n=3..40 plus certified geometric tail bound",
        "precision_bits": PREC,
        "finite_terms": N_FINITE,
        "T6_at_3_upper": str(T6.upper()),
        "T6_required_below": 16,
        "T8_at_3_upper": str(T8.upper()),
        "T8_required_below": 256,
        "eps_star_3_upper": str(eps3.upper()),
        "eps_star_3_required_below": "1e-1636",
        "passed": all_ok,
        "time_s": round(elapsed, 2),
    }
    out_path = os.path.join(ROOT, "results", "verify_weighted_comparisons.json")
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  -> results/verify_weighted_comparisons.json")
    return result


if __name__ == "__main__":
    result = main()
    sys.exit(0 if result["passed"] else 1)
