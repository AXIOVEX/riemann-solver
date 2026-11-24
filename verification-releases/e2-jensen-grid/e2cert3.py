# E2 CERTIFICATION -- OPERATIVE SCRIPT (e2cert3.py, 2026-11-23).
#
# ROUND-4 REVISION (2026-11-24, Spec 013 round-4 adjudication,
# items E2-R4-1..3): this copy corrects two header comments (the
# STEP 1 widening placement, and the noncontainment scope
# sentence below) and two payload-metadata constructions (the
# venue field now records the actual execution host instead of a
# hardcoded string, and the selftest grid metadata now describes
# the actual seven-point grid). No certification code, setting,
# or numeric path changed. The script as executed for the
# recorded full-grid replay (sha256 1e77032b...) is preserved
# byte-identical at executed-original/e2cert3.py in this package;
# see SOURCE-COMPARISON.md for the exact token-level comparison,
# and coscientist-adjudication-round4.md for the validation
# re-run of this corrected copy.
#
# This is the operative E2 certification script for any future run.
# It supersedes, as the run-of-record going forward, the historical
# scripts e2cert.py, e2cert2.py and e2gamma.py, which are ARTIFACTS
# OF RECORD and remain byte-untouched (their hashes are citations of
# record: e2cert.py sha256 97ddfd85..., e2cert2.py sha256
# c54d8079..., e2gamma.py sha256 09850a70...). The Sturm / ball
# machinery below (Phi_acb, tail_ball, poly_rem, discriminant_ball,
# sturm_attempt and helpers) is carried over from e2cert2.py
# unchanged in method; the gamma family coefficients are those of
# e2gamma.py (gamma_k = k! * a_k, the standard xi Jensen family).
# The Borel family (coefficients a_k) and the standard gamma family
# are both certified here, on both registered grids:
#   base grid      d = 2..8,  n = 0..10   (77 points, the control
#                                          block)
#   extension grid d = 2..12, n = 0..20   (231 points; includes the
#                                          base grid)
# Definitions (paper's normalization, identical to E1/E2/E2-cert):
#   Phi(u) = sum_{n>=1} (2 pi^2 n^4 e^{4.5u} - 3 pi n^2 e^{2.5u})
#                          * e^{-pi n^2 e^{2u}}
#   a_k = int_0^inf Phi(u) u^{2k} du / (2k)!
#   Borel  J^{d,n}(x) = sum_{k=0}^{d} C(d,k) a_{n+k} x^k
#   gamma  J^{d,n}(x) = sum_{k=0}^{d} C(d,k) gamma_{n+k} x^k,
#          gamma_k = k! * a_k.
#
# THE SIX-TERM TRUNCATION AND THE n >= 7 TAIL -- CONTAINMENT NOW
# PROVED, NOT ASSERTED. The callback Phi_acb evaluates the finite
# six-term truncation sum (n = 1..6), which is entire in u; the full
# Phi series is analytic only in the strip |Im u| < pi/4 in this
# e^{2u} convention. The historical scripts enclosed only the
# truncated moments I_k = (int_0^8 Phi_acb(u) u^{2k} du + T) / (2k)!,
# where T is the integration-domain tail ball for u >= 8 (tail_ball()
# below), and contained NO summation-index tail. The round-3
# adjudication (coscientist-adjudication-round3.md, R3-1 SUSTAINED,
# section 1) proved position-independent noncontainment at
# k = 0, 1, 2; the test does not establish containment in the old
# balls at the remaining indices. The
# omitted contributions delta a_k (Arb-enclosed at 512 bits) exceed
# the accepted balls' diameters there (delta a_0 = 2.148808e-65 vs
# rad(I_0) = 5.421e-76), so the true moments lay OUTSIDE I_0, I_1,
# I_2. The repair, executed in that adjudication and made first-class
# here, widens every moment ball by the proved n >= 7 tail bound:
#   W_k = I_k + [0, D_k],
#   D_k = 2 (49 pi)^2 e^{-49 pi}
#         / [ (98 pi - 9/2)^{2k+1} (1 - r*) ],
#   r*  = (8/7)^4 e^{-15 pi} = 5.8387849549...e-21,
# derived in coscientist-adjudication-round3.md section 1.5 (steps
# 1-5: term majorant 2 q_n^2 e^{-q_n} e^{-(2 q_n - 9/2) u} with
# q_n = pi n^2, integrated against u^{2k} over [0, inf), successive-
# term ratio at most r*, geometric sum from n = 7) and evaluated in
# Arb there: D_0 = 2.184234e-65, D_1 = 2.373212e-70,
# D_2 = 2.578540e-75, with delta a_k <= D_k confirmed against the
# independent Arb enclosures for every k = 0..32. Since
# trunc(a_k) in I_k by the acb_integral enclosure and
# 0 <= delta a_k <= D_k by the bound, the TRUE a_k lies in W_k for
# every k: containment is proved, not asserted. The widening code
# below (widening_bounds) is promoted verbatim from the round-3
# adjudication's executed script adj3_r3_1.py (tower run dir
# ~/riemann-runs/adj3-20261123T0725, sha256 3f9881f5...), whose
# widened re-runs are the certificate of record: payload
# results_adj3.json (sha256 981078c5...), Borel 231/231 and standard
# gamma 231/231 CERTIFIED from W_k. Spec 006
# (specs/006-round3-repair-pass/REPORT.md, section 5 flag) recorded
# the remaining hygiene item -- making this widening a first-class
# part of the certification code -- which this script closes.
#
# STEP 1 -- INTERVAL MOMENTS, WIDENED. a_0..a_32 enclosed as in
# e2cert2.py: (a) int_0^8 Phi(u) u^{2k} du by FLINT's acb_integral
# (adaptive Gauss-Legendre in ball arithmetic; guaranteed enclosure
# of the integral of the six-term truncation as evaluated); (b) the
# u >= 8 integration tail added as tail_ball() (analytic bound,
# ~10^{-1.21e7}, derivation in e2cert2.py's header, unchanged);
# (c) the n >= 7 summation tail added as the ball [0, D_k] above:
# the integration-domain tail in (b) is added before dividing by
# (2k)!, and the normalized summation-tail enclosure [0, D_k]
# (D_k bounds the contribution to a_k, already normalized by
# (2k)!) is added afterward, after the division (see
# compute_moments_widened below).
#
# STEP 2 -- STURM CERTIFICATION WITH BALL COEFFICIENTS, identical
# method to e2cert2.py / e2gamma.py: fixed straight-line Sturm
# sequence in ball arithmetic; count == d plus a separable final
# constant => CERTIFIED; count < d with all signs determinate =>
# FAILED (falsifier); any indeterminate sign => INDETERMINATE at
# that level. Discriminant ball (Sylvester determinant) recorded per
# polynomial as the margin dataset; a determinate disc <= 0 with a
# Sturm count of d demotes to INDETERMINATE (internal inconsistency).
#
# PRECISION ESCALATION: levels 256, 512, 1024, 2048 bits (cap 2048),
# identical ladder to the historical runs. Moments (and D_k) are
# recomputed at each level. Every attempt is logged.
#
# CROSS-CHECKS vs numerical records (as in e2cert2.py): the
# unwidened enclosures I_k are checked for intersection against the
# E2 payload moments (40 significant digits) and the E1 moments
# (80 dps); consistency = the stored ball intersects I_k.
#
# MODES:
#   python3 e2cert3.py           full run: extension grid, both
#                                families -> results_e2cert3.json
#   python3 e2cert3.py smoke     control block only (77 points,
#                                both families)
#                                -> results_e2cert3_smoke.json
#   python3 e2cert3.py selftest  7-point spot grid, level 0 only
#                                -> results_e2cert3_selftest.json
# The verdict comes from the payload file, never from the wrapper
# (Spec 002 TC-03).

import json, math, os, socket, sys, time, traceback
from flint import arb, acb, arb_mat, ctx

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
E2_PAYLOAD = os.path.expanduser(
    "~/riemann-runs/e2-20261122T143947/results_e2.json")
E1_MOMENTS = os.path.expanduser(
    "~/riemann-runs/e1-20261122/moments_partial.json")
ADJ3_PAYLOAD = os.path.join(RUN_DIR, "results_adj3.json")
KMAX = 32
GRID = [(d, n) for d in range(2, 13) for n in range(0, 21)]
CONTROL_BLOCK = [(d, n) for d in range(2, 9) for n in range(0, 11)]
CONTROLS_TC02 = [(2, 0), (3, 0), (3, 1), (4, 0)]
LEVELS = [256, 512, 1024, 2048]
FAMILIES = ("borel", "gamma")

def Phi_acb(u):
    s = acb(0)
    pi = acb.pi()
    for n in range(1, 7):
        s += (2 * pi**2 * n**4 * (acb("4.5") * u).exp()
              - 3 * pi * n**2 * (acb("2.5") * u).exp()) \
            * (-pi * n**2 * (2 * u).exp()).exp()
    return s

def tail_ball():
    """Ball T with [0, 4550 pi^2 e^{g(8)} / lambda] subset of T, plus
    the bound ball itself (see header derivation, KMAX = 32:
    g(u) = (9/2) u + 64 ln u - pi e^{2u}, lambda = 2 pi e^{16} - 12.5,
    g(8) = 36 + 64 ln 8 - pi e^{16})."""
    pi = arb.pi()
    e16 = arb(16).exp()
    g8 = arb(36) + arb(2 * KMAX) * arb(8).log() - pi * e16
    lam = 2 * pi * e16 - (arb(9) / 2 + arb(2 * KMAX) / 8)
    assert lam > 0
    bound = arb(4550) * pi**2 * g8.exp() / lam
    assert bound > 0
    T = arb(0).union(bound)  # contains 0 and all of bound => contains [0, sup]
    return T, bound


def widening_bounds():
    """The proved n >= 7 summation-tail bounds D_k, k = 0..KMAX, in
    Arb at the current precision, together with r*. PROMOTED
    VERBATIM (formula and evaluation) from the round-3
    adjudication's executed script adj3_r3_1.py, Part 2 (tower run
    dir ~/riemann-runs/adj3-20261123T0725, sha256 3f9881f5...);
    derivation: coscientist-adjudication-round3.md section 1.5.
    D_k bounds delta a_k, the full omitted n >= 7 contribution to
    a_k, for every k."""
    pi = arb.pi()
    r_star = (arb(8) / 7) ** 4 * (-15 * pi).exp()
    Dk = []
    for k in range(KMAX + 1):
        D = (2 * (49 * pi) ** 2 * (-49 * pi).exp()
             / ((98 * pi - arb(9) / 2) ** (2 * k + 1) * (1 - r_star)))
        Dk.append(D)
    return Dk, r_star


def compute_moments_widened(prec):
    """Interval moments at precision prec: the unwidened enclosures
    I_k exactly as e2cert2.py's compute_moments (six-term truncation
    integral over [0,8] plus the u >= 8 tail ball, divided by (2k)!),
    the widening bounds D_k, and the widened balls
    W_k = I_k + [0, D_k] (the ball [0, sup D_k] added in ball
    arithmetic; equivalent containment to the adjudication's
    endpoint construction W = [lo(I), hi(I) + sup D_k])."""
    ctx.prec = prec
    rel = arb(2) ** (-(prec - 64))
    abst = arb(2) ** (-(prec - 64)) * arb("1e-12")
    T, bound = tail_ball()
    a = []
    for k in range(KMAX + 1):
        def f(u, analytic, k=k):
            return Phi_acb(u) * u ** (2 * k)
        r = acb.integral(f, acb(0), acb(8), rel_tol=rel, abs_tol=abst,
                         eval_limit=10**7, depth_limit=1000)
        assert r.imag.contains(arb(0)) or r.imag == 0
        a.append((r.real + T) / arb(math.factorial(2 * k)))
        print(f"  moment a_{k} done (prec {prec})", flush=True)
    Dk, r_star = widening_bounds()
    W = [a[k] + arb(0).union(Dk[k].upper()) for k in range(KMAX + 1)]
    return a, W, Dk, r_star, bound


def family_coeffs(family, W, d, n):
    """Coefficient balls of J^{d,n} from the widened moment balls W:
    borel: C(d,k) * W[n+k]; gamma: C(d,k) * (n+k)! * W[n+k]
    (e2gamma.py's coefficient construction)."""
    if family == "borel":
        return [arb(math.comb(d, k)) * W[n + k]
                for k in range(d + 1)]
    return [arb(math.comb(d, k)) * arb(math.factorial(n + k)) * W[n + k]
            for k in range(d + 1)]


def sign_of(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0

def poly_rem(A, B):
    """Enclosure of A mod B (lists of arb, ascending). Requires the
    leading coefficient ball of B to exclude 0 (checked by caller)."""
    dB = len(B) - 1
    R = list(A)
    if len(R) < len(B):
        R = R + [arb(0)] * (len(B) - len(R))
    leadB = B[-1]
    for k in range(len(A) - 1 - dB, -1, -1):
        factor = R[dB + k] / leadB
        for j in range(dB + 1):
            R[k + j] = R[k + j] - factor * B[j]
    return R[:dB]

def discriminant_ball(coeffs):
    """disc(p) as an arb ball via the Sylvester determinant in ball
    arithmetic. Returns the ball, or None if unavailable."""
    d = len(coeffs) - 1
    if d < 2:
        return None
    try:
        P = [coeffs[d - i] for i in range(d + 1)]            # descending
        Q = [(d - j) * coeffs[d - j] for j in range(d)]      # p' descending
        size = 2 * d - 1
        rows = []
        for i in range(d - 1):
            rows.append([arb(0)] * i + P + [arb(0)] * (d - 2 - i))
        for i in range(d):
            rows.append([arb(0)] * i + Q + [arb(0)] * (d - 1 - i))
        res = arb_mat(rows).det()
        disc = res / coeffs[-1]
        if (d * (d - 1) // 2) % 2 == 1:
            disc = -disc
        return disc
    except Exception:
        return None

def disc_record(disc):
    """Serializable record of a discriminant ball: sign, and the
    margin dataset (mid, radius, relative width, sign margin =
    (|mid| - rad)/|mid|, log10 |mid|)."""
    if disc is None:
        return {"sign": 0, "available": False}
    rec = {"sign": sign_of(disc), "available": True, "ball": str(disc)}
    try:
        mid = disc.mid()
        rad = disc.rad()
        rec["mid"] = float(mid)
        rec["rad"] = float(rad)
        if mid != 0:
            rec["rel_width"] = float(rad / abs(mid))
            rec["sign_margin"] = float((abs(mid) - rad) / abs(mid))
            rec["log10_abs_mid"] = float(
                abs(mid).log() / arb(10).log())
    except Exception:
        pass
    return rec

def sturm_attempt(coeffs):
    """One Sturm certification attempt at the current precision.
    coeffs: ascending arb list for J^{d,n}. Identical logic to the
    cert run, plus the discriminant ball record."""
    d = len(coeffs) - 1
    P0 = list(coeffs)
    P1 = [(j + 1) * coeffs[j + 1] for j in range(d)]
    for P in (P0, P1):
        if sign_of(P[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "leading coefficient sign of J or J'"}
    seq = [P0, P1]
    while True:
        A, B = seq[-2], seq[-1]
        if len(B) == 1:
            break
        if sign_of(B[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "divisor leading coefficient sign"}
        R = poly_rem(A, B)
        if len(R) == 0:
            return {"status": "indeterminate",
                    "reason": "empty remainder"}
        if sign_of(R[-1]) == 0:
            return {"status": "indeterminate",
                    "reason": "remainder leading coefficient sign",
                    "sturm_length_so_far": len(seq)}
        seq.append([-c for c in R])
        if len(R) == 1:
            break
    final = seq[-1]
    if len(final) != 1 or sign_of(final[0]) == 0:
        return {"status": "indeterminate",
                "reason": "final Sturm constant sign (separability)"}

    def variations(at_minus):
        signs = []
        for P in seq:
            s = sign_of(P[-1])
            if s == 0:
                return None
            if at_minus and ((len(P) - 1) % 2 == 1):
                s = -s
            signs.append(s)
        return sum(1 for i in range(len(signs) - 1)
                   if signs[i] != signs[i + 1])

    vm, vp = variations(True), variations(False)
    if vm is None or vp is None:
        return {"status": "indeterminate",
                "reason": "sign variation evaluation"}
    count = vm - vp
    drec = disc_record(discriminant_ball(coeffs))
    rec = {"sturm_length": len(seq), "real_root_count": count,
           "final_sturm_constant": str(final[0]),
           "discriminant_sign": drec["sign"],
           "discriminant": drec}
    if count == d:
        if drec["sign"] < 0:
            rec.update({"status": "indeterminate",
                        "reason": "INTERNAL INCONSISTENCY: Sturm count "
                                  "= d but discriminant <= 0"})
            return rec
        rec["status"] = "certified"
        return rec
    if count < d:
        rec["status"] = "failed"
        return rec
    rec.update({"status": "indeterminate",
                "reason": f"root count {count} > degree {d} (impossible; "
                          "enclosure failure)"})
    return rec

def stored_ball(s):
    mant = s.lower().split('e')[0].replace('-', '').replace('+', '') \
        .replace('.', '').lstrip('0')
    sig = len(mant) if mant else 1
    v = arb(s)
    rho = v * arb(10) ** (-(sig - 1)) / 2
    return (v - rho).union(v + rho)

def ballrec(x):
    return {"ball": str(x), "mid": float(x.mid()), "rad": float(x.rad()),
            "rel_width": float(x.rad() / abs(x.mid())) if x.mid() != 0
            else None}

def tally(recs):
    return {
        "certified": sum(1 for r in recs if r["status"] == "certified"),
        "failed": sum(1 for r in recs if r["status"] == "failed"),
        "indeterminate": sum(1 for r in recs
                             if r["status"] == "indeterminate"),
        "total": len(recs),
    }


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    assert mode in ("full", "smoke", "selftest"), mode
    t_start = time.time()
    with open(E2_PAYLOAD) as fh:
        e2 = json.load(fh)
    e2_moments = e2["moments_a"]
    with open(E1_MOMENTS) as fh:
        e1_moments = json.load(fh)["a"]

    if mode == "selftest":
        grids = {fam: [(2, 0), (3, 0), (3, 1), (4, 0), (8, 10),
                       (12, 0), (12, 20)] for fam in FAMILIES}
        levels = LEVELS[:1]
    elif mode == "smoke":
        grids = {fam: list(CONTROL_BLOCK) for fam in FAMILIES}
        levels = LEVELS
    else:
        grids = {fam: list(GRID) for fam in FAMILIES}
        levels = LEVELS

    decided = {fam: {} for fam in FAMILIES}
    attempts = {fam: [] for fam in FAMILIES}
    level_info = []
    moments_I = None
    moments_W = None
    Dk_final = None
    tail_bound_rec = None

    for level, prec in enumerate(levels):
        if all(all(pt in decided[fam] for pt in grids[fam])
               for fam in FAMILIES):
            break
        t0 = time.time()
        print(f"LEVEL {level}: prec {prec} bits", flush=True)
        a, W, Dk, r_star, bound = compute_moments_widened(prec)
        moments_I, moments_W, Dk_final = a, W, Dk
        tail_bound_rec = {"ball": str(bound),
                          "log10_approx": float(bound.log()
                                                / arb(10).log())}
        # cross-checks on the UNWIDENED enclosures I_k, as in
        # e2cert2.py: k=0..18 vs E2 (40 sig digits) AND E1 (80 dps);
        # k=19..32 vs E1 only.
        cross_e2, cross_e1 = {}, {}
        max_rel_width = 0.0
        for k in range(KMAX + 1):
            if k <= 18:
                cross_e2[str(k)] = bool(
                    a[k].overlaps(stored_ball(e2_moments[str(k)])))
            cross_e1[str(k)] = bool(
                a[k].overlaps(stored_ball(e1_moments[k])))
            rw = float(a[k].rad() / abs(a[k].mid()))
            max_rel_width = max(max_rel_width, rw)
        level_info.append({
            "level": level, "prec_bits": prec,
            "seconds": time.time() - t0,
            "max_moment_rel_width_unwidened": max_rel_width,
            "cross_check_vs_E2_all_intersect": all(cross_e2.values()),
            "cross_check_vs_E2": cross_e2,
            "cross_check_vs_E1_all_intersect": all(cross_e1.values()),
            "cross_check_vs_E1": cross_e1,
        })
        print(f"  cross-checks: E2 {all(cross_e2.values())}, "
              f"E1 {all(cross_e1.values())}, "
              f"max rel width (unwidened) {max_rel_width:.3e}",
              flush=True)
        for fam in FAMILIES:
            for (d, n) in grids[fam]:
                if (d, n) in decided[fam]:
                    continue
                coeffs = family_coeffs(fam, W, d, n)
                try:
                    rec = sturm_attempt(coeffs)
                except Exception as e:
                    rec = {"status": "indeterminate",
                           "reason": f"exception: {e}"}
                rec.update({"d": d, "n": n, "level": level,
                            "prec_bits": prec, "family": fam})
                attempts[fam].append(rec)
                print(f"  {fam} J^{d},{n}: {rec['status']}"
                      + (f" ({rec.get('reason')})"
                         if rec["status"] != "certified" else ""),
                      flush=True)
                if rec["status"] in ("certified", "failed"):
                    decided[fam][(d, n)] = rec

    fam_out = {}
    for fam in FAMILIES:
        records = {}
        for (d, n) in grids[fam]:
            if (d, n) in decided[fam]:
                records[f"{d},{n}"] = decided[fam][(d, n)]
            else:
                last = [r for r in attempts[fam]
                        if r["d"] == d and r["n"] == n][-1]
                last = dict(last)
                last["status"] = "indeterminate"
                last["cap_bits"] = levels[-1]
                records[f"{d},{n}"] = last
        entry = {"tallies": tally(list(records.values())),
                 "polynomials": records}
        if mode == "full":
            ctrl = [records[f"{d},{n}"] for (d, n) in CONTROL_BLOCK]
            entry["control_block_77"] = tally(ctrl)
            entry["control_block_all_certified"] = all(
                r["status"] == "certified" for r in ctrl)
            tt = entry["tallies"]
            if not entry["control_block_all_certified"]:
                entry["label_outcome"] = (
                    "NO VERDICT: the 77-point control block did not "
                    "certify in full from the widened balls")
            elif tt["failed"] > 0:
                failed_named = sorted(
                    k for k, r in records.items()
                    if r["status"] == "failed")
                entry["label_outcome"] = (
                    f"FALSIFIER: {tt['failed']} polynomial(s) proved "
                    f"NOT hyperbolic ({', '.join(failed_named)}); "
                    f"{tt['certified']} CERTIFIED, "
                    f"{tt['indeterminate']} INDETERMINATE")
            elif tt["certified"] == len(GRID):
                entry["label_outcome"] = (
                    f"CERTIFIED for the extended grid d=2..12, "
                    f"n=0..20 (231/231), {fam} family, from the "
                    f"widened balls W_k = I_k + [0, D_k]")
            else:
                entry["label_outcome"] = (
                    f"CERTIFIED for a subset ({tt['certified']}/231); "
                    f"INDETERMINATE for {tt['indeterminate']}")
        fam_out[fam] = entry

    out = {
        "experiment": "E2 certification, OPERATIVE script e2cert3.py "
                      "(Sturm/interval, widened moment balls "
                      "W_k = I_k + [0, D_k]), Spec 002 R5; widening "
                      "per coscientist-adjudication-round3.md R3-1 "
                      "and Spec 006",
        "classification": "CERTIFIED (ball arithmetic) where stated; "
                          "see per-family label_outcome / tallies",
        "venue": socket.gethostname(),
        "mode": mode,
        "families": list(FAMILIES),
        "grid_registered": ({"d": "2..12", "n": "0..20", "count": 231}
                            if mode == "full" else
                            {"control_block": {"d": "2..8",
                                               "n": "0..10",
                                               "count": 77}}
                            if mode == "smoke" else
                            {"selftest_points":
                                 [f"{d},{n}" for (d, n)
                                  in grids["borel"]],
                             "count": len(grids["borel"]),
                             "levels": "level 0 only (256 bits)"}),
        "control_block": {"d": "2..8", "n": "0..10", "count": 77},
        "controls_tc02": [f"{d},{n}" for (d, n) in CONTROLS_TC02],
        "widening": {
            "formula": "W_k = I_k + [0, D_k]; D_k = 2 (49 pi)^2 "
                       "e^{-49 pi} / ((98 pi - 9/2)^{2k+1} "
                       "(1 - r*)); r* = (8/7)^4 e^{-15 pi}",
            "source": "coscientist-adjudication-round3.md s1.5; "
                      "code promoted verbatim from adj3_r3_1.py "
                      "(sha256 3f9881f5...)",
            "certificate_of_record": "results_adj3.json (sha256 "
                                     "981078c5...)",
            "D_k": {str(k): ballrec(Dk_final[k])
                    for k in range(KMAX + 1)},
            "r_star_sup": float(r_star.upper()),
        },
        "tail_bound_u_ge_8": tail_bound_rec,
        "levels": level_info,
        "moments_unwidened_final_level": {
            str(k): ballrec(moments_I[k]) for k in range(KMAX + 1)},
        "moments_widened_final_level": {
            str(k): ballrec(moments_W[k]) for k in range(KMAX + 1)},
        "results_by_family": fam_out,
        "attempts_by_family": attempts,
        "seconds_total": time.time() - t_start,
    }
    name = {"full": "results_e2cert3.json",
            "smoke": "results_e2cert3_smoke.json",
            "selftest": "results_e2cert3_selftest.json"}[mode]
    with open(os.path.join(RUN_DIR, name), "w") as fh:
        json.dump(out, fh, indent=2)
    for fam in FAMILIES:
        print(f"{fam.upper()} TALLIES:", fam_out[fam]["tallies"])
        if "label_outcome" in fam_out[fam]:
            print(f"{fam.upper()} LABEL:", fam_out[fam]["label_outcome"])
    print("PAYLOAD:", name)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        err = traceback.format_exc()
        print(err, flush=True)
        try:
            with open(os.path.join(RUN_DIR,
                                   "results_e2cert3_ERROR.json"),
                      "w") as fh:
                json.dump({"error": err}, fh, indent=2)
        except Exception:
            pass
        sys.exit(1)

