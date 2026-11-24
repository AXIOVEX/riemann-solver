# E2: Jensen polynomial hyperbolicity grid + E1 rung 2 (Newton cross-check).
# Spec 002 (conventional program), registered 2026-11-22 BEFORE computing.
#
# Definitions (paper's normalization, identical to E1 rung 1 / TC-02):
#   Phi(u) = sum_{n=1}^{6} (2 pi^2 n^4 e^{4.5u} - 3 pi n^2 e^{2.5u})
#                          * e^{-pi n^2 e^{2u}}
#   a_k = int_0^inf Phi(u) u^{2k} du / (2k)!   (Taylor coefficients of Xi;
#          quadrature on [0,8], tail double-exponentially small, as E1)
#   J^{d,n}(x) = sum_{k=0}^{d} C(d,k) a_{n+k} x^k
#
# REGISTERED GRID (fixed here, before any computation; NOT extended after
# seeing results): d = 2..8, n = 0..10  ->  63 polynomials.
# CONTROLS (run first): J^{2,0}, J^{3,0}, J^{3,1}, J^{4,0} must be
# hyperbolic, reproducing TC-02 control C-D, or the batch is NO VERDICT.
#
# NEWTON INEQUALITIES -- exact form implemented:
#   For J(x) = sum_{k=0}^{d} c_k x^k with c_k = C(d,k) a_{n+k}, define the
#   binomial-normalized coefficients b_k = c_k / C(d,k) = a_{n+k}.
#   Newton: for 1 <= k <= d-1,
#       b_k^2 - b_{k-1} * b_{k+1} >= 0.
#   Reported per polynomial: the minimum over k of the NORMALIZED margin
#       (b_k^2 - b_{k-1} b_{k+1}) / b_k^2  =  1 - a_{n+k-1} a_{n+k+1} / a_{n+k}^2
#   together with the raw margin and the argmin k.
#
# HYPERBOLICITY: roots by mpmath.polyroots at 100 dps (extraprec 50).
#   Reported per polynomial: max |Im(root)|, min pairwise root separation,
#   max residual |J(root)|. NUMERICAL evidence at stated precision only
#   (Spec 002 R1); this is NOT a Sturm/interval certification.
#
# FALSIFIER THRESHOLDS (fixed here): relative max|Im| > 1e-40 (relative to
#   max(1, max|root|)), or normalized Newton margin < -1e-40, is a FALSIFIER
#   CANDIDATE: flagged in the payload, grid NOT extended, interval
#   certification required before any claim.
#
# MOMENT CROSS-CHECK: a_0..a_18 recomputed fresh at 100 dps and compared
#   against the E1 payload (~/riemann-runs/e1-20261122/moments_partial.json,
#   80 dps); max relative deviation reported.
#
# Payload: results_e2.json in the run directory. The verdict comes from
# this file, never from the wrapper.
import json, sys, os
import mpmath as mp

mp.mp.dps = 100

RUN_DIR = os.path.dirname(os.path.abspath(__file__))
E1_MOMENTS = os.path.expanduser("~/riemann-runs/e1-20261122/moments_partial.json")

IMAG_RTOL_FLOOR = mp.mpf("1e-40")
NEWTON_MARGIN_FLOOR = mp.mpf("-1e-40")

def Phi(u):
    s = mp.mpf("0")
    for n in range(1, 7):
        s += (2*mp.pi**2 * n**4 * mp.e**(mp.mpf("4.5")*u)
              - 3*mp.pi * n**2 * mp.e**(mp.mpf("2.5")*u)) * mp.e**(-mp.pi * n**2 * mp.e**(2*u))
    return s

def s(x):
    return mp.nstr(x, 40, strip_zeros=False)

# ---- moments a_0..a_18, fresh at 100 dps ---------------------------------
KMAX = 18
a = []
for k in range(0, KMAX + 1):
    m = mp.quad(lambda u: Phi(u) * u**(2*k), [0, 8]) / mp.factorial(2*k)
    a.append(m)
    print(f"moment a_{k} done", flush=True)

# ---- cross-check against E1 payload --------------------------------------
cross = {"available": False}
try:
    with open(E1_MOMENTS) as f:
        e1 = json.load(f)["a"]
    devs = []
    for k in range(0, KMAX + 1):
        ref = mp.mpf(e1[k])
        devs.append(abs(a[k] - ref) / abs(ref))
    cross = {
        "available": True,
        "source": E1_MOMENTS,
        "max_rel_dev": s(max(devs)),
        "argmax_k": int(max(range(len(devs)), key=lambda i: devs[i])),
        "per_k_rel_dev": {str(k): s(devs[k]) for k in range(len(devs))},
    }
except Exception as e:
    cross = {"available": False, "error": str(e)}
print("cross-check:", cross.get("max_rel_dev", cross), flush=True)

# ---- Jensen machinery -----------------------------------------------------
def jensen(d, n):
    c = [mp.binomial(d, k) * a[n + k] for k in range(d + 1)]  # c_k, ascending
    # Newton margins on b_k = a_{n+k}
    margins = []
    for k in range(1, d):
        b0, b1, b2 = a[n + k - 1], a[n + k], a[n + k + 1]
        raw = b1**2 - b0 * b2
        margins.append((raw / b1**2, raw, k))
    nm_min = min(margins, key=lambda t: t[0])
    # roots: polyroots wants coefficients highest degree first
    coeffs_desc = [c[k] for k in range(d, -1, -1)]
    roots = mp.polyroots(coeffs_desc, maxsteps=1000, extraprec=50)
    scale = max(mp.mpf(1), max(abs(r) for r in roots))
    max_im = max(abs(mp.im(r)) for r in roots)
    min_sep = min(abs(roots[i] - roots[j])
                  for i in range(len(roots)) for j in range(i + 1, len(roots)))
    def Jeval(x):
        return sum(c[k] * x**k for k in range(d + 1))
    max_resid = max(abs(Jeval(r)) for r in roots)
    hyp = bool(max_im / scale <= IMAG_RTOL_FLOOR)
    newt_ok = bool(nm_min[0] >= NEWTON_MARGIN_FLOOR)
    return {
        "d": d, "n": n,
        "hyperbolic_numerical": hyp,
        "max_abs_imag": s(max_im),
        "max_rel_imag": s(max_im / scale),
        "root_scale": s(scale),
        "min_root_separation": s(min_sep),
        "max_residual": s(max_resid),
        "roots_real_parts": [s(mp.re(r)) for r in roots],
        "newton_min_normalized_margin": s(nm_min[0]),
        "newton_min_raw_margin": s(nm_min[1]),
        "newton_argmin_k": nm_min[2],
        "newton_ok": newt_ok,
        "falsifier_candidate": bool((not hyp) or (not newt_ok)),
    }

CONTROLS = [(2, 0), (3, 0), (3, 1), (4, 0)]
GRID = [(d, n) for d in range(2, 9) for n in range(0, 11)]

results = {}
controls_ok = True
for (d, n) in CONTROLS:
    r = jensen(d, n)
    results[f"{d},{n}"] = r
    print(f"control J^{d},{n}: hyp={r['hyperbolic_numerical']} "
          f"maxIm={r['max_abs_imag']} newton={r['newton_min_normalized_margin']}",
          flush=True)
    if not r["hyperbolic_numerical"]:
        controls_ok = False

verdict = "CONTROLS PASS"
if not controls_ok:
    verdict = "NO VERDICT (control failure)"

falsifiers = []
if controls_ok:
    for (d, n) in GRID:
        if f"{d},{n}" in results:
            continue
        r = jensen(d, n)
        results[f"{d},{n}"] = r
        print(f"J^{d},{n}: hyp={r['hyperbolic_numerical']} "
              f"maxIm={r['max_abs_imag']} newton={r['newton_min_normalized_margin']}",
              flush=True)
    falsifiers = [k for k, r in results.items() if r["falsifier_candidate"]]
    if falsifiers:
        verdict = "FALSIFIER CANDIDATE(S) PRESENT"
    else:
        verdict = "ALL HYPERBOLIC ON GRID (NUMERICAL, 100 dps)"

all_hyp = all(r["hyperbolic_numerical"] for r in results.values())
all_newt = all(r["newton_ok"] for r in results.values())
worst_im = max(results.values(), key=lambda r: mp.mpf(r["max_rel_imag"]))
worst_nm = min(results.values(), key=lambda r: mp.mpf(r["newton_min_normalized_margin"]))
smallest_sep = min(results.values(), key=lambda r: mp.mpf(r["min_root_separation"]))

out = {
    "experiment": "E2 Jensen hyperbolicity grid + E1 rung 2 Newton cross-check",
    "spec": "002-conventional-program",
    "classification": "NUMERICAL (mpmath floating point at 100 dps; NOT a Sturm/interval certification)",
    "dps": 100,
    "grid_registered": {"d": "2..8", "n": "0..10", "count": 63},
    "controls": [f"{d},{n}" for (d, n) in CONTROLS],
    "controls_ok": controls_ok,
    "verdict": verdict,
    "moments_a": {str(k): s(a[k]) for k in range(len(a))},
    "moments_cross_check_vs_E1": cross,
    "all_hyperbolic": bool(all_hyp),
    "all_newton_ok": bool(all_newt),
    "falsifier_candidates": falsifiers,
    "worst_max_rel_imag": {"key": f"{worst_im['d']},{worst_im['n']}",
                           "value": worst_im["max_rel_imag"]},
    "smallest_newton_normalized_margin": {"key": f"{worst_nm['d']},{worst_nm['n']}",
                                          "value": worst_nm["newton_min_normalized_margin"]},
    "smallest_root_separation": {"key": f"{smallest_sep['d']},{smallest_sep['n']}",
                                 "value": smallest_sep["min_root_separation"]},
    "polynomials": results,
}
with open(os.path.join(RUN_DIR, "results_e2.json"), "w") as f:
    json.dump(out, f, indent=2)
print("VERDICT:", verdict)
print(json.dumps({k: out[k] for k in
                  ("controls_ok", "all_hyperbolic", "all_newton_ok",
                   "falsifier_candidates", "worst_max_rel_imag",
                   "smallest_newton_normalized_margin",
                   "smallest_root_separation")}, indent=2))
sys.exit(0)
