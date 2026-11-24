# Inventory — Jensen Grid Verification Release (E2), package contents

Every file in this package, with its size and SHA-256 hash.
Paths are package-relative. Verify any file by hashing it and
comparing with its row below.

| File | Bytes | SHA256 | Role |
|---|---|---|---|
| README.md | 18,711 | 38292aaa8349d9cd51ae76e73490b29eb7c2d51eb4fb22c42f5733d53d1554fb | Start here: claim, provenance, replay instructions |
| e2cert3.py | 25,346 | 5a080410456bb1f11e501cf8b0c739a9accc79ca904b79be9f8ef8a94ce428c6 | Operative certification script, round-4 corrected copy (both families, widening first-class; two header comments and two payload-metadata constructions corrected — see SOURCE-COMPARISON.md) |
| executed-original/e2cert3.py | 23,872 | 1e77032b48e6544b22633ddea293772d65aac0681ed291d040610e280da02db1 | The operative script's original bytes, as executed for the recorded full-grid replay (artifact of record, byte-untouched) |
| SOURCE-COMPARISON.md | 2,656 | 9139b4d02f73c38b1c0cc679a914a93c1f46792ea954399919aa0974ab59eada | Token-level comparison: executed original vs round-4 corrected operative copy, and the corrected copy's validation re-run |
| results_adj3.json | 21,667 | 981078c5f7f25177b9791f717e97af4772d7babbc770ce866aa5cd7a1d399769 | Certificate of record: both families 231/231 from W_k (generator: adj3_r3_1.py — see README) |
| results_e2cert3.json | 700,228 | 66a94e9a4e2464a677bfa013846eec08d1f0cd39f8db910d7bdf8e1c3942c51f | Recorded full-grid replay payload (generator: the executed-original bytes at executed-original/e2cert3.py; exit 0; see FULL-GRID-REPLAY.md) |
| results_e2cert3_smoke.json | 261,287 | 3fa786539b9966b6b2bcd9af87a77d27da85d26dee8e04a85c36174f67dbc14f | Smoke payload: control block 77/77 per family, exit 0 (partial coverage) |
| FULL-GRID-REPLAY.md | 8,815 | c0d302ceaae6d9c3b414c9adf7ada57edd587a1119dc1a8a46e92e09652461d1 | Replay record and field-by-field comparison vs the certificate of record (with the round-4 and round-5 notes; generator descriptions corrected in place, round 5) |
| SCHEMA-CROSSWALK.md | 9,045 | 7964bf86ebe42a7f7c552b19101b05db824c0e4199ebbc72041480471500e5da | Schema crosswalk: certificate vs replay payload (direct comparisons, with the radius rows split by extraction procedure; replay-only assertions, external expectations); bound-like float fields documented as diagnostics; replay-payload generator corrected to the executed-original script (round 5) |
| CERTIFICATE-ENDPOINTS.json | 20,583 | 3679f1ae5e61c8afc1ae174d084169da53fdb2f823faec0e00e321f41f536ac5 | Exact decimal endpoints of every ball in the certificate of record — all 68 balls (r*, 33 D_k, 33 failure-set deltas, W_0), computed from its printed ball strings (round-3 companion, W_0 added round 4; the certificate is byte-untouched) |
| DK-DERIVATION.md | 3,393 | 4e34b53463e47fde2cac5f8c257bc6553cb64c1c890d440135bf0439743c9af8 | The D_k tail-bound derivation, in full |
| CERTIFICATION.md | 12,865 | 670b366912f18c8e1583de9ff8775fc43f0cae75612c8ce3e1e550787515ce9f | Base-grid certification record with round-3 completion and R4-1 notes |
| EXTENSION.md | 15,323 | 1f7f072bd62758e7e15fe56927755d310ad5b582673eecf4dd068a2824e9b6ad | Extension-grid record, γ-run addendum, round-3 and R4-1 notes |
| REPORT.md | 6,249 | da05f9463f9909ae54ff4a0f7a4abfd41ade4bda57ec5184acf87049bdd971ae | Numerical precursor report (100-dps values; no enclosure claim) |
| adj2_sidechecks.json | 2,050 | 643e9a7f373e6fe2588044637c3dc3bf5c8299174131f8b882851070c99a2f80 | Side checks behind the original summation-tail bound |
| crosscheck-inputs/results_e2.json | 73,521 | 9d9d3815a6e882496e4ec4e5b754448b953200953d4af8ae9f0525090e96e97b | Cross-check input: E2 numerical payload moments (read by e2cert3.py at startup) |
| crosscheck-inputs/moments_partial.json | 3,842 | 591109b41fb16fc5e546a8a5bc58c846aa75f1e77f4619528f8c0cb83c673892 | Cross-check input: E1 moments (read by e2cert3.py at startup) |

## historical/ — artifacts of record, byte-untouched

Superseded as true-kernel certificates (see the README); retained
so the audit trail is replayable.

| File | Bytes | SHA256 | Role |
|---|---|---|---|
| historical/e2cert.py | 17,393 | 97ddfd8515aa67222ed4012f533abc705a72f78deeeff788a50c0dd62a7555f9 | First certification script (Borel family, base grid) |
| historical/e2cert2.py | 20,988 | c54d807933b3ba8ff0d9a290b23f6f3ce9c078a38e195dc6dd426997719d0282 | Extension script (Borel family, 231 grid); source of e2cert3.py's machinery |
| historical/e2gamma.py | 20,184 | 09850a704a0447e159b775c0567b7cc24145ddc092939c1e48b0868f558aed61 | Standard-grid consequence run script |
| historical/e2_jensen.py | 8,319 | f02d71cab2970097c5c8749b23ed874497b2ea10ae006f1926a180fe694432bf | Numerical precursor script |
| historical/results_e2gamma.json | 303,801 | 438d4e7f4b0f151aa99ba803099421e6be860c5e7ed98afff8a3906f0fa55acc | γ-run payload (superseded as a true-kernel certificate by results_adj3.json) |

## How this inventory is authenticated

This file deliberately does **not** carry a hash of itself, and
no self-hash would mean anything. The authentication route is
external: this inventory's own SHA-256 is recorded outside the
package, in the Spec 013 release record
(`specs/013-derivative-works/e2-release/INVENTORY.md` in the
program's records, which cites this file's hash), and at
publication the deposited archive's hash is recorded in the
deposit record (the Zenodo record for
`riemann-hypothesis-jensen-grid-verification-v1.zip`). Verify the
files against the table above; verify the table against one of
those external records.

## Not packaged (by design)

- `results_e2cert.json` and `results_e2cert2.json`: the original
  Borel payloads. They live only in the tower run dirs of record;
  their supersession is documented in CERTIFICATION.md and
  EXTENSION.md. The release does not depend on them.
- Tower run dirs (~/riemann-runs/e2cert-*, e2gamma-adj2-*,
  adj3-*, e2cert3-smoke-*, e2cert3-fullreplay-*): venues of
  record, cited in the method records, not redistributed.
