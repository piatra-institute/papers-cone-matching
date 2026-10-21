# Audit

Dated log of editorial passes and verification runs. Newest first.
See the workspace docs (run `papers docs`): writing-pipeline.md §7 and refresh-pipeline.md.

## 2026-10-21 — publication

Status set to published and PDF rebuilt (zero missing-character warnings). All 17 pages rendered at 100 dpi are pixel-identical to the inspected build; visual-review.json moved from 8cc4c7de2e4e to 3c63ecf046c5 (rebuilt twice; the second rebuild followed a web_blurb edit using subscripts, and was also pixel-identical). Manuscript unchanged, so the claims ledger and editorial review still cover it. PDF synced to the website repository (public/papers), ownPapers entry added to app/papers/page.tsx (the website repository is left uncommitted for the editor), published-stage check run, GitHub repository created and initial signed commit pushed.

## 2026-10-21 — polish pass (figure legend, Table 2 placement, reference subscripts)

Scope: polish fixes requested in the main review.
Changes:
  - simulation/figures.py: Figure 2(b) legend "automorphic support" / "spectral support with the integrality condition of Proposition 5" (two lines); recorded run repeated (verification/model.json; previous receipts archived), 26 checks passed, no numerical change.
  - paper/PAPER.md: Table 2 caption shortened; code derivation moved to the text after the table.
  - references.yaml: inline math in titles (Ginzburg GL_3, Cauchi et al. GSO_6, Gaiotto-Witten N = 4, Nakajima G, L-functions/L-factors in Borel, Sakellaridis 2012, Sakellaridis-Wang).
Verification:
  - voice: OK
  - refs: PASS (23/23)
  - claims: PASS (83 claims, 29 computation, rebound to the recorded results and the new manuscript hash)
  - build: 17 pages, 0 missing-character warnings; all pages inspected
  - check --stage local: PASS

## 2026-10-21 — initial write and local completion

Scope: new paper written from the origin material (technical note on cone matching in Chen-Venkatesh singular examples).
Changes:
  - brief, research (tiered, corrections to origin material), sources, references.yaml (23 entries), PAPER.md, claims.yaml, editorial.md, README.
  - simulation/: exact checks (rootdata, cones, invariants, lattices, examples, figures, run_all), uv.lock; recorded run verification/model.json (26 checks passed); earlier runs archived.
Verification:
  - voice: OK (0 candidates)
  - refs: PASS (23/23)
  - claims: PASS (83 claims, 29 computation)
  - build: 17 pages, 0 missing-character warnings
  - check --stage local: PASS
Limits: textbook locators (Humphreys, Springer, Borel, Bourbaki, Kostant) not re-inspected and not entered as claims; novelty unestablished beyond a targeted search.

<!--
## YYYY-MM-DD — <pass name>
Scope: <the bounded list of changes>
Changes:
  - ...
Verification:
  - voice: <errors / warns>
  - refs: <missing / unused>
  - claims: <reconciliation result>
  - build: <pages, missing-char warnings>
  - check => PASS/FAIL
-->
