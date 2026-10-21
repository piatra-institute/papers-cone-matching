# Research

Findings, tiered by source proximity. See the workspace docs (`papers docs`): research-pipeline.md §2.
T1 primary · T2 authoritative secondary · T3 reference · T4 general web (leads only).
A claim that reaches the paper rests on a T1 or T2 source.

## Findings

### Corrections to the origin material

- The origin material cited Chen–Venkatesh by the journal DOI (Selecta 32, article 97, 2026). The publisher PDF is behind a bot challenge; arXiv v1 (2405.18212) was read in full and the journal HTML was read through the publisher page on 2026-10-21 for Proposition 5.1, §5.1, (57), (66), (70), §6.5 and Remark 7.4. Numbering of the cited items agrees between versions; the journal version writes $\kappa$ where arXiv v1 writes $a$, and adds in §6.5 that $\check\rho'$ descends to a genuine cocharacter for the dual Garrett group.
- The origin material presented the convex-redundancy property as true for $B_2$ and $G_2$ only as a diagnostic. Re-derivation shows more: Theorem 3 classifies the property completely ($A_2$, $B_2$, $G_2$, $A_1^3$ and no others), so the literal reading of CV Remark 7.4 is false, while the full cone matching (Corollary 4) agrees with it.
- The origin material described the integrality condition only through two examples and an $\mathbb{F}_2^4$ count. The paper adds the general criterion (Theorem 6: equal cones give equal supports iff $\check\rho'$ of the dual group is a cocharacter), the count of admissible codes (35, with 12 one-sided and 18 failing both), the identification of the three self-dual codes with $\mathrm{SO}_4\times\mathrm{GL}_2$, and the degree-by-degree discrepancy for $\mathrm{GL}_2\times\mathrm{PGL}_2\times\mathrm{PGL}_2$ (first at $z = 3$: 12 versus 10 labels).
- The origin material called Garrett's groups "$\mathrm{SL}_2^3\times\mathbb{G}_m$ or $\mathrm{PGL}_2^3\times\mathbb{G}_m$" in places; these are the isogeny classes named in CV Theorem 7.1. The groups actually used are (51) and (53), with codes $\langle u\rangle$ and $u^\perp$.
- A claimed "physical parallel" (line operators) and Nakajima's composition theorem were used in the origin material as support for a philosophical conclusion. The paper uses Aharony–Seiberg–Tachikawa only at abstract level as an analogy, says that it is not an identification, and does not cite the composition theorem (not inspected).
- All computations in the origin packages were reimplemented independently; no code from the archives was executed.

### Primary mathematics (T1)

- [T1] Chen, Venkatesh 2026 (arXiv v1 read; journal checked) — §1.1 (coincidences specific to the groups), §1.4 (naive basic function), §5.1 (rational $\check\theta$, finite covers), Proposition 5.1 (32), proof with monomial $x^{\langle\check\rho',\lambda\rangle}$, Lemma 5.2, Proposition 5.5 (45)–(46), §6.1 ($\check\rho'$ for $\mathrm{PGL}_3\times\mathbb{G}_m$), §6.4–6.5 (Garrett groups (51), (53)), Theorem 7.1, (56)–(58), Lemma 7.2, (62)–(66), Proposition 7.3, (67)–(70), Remark 7.4.
- [T1] Chen 2024 (arXiv 2405.18231v2 read in part) — §1.2 (integrality on $X$ ↔ weights on $\check X$), Definition 2.3 (toric dual varieties by dual cones), Theorem 4.1, §6 and Theorem 6.10 (Deligne–Mumford stacks; expectation for central isogenies of reductive groups).
- [T1] Ben-Zvi, Sakellaridis, Venkatesh 2024 — abstract and framing only.
- [T1] Cauchi, Chen, Gutierrez Terradillos 2026 — abstract: GSO6 singular example, stacky spectral side.
- [T1] Nakajima 2024 — abstract: S-dual under cotangent-type assumption, from Coulomb branches.
- [T1] Aharony, Seiberg, Tachikawa 2013 — abstract: line operators distinguish theories with the same local correlators.
- [T1] Feng, Wang 2025; Knop, Schalke 2017; Sakellaridis 2012; Sakellaridis, Wang 2022; Bouthier, Ngô, Sakellaridis 2016; Braverman, Kazhdan 1999; Sakellaridis, Venkatesh 2017 — abstracts read (Crossref/OpenAlex/Semantic Scholar).
- [T1] Ginzburg 1991; Garrett 1987 — identities checked in Crossref; roles as described in CV §1.1.
- [T1] Kostant 1963 — standard result that the nilpotent cone is the complete intersection of the basic invariants; used only to name the Hilbert series of the $\mathfrak{sl}_3$ nilcone. Passage not re-inspected.
- [T2] Beuzart-Plessis 2025 (introduction to the relative Langlands program) — abstract.
- [T2] Humphreys 1972 §10.2; Springer 1998; Borel 1979; Bourbaki 2002 Plates I–IX — standard facts on root systems, root data and dual groups. Not re-inspected; every fact used is also proved in the text or checked by the program.

### Results of this paper (proved and computed)

- Lemma 1 ($d = |Z^0\cap G_{\mathrm{der}}|$, $\theta_G = d\check\theta_{\check G}$): proof in text.
- Lemma 2, Theorem 3, Corollary 4: proofs in text; exact certificates for all 470 types of rank at most 8 with $\kappa > 1$ (recorded run verification/model.json).
- Proposition 5: proof in text; 799 + 144 direct invariant computations; Hilbert series of the two published $X$ reproduced through degree 12.
- Theorem 6, Proposition 7: proofs in text; enumeration of 67 subspaces; support comparison for $\mathrm{GL}_2\times\mathrm{PGL}_2\times\mathrm{PGL}_2$ through $z = 8$ in both directions.

### Open (stated in the paper)

- Orbit condition of CV §7.3 for $\mathrm{SO}_4\times\mathrm{GL}_2$; identities for the 12 one-sided forms; modifications for $B_2$, $G_2$.
- No priority search can certify novelty; a targeted search (arXiv listing of relative Langlands papers through 2026-09, two web searches) found no erratum to CV Proposition 5.1 and no classification of the convex reduction.
