# Cone Matching and Integral Constraints in Relative Langlands Duality

In two singular examples of relative Langlands duality due to Chen and Venkatesh, the inequalities describing integral points of $X = \overline{T'U'\backslash G}$ coincide with the inequalities describing the weights of functions on the dual space. This paper classifies when that happens for split reductive groups with one-dimensional centre. The automorphic inequalities reduce for root types $A_2$, $B_2$, $G_2$ and $A_1^3$; the cones coincide only for $(A_2, d=1)$ and $(A_1^3, d=2)$, where $d$ is the primitive central pairing; and equal cones give equal supports only when the shifted cocharacter of the dual group is integral, a condition missing from the source's description of the coordinate ring. Of 35 labelled global forms of type $A_1^3$ with $d = 2$, five pass both integrality tests. Every finite claim is checked in exact arithmetic.

## Build

```bash
uv run build.py          # -> paper/PAPER.pdf  (vendored canonical recipe)
```

Requires `pandoc` and `xelatex` on PATH. From the workspace you can also run
`papers build cone-matching`.

## Computations

```bash
cd simulation && uv run python run_all.py   # -> simulation/output/results.json and figures, 26 checks
```

The recorded run is `verification/model.json`.

Part of [piatra-papers](https://github.com/piatra-institute). See the workspace
docs for the research and writing pipelines.
