# ThreePCF triangular plot

Reproduces the *structure* of **Fig. 2 of Hahn et al. 2019**
([arXiv:1909.11107](https://arxiv.org/abs/1909.11107)) for the three-point
correlation function (3PCF).

The original figure draws the redshift-space halo bispectrum monopole
`B0(k1, k2, k3)` as a single jagged curve against a *triangle-configuration
index*: every triangle with `k1 >= k2 >= k3` (and satisfying the triangle
inequality) is enumerated by a fixed nested loop over the three sides, and its
amplitude is plotted at the corresponding index.

This repository reproduces the same structure with two changes:

- the configuration-space **3PCF `zeta(r1, r2, r3)`** is plotted instead of the
  bispectrum `B0`, and
- the side ordering is **reversed**: triangles are enumerated with
  `r1 <= r2 <= r3` (the original uses `k1 >= k2 >= k3`).

The `zeta` values are **synthetic/illustrative**, generated from the
hierarchical (Groth–Peebles) ansatz on top of a power-law two-point function.
This is enough to reproduce the characteristic configuration dependence; to use
real data, replace `synthetic_zeta` with your measurements.

## Files

- `three_pcf_triangular.py` — triangle-configuration enumeration
  (`build_triangle_configurations`) and the synthetic 3PCF model
  (`synthetic_zeta`).
- `plot_3pcf_triangular.py` — produces `zeta_triangular.png`.

## Usage

```bash
pip install -r requirements.txt
python plot_3pcf_triangular.py
```

This writes `zeta_triangular.png`: `zeta(r1, r2, r3)` on a log scale versus the
triangle-configuration index. The configurations are ordered by looping `r1`
(outermost) → `r2` → `r3` (innermost) with `r1 <= r2 <= r3`. Light grey bands
mark blocks of constant smallest side `r1`, which exposes the nested ordering
(the analogue of the colour blocks in the original figure).

### Plugging in real data

`build_triangle_configurations(r_bins)` returns aligned `r1`, `r2`, `r3` arrays
in the plot order. Evaluate your measured `zeta` on those triplets (instead of
calling `synthetic_zeta`) and pass the result to the plotting code.
