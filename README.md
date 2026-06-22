# ThreePCF triangular plot

Reproduces the *structure* of **Fig. 2 of Hahn et al. 2019**
([arXiv:1909.11107](https://arxiv.org/abs/1909.11107)) for the three-point
correlation function (3PCF).

Figure 2 of the paper shows the redshift-space halo bispectrum monopole
`B0(k1, k2, k3)` as a function of **triangle-configuration shape** in a 2×4 grid
of panels. Each panel is a 2D colour map over the side-length ratios
`(k3/k1, k2/k1)`: every triangle (with `k1 >= k2 >= k3`) lands at a point in that
plane and the colour encodes `B0`. The valid region is an inverted triangle with

- **squeezed** configurations at the top-left,
- **equilateral** at the top-right, and
- **folded** at the bottom-centre.

The top row varies the neutrino mass `Mnu`; the bottom row varies `sigma8` (with
a schematic in the lower-left panel). The three right-most columns share
`sigma8`, which is how the figure illustrates the `Mnu`–`sigma8` degeneracy.

This repository reproduces the same structure with two changes:

- the configuration-space **3PCF `zeta(r1, r2, r3)`** is plotted instead of the
  bispectrum `B0`, and
- the side ordering is **reversed**: triangles are enumerated with
  `r1 <= r2 <= r3` (the original uses `k1 >= k2 >= k3`). The largest side `r3`
  plays the role of the reference side `k1`, so the shape axes become
  **`(r1/r3, r2/r3)`** — the analogue of `(k3/k1, k2/k1)`.

The `zeta` values are **synthetic/illustrative**, from the hierarchical
(Groth–Peebles) ansatz `zeta = Q [xi1 xi2 + xi2 xi3 + xi3 xi1]` on a power-law
`xi(r)`, with the amplitude scaled by `(sigma8/sigma8_fid)^4` (since
`zeta ~ xi^2 ~ sigma8^4`). To use real data, replace `synthetic_zeta` with your
measurements.

## Files

- `three_pcf_triangular.py` — triangle enumeration
  (`build_triangle_configurations`, `r1 <= r2 <= r3`), the shape ratios
  (`TriangleConfigurations.shape_ratios`), the synthetic model
  (`synthetic_zeta`), the angle cosine (`cosine_r1_r2`), and the shape-plane
  binning (`bin_on_shape_plane`).
- `plot_3pcf_triangular.py` — builds the 2×4 grid and writes
  `zeta_triangular.png`.
- `plot_3pcf_cos_angle.py` — single triangle-shape panel coloured by the cosine
  of the angle between `r1` and `r2`; writes `cos_angle_triangular.png`.

## Usage

```bash
pip install -r requirements.txt
python plot_3pcf_triangular.py   # 2x4 grid of zeta(r1,r2,r3) shape maps
python plot_3pcf_cos_angle.py    # single panel coloured by cos(theta_12)
```

`plot_3pcf_triangular.py` writes `zeta_triangular.png`: a 2×4 grid of
triangle-shape colour maps of `zeta`, with a shared `RdBu` log colour scale
(blue = high ≈ squeezed, red = low ≈ equilateral).

`plot_3pcf_cos_angle.py` writes `cos_angle_triangular.png`: one triangle-shape
panel over the same `(r1/r3, r2/r3)` plane, coloured by

```
cos(theta_12) = (r1^2 + r2^2 - r3^2) / (2 r1 r2),
```

the cosine of the angle between sides `r1` and `r2` (opposite the longest side
`r3`). It runs from `-1` on the folded edge through `~0` at squeezed up to
`0.5` at the equilateral corner — a purely geometric quantity, independent of
the 3PCF model.

### Plugging in real data

`build_triangle_configurations(r_bins)` returns the `(r1, r2, r3)` triplets;
`cfg.shape_ratios()` gives the `(r1/r3, r2/r3)` coordinates. Evaluate your
measured `zeta` on those triplets (instead of `synthetic_zeta`) and pass it to
`bin_on_shape_plane` / the plotting code.
