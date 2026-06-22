"""Reproduce the structure of Fig. 2 of Hahn et al. 2019 (arXiv:1909.11107)
for the three-point correlation function (3PCF).

The original figure plots the redshift-space halo bispectrum monopole
B0(k1, k2, k3) as a function of a *triangle-configuration index*: every
triangle (k1, k2, k3) with k1 >= k2 >= k3 (and satisfying the triangle
inequality) is assigned a single index by looping over the three sides in a
fixed nested order, and the amplitude is drawn as one jagged curve.

Here we reproduce the same *structure* for the configuration-space 3PCF
zeta(r1, r2, r3). Two differences from the original:

  * we plot the 3PCF zeta instead of the bispectrum B0, and
  * the side ordering is *reversed*: triangles are enumerated with
    r1 <= r2 <= r3 (the original uses k1 >= k2 >= k3).

The zeta values are synthetic/illustrative: they come from the hierarchical
("Groth-Peebles") ansatz built on a power-law two-point function, which is
enough to reproduce the characteristic configuration dependence of the curve.
Swap `synthetic_zeta` for measured values to use real data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# --------------------------------------------------------------------------- #
# Triangle configurations
# --------------------------------------------------------------------------- #
@dataclass
class TriangleConfigurations:
    """Ordered list of triangle side bins (r1, r2, r3).

    Attributes are 1-D arrays of equal length ``n_config``; element ``i`` is
    the ``i``-th triangle in the enumeration order used along the x-axis.
    """

    r1: np.ndarray
    r2: np.ndarray
    r3: np.ndarray

    def __len__(self) -> int:
        return self.r1.size


def build_triangle_configurations(r_bins: np.ndarray) -> TriangleConfigurations:
    """Enumerate all valid triangles from a set of radial bins.

    Mirrors the nested-loop ordering of Hahn et al. (2019) but with the side
    ordering *reversed*: r1 <= r2 <= r3.  The loops run with r1 outermost and
    r3 innermost, so the index increases as the smallest side r1 grows.

    A triangle is kept only if it can close, i.e. it satisfies the triangle
    inequality r3 <= r1 + r2 (degenerate/flat triangles are allowed).

    Parameters
    ----------
    r_bins:
        1-D array of bin centres (e.g. in Mpc/h), assumed sorted ascending.

    Returns
    -------
    TriangleConfigurations
        The ordered (r1, r2, r3) triplets.
    """
    r_bins = np.asarray(r_bins, dtype=float)

    r1_list, r2_list, r3_list = [], [], []
    # r1 outermost (smallest side), r3 innermost (largest side): r1 <= r2 <= r3.
    for r1 in r_bins:
        for r2 in r_bins[r_bins >= r1]:
            for r3 in r_bins[r_bins >= r2]:
                if r3 <= r1 + r2:  # triangle inequality
                    r1_list.append(r1)
                    r2_list.append(r2)
                    r3_list.append(r3)

    return TriangleConfigurations(
        r1=np.array(r1_list),
        r2=np.array(r2_list),
        r3=np.array(r3_list),
    )


# --------------------------------------------------------------------------- #
# Synthetic 3PCF model
# --------------------------------------------------------------------------- #
def power_law_xi(r: np.ndarray, r0: float = 5.0, gamma: float = 1.8) -> np.ndarray:
    """Power-law two-point correlation function xi(r) = (r/r0)^(-gamma)."""
    return (np.asarray(r, dtype=float) / r0) ** (-gamma)


def synthetic_zeta(
    cfg: TriangleConfigurations,
    Q: float = 0.7,
    r0: float = 5.0,
    gamma: float = 1.8,
    scatter: float = 0.05,
    seed: int | None = 42,
) -> np.ndarray:
    """Illustrative 3PCF from the hierarchical (Groth-Peebles) ansatz.

        zeta = Q [ xi(r1) xi(r2) + xi(r2) xi(r3) + xi(r3) xi(r1) ]

    A small log-normal scatter is added to mimic measurement noise so the
    curve looks "measured" while staying reproducible via ``seed``.
    """
    xi1 = power_law_xi(cfg.r1, r0, gamma)
    xi2 = power_law_xi(cfg.r2, r0, gamma)
    xi3 = power_law_xi(cfg.r3, r0, gamma)

    zeta = Q * (xi1 * xi2 + xi2 * xi3 + xi3 * xi1)

    if scatter > 0:
        rng = np.random.default_rng(seed)
        zeta = zeta * np.exp(rng.normal(0.0, scatter, size=zeta.shape))

    return zeta
