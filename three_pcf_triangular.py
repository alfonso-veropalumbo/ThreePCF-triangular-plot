"""Reproduce the structure of Fig. 2 of Hahn et al. 2019 (arXiv:1909.11107)
for the three-point correlation function (3PCF).

Figure 2 of the paper shows the redshift-space halo bispectrum monopole
B0(k1, k2, k3) as a function of *triangle-configuration shape* in a 2x4 grid of
panels.  Each panel is a 2D map over the side-length ratios (k3/k1, k2/k1):
every triangle (with k1 >= k2 >= k3) lands at a point in that plane and the
colour encodes B0.  The plane is the inverted-triangle region with the squeezed
configurations at the top-left, equilateral at the top-right and folded at the
bottom-centre.  The top row varies the neutrino mass Mnu and the bottom row
varies sigma8, illustrating the Mnu-sigma8 degeneracy.

Here we reproduce the same *structure* for the configuration-space 3PCF
zeta(r1, r2, r3).  Two differences from the original:

  * we plot the 3PCF zeta instead of the bispectrum B0, and
  * the side ordering is *reversed*: triangles are enumerated with
    r1 <= r2 <= r3 (the original uses k1 >= k2 >= k3).  The largest side r3
    plays the role of the reference side k1, so the shape axes become
    (r1/r3, r2/r3) -- the analogue of (k3/k1, k2/k1).

The zeta values are synthetic/illustrative: they come from the hierarchical
("Groth-Peebles") ansatz built on a power-law two-point function, scaled by
(sigma8/sigma8_fid)^4 since zeta ~ xi^2 ~ sigma8^4.  Swap `synthetic_zeta` for
measured values to use real data.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# Fiducial sigma8 of the massless-neutrino HADES cosmology, and the mapping
# between neutrino mass and the (matched) sigma8 used in Fig. 2: the massive-Mnu
# cosmologies share sigma8 with the lower-sigma8 / massless cosmologies, which is
# exactly the degeneracy the figure highlights.
SIGMA8_FID = 0.833
MNU_TO_SIGMA8 = {0.0: 0.833, 0.06: 0.822, 0.10: 0.818, 0.15: 0.807}


# --------------------------------------------------------------------------- #
# Triangle configurations
# --------------------------------------------------------------------------- #
@dataclass
class TriangleConfigurations:
    """Ordered list of triangle side bins (r1, r2, r3) with r1 <= r2 <= r3.

    Attributes are 1-D arrays of equal length ``n_config``.
    """

    r1: np.ndarray
    r2: np.ndarray
    r3: np.ndarray

    def __len__(self) -> int:
        return self.r1.size

    def shape_ratios(self) -> tuple[np.ndarray, np.ndarray]:
        """Return the shape coordinates (x, y) = (r1/r3, r2/r3).

        With r1 <= r2 <= r3 these satisfy 0 <= x <= y <= 1 and, for closed
        triangles, x + y >= 1.  The corners are:
          * squeezed    (r2 = r3 >> r1): (x, y) -> (0, 1)   top-left
          * equilateral (r1 = r2 = r3) : (x, y) -> (1, 1)   top-right
          * folded      (r3 = 2 r1 = 2 r2): (x, y) -> (0.5, 0.5) bottom-centre
        This mirrors the (k3/k1, k2/k1) plane of Fig. 2 with the reversed
        side ordering.
        """
        return self.r1 / self.r3, self.r2 / self.r3


def build_triangle_configurations(r_bins: np.ndarray) -> TriangleConfigurations:
    """Enumerate all valid triangles from a set of radial bins.

    The side ordering is *reversed* relative to Hahn et al. (2019): triangles
    are kept with r1 <= r2 <= r3.  A triangle is retained only if it closes,
    i.e. it satisfies the triangle inequality r3 <= r1 + r2 (degenerate/flat
    triangles are allowed).

    Parameters
    ----------
    r_bins:
        1-D array of bin centres (e.g. in Mpc/h), assumed sorted ascending.

    Returns
    -------
    TriangleConfigurations
        The (r1, r2, r3) triplets with r1 <= r2 <= r3.
    """
    r_bins = np.asarray(r_bins, dtype=float)

    r1_list, r2_list, r3_list = [], [], []
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
    sigma8: float = SIGMA8_FID,
    Q: float = 1.0,
    r0: float = 5.0,
    gamma: float = 1.8,
    scatter: float = 0.0,
    seed: int | None = 42,
) -> np.ndarray:
    """Illustrative 3PCF from the hierarchical (Groth-Peebles) ansatz.

        zeta = Q [ xi(r1) xi(r2) + xi(r2) xi(r3) + xi(r3) xi(r1) ]

    The amplitude is scaled by (sigma8 / sigma8_fid)^4, since zeta ~ xi^2 and
    xi ~ sigma8^2.  Lowering sigma8 (or, equivalently, raising the neutrino
    mass) suppresses zeta -- the effect that drives Fig. 2.

    An optional small log-normal scatter mimics measurement noise while staying
    reproducible through ``seed``.
    """
    xi1 = power_law_xi(cfg.r1, r0, gamma)
    xi2 = power_law_xi(cfg.r2, r0, gamma)
    xi3 = power_law_xi(cfg.r3, r0, gamma)

    zeta = Q * (xi1 * xi2 + xi2 * xi3 + xi3 * xi1)
    zeta = zeta * (sigma8 / SIGMA8_FID) ** 4

    if scatter > 0:
        rng = np.random.default_rng(seed)
        zeta = zeta * np.exp(rng.normal(0.0, scatter, size=zeta.shape))

    return zeta


def cosine_r1_r2(cfg: TriangleConfigurations) -> np.ndarray:
    """Cosine of the angle between sides r1 and r2.

    That angle sits opposite the third side r3, so by the law of cosines

        cos(theta_12) = (r1^2 + r2^2 - r3^2) / (2 r1 r2).

    With r1 <= r2 <= r3, r3 is the longest side, so theta_12 is the largest
    angle: cos ranges from -1 (folded/degenerate, r3 = r1 + r2) up to 0.5
    (equilateral); squeezed configurations give cos -> 0 (right angle).
    """
    return (cfg.r1 ** 2 + cfg.r2 ** 2 - cfg.r3 ** 2) / (2.0 * cfg.r1 * cfg.r2)


def bin_on_shape_plane(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    nx: int = 20,
    ny: int = 11,
    x_range: tuple[float, float] = (0.0, 1.0),
    y_range: tuple[float, float] = (0.5, 1.0),
    min_count: int = 2,
    log: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    """Average ``values`` over a 2D grid in the shape plane (x, y).

    With ``log=True`` (the default) a *geometric* mean is used, appropriate for
    a positive quantity spanning several decades such as the 3PCF: it gives the
    smooth gradient seen in the original figure rather than letting a single
    small-scale triangle saturate a cell.  With ``log=False`` an ordinary linear
    mean is used, suitable for signed quantities such as cos(theta_12).  Cells
    with fewer than ``min_count`` triangles are masked (drawn white).

    Returns the bin edges ``xedges``, ``yedges`` and a masked array of the mean
    value per cell, ready for ``pcolormesh``.
    """
    xedges = np.linspace(*x_range, nx + 1)
    yedges = np.linspace(*y_range, ny + 1)

    weights = np.log10(values) if log else values
    total, _, _ = np.histogram2d(x, y, bins=[xedges, yedges], weights=weights)
    count, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
        if log:
            mean = 10 ** mean
    mean = np.ma.masked_array(mean, mask=~np.isfinite(mean) | (count < min_count))

    return xedges, yedges, mean
