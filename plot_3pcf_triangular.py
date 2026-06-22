"""Make the Fig. 2-style 3PCF triangle-shape plot.

Run:

    python plot_3pcf_triangular.py

Produces ``zeta_triangular.png``: a 2x4 grid of triangle-shape colour maps of
the synthetic 3PCF zeta(r1, r2, r3), reproducing the *structure* of Fig. 2 of
Hahn et al. 2019 (arXiv:1909.11107).  The top row varies the neutrino mass Mnu,
the bottom row varies sigma8 (with a schematic in the lower-left panel), and the
three right-most columns share sigma8 -- the Mnu-sigma8 degeneracy the figure
illustrates.  Side ordering is reversed (r1 <= r2 <= r3), so the shape axes are
(r1/r3, r2/r3) instead of (k3/k1, k2/k1).
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon

from three_pcf_triangular import (
    MNU_TO_SIGMA8,
    build_triangle_configurations,
    bin_on_shape_plane,
    synthetic_zeta,
)

# Shape-plane binning (mirrors the pixelation of the original figure).
NX, NY = 20, 11
X_RANGE, Y_RANGE = (0.0, 1.0), (0.5, 1.0)

CMAP = plt.get_cmap("RdBu").copy()
CMAP.set_bad("white")  # empty (non-triangle) cells stay white


def draw_schematic(ax) -> None:
    """Draw the squeezed / equilateral / folded schematic (lower-left panel)."""
    # Boundary of the valid shape region: an inverted triangle with vertices
    # (0,1) squeezed, (1,1) equilateral, (0.5,0.5) folded.
    ax.plot([0, 1], [1, 1], color="k", lw=1.5)        # top edge
    ax.plot([1, 0.5], [1, 0.5], color="k", lw=1.5)    # right edge
    ax.plot([0.5, 0], [0.5, 1], color="k", lw=1.5)    # left edge

    arrow = dict(arrowstyle="->", lw=1.2, color="k")

    # squeezed: a thin sliver triangle, arrow to the top-left corner.
    ax.add_patch(Polygon([(0.20, 0.92), (0.36, 0.90), (0.36, 0.905)],
                         closed=True, fill=False, ec="k", lw=1.0))
    ax.annotate("", xy=(0.02, 0.99), xytext=(0.20, 0.92), arrowprops=arrow)

    # equilateral: a small equilateral triangle, arrow to the top-right corner.
    ax.add_patch(Polygon([(0.70, 0.86), (0.80, 0.86), (0.75, 0.95)],
                         closed=True, fill=False, ec="k", lw=1.0))
    ax.annotate("", xy=(0.98, 0.99), xytext=(0.80, 0.92), arrowprops=arrow)

    # folded: a flat triangle, arrow down to the bottom-centre corner.
    ax.add_patch(Polygon([(0.44, 0.66), (0.56, 0.66), (0.50, 0.645)],
                         closed=True, fill=False, ec="k", lw=1.0))
    ax.annotate("", xy=(0.50, 0.52), xytext=(0.50, 0.63), arrowprops=arrow)


def make_plot(outfile: str = "zeta_triangular.png") -> None:
    # Radial bins: fine enough to populate the shape plane densely (analogue of
    # "all triangles with k <= kmax"). Units are Mpc/h; values are illustrative.
    r_bins = np.arange(1, 41) * 4.0  # 4 .. 160 Mpc/h, 40 bins
    cfg = build_triangle_configurations(r_bins)
    x, y = cfg.shape_ratios()

    # Column cosmologies. Top row: neutrino masses. Bottom row: schematic then
    # the matched sigma8 cosmologies (so columns 2-4 mirror the top row).
    top = [("0.0eV", MNU_TO_SIGMA8[0.0]),
           ("0.06eV", MNU_TO_SIGMA8[0.06]),
           ("0.1eV", MNU_TO_SIGMA8[0.10]),
           ("0.15eV", MNU_TO_SIGMA8[0.15])]
    bottom = [None,  # schematic
              (r"$\sigma_8 = 0.822$", 0.822),
              (r"$\sigma_8 = 0.818$", 0.818),
              (r"$\sigma_8 = 0.807$", 0.807)]

    # Shared colour scale from the fiducial map.
    _, _, mean_fid = bin_on_shape_plane(
        x, y, synthetic_zeta(cfg, sigma8=top[0][1]), NX, NY, X_RANGE, Y_RANGE)
    vmin = 10 ** np.floor(np.log10(mean_fid.min()))
    vmax = 10 ** np.ceil(np.log10(mean_fid.max()))
    norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(
        2, 4, figsize=(15, 6.4), sharex=True, sharey=True,
        gridspec_kw=dict(wspace=0.08, hspace=0.08, right=0.9),
    )

    mesh = None
    for row, specs in enumerate((top, bottom)):
        for col, spec in enumerate(specs):
            ax = axes[row, col]
            ax.set_xlim(*X_RANGE)
            ax.set_ylim(*Y_RANGE)
            ax.set_aspect("auto")

            if spec is None:
                draw_schematic(ax)
                continue

            label, sigma8 = spec
            zeta = synthetic_zeta(cfg, sigma8=sigma8)
            xe, ye, mean = bin_on_shape_plane(x, y, zeta, NX, NY, X_RANGE, Y_RANGE)
            mesh = ax.pcolormesh(xe, ye, mean.T, cmap=CMAP, norm=norm)

            # Panel labels, echoing the original placement: neutrino mass at
            # lower-left (top row); "0.0eV" lower-left and sigma8 lower-right
            # (bottom row).
            if row == 0:
                ax.text(0.05, 0.08, label, transform=ax.transAxes, fontsize=12)
            else:
                ax.text(0.05, 0.08, "0.0eV", transform=ax.transAxes, fontsize=12)
                ax.text(0.95, 0.08, label, transform=ax.transAxes,
                        fontsize=12, ha="right")

            ax.tick_params(direction="in", top=True, right=True)

    # Axis labels only on the outer panels.
    for ax in axes[1, :]:
        ax.set_xlabel(r"$r_1 / r_3$", fontsize=14)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$r_2 / r_3$", fontsize=14)

    # Shared colourbar.
    cax = fig.add_axes([0.915, 0.12, 0.015, 0.76])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label(r"$\zeta(r_1, r_2, r_3)$", fontsize=14)

    fig.suptitle(
        r"3PCF over triangle-configuration shape "
        r"(reversed ordering $r_1 \leq r_2 \leq r_3$)",
        fontsize=14, y=0.95,
    )
    fig.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Wrote {outfile} ({len(cfg)} triangle configurations).")


if __name__ == "__main__":
    make_plot()
