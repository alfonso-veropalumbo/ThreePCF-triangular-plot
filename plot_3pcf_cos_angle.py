"""Single-panel triangle-shape plot coloured by cos(theta_12).

Run:

    python plot_3pcf_cos_angle.py

Produces ``cos_angle_triangular.png``: one triangle-shape panel over the
(r1/r3, r2/r3) plane, coloured by the cosine of the angle between sides r1 and
r2.  Unlike the 3PCF map, this is a purely geometric quantity, fixed by the
triangle shape:

    cos(theta_12) = (r1^2 + r2^2 - r3^2) / (2 r1 r2).
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from three_pcf_triangular import (
    build_triangle_configurations,
    bin_on_shape_plane,
    cosine_r1_r2,
)

NX, NY = 20, 11
X_RANGE, Y_RANGE = (0.0, 1.0), (0.5, 1.0)

CMAP = plt.get_cmap("viridis").copy()
CMAP.set_bad("white")  # empty (non-triangle) cells stay white


def make_plot(outfile: str = "cos_angle_triangular.png") -> None:
    r_bins = np.arange(1, 41) * 4.0  # 4 .. 160 Mpc/h, 40 bins
    cfg = build_triangle_configurations(r_bins)
    x, y = cfg.shape_ratios()
    cos12 = cosine_r1_r2(cfg)

    # cos(theta_12) is exact per shape, so a linear mean with min_count=1 keeps
    # every populated cell (including the folded edge where cos -> -1).
    xe, ye, mean = bin_on_shape_plane(
        x, y, cos12, NX, NY, X_RANGE, Y_RANGE, min_count=1, log=False)

    fig, ax = plt.subplots(figsize=(6.8, 5.4))
    mesh = ax.pcolormesh(xe, ye, mean.T, cmap=CMAP, vmin=-1.0, vmax=0.5)

    ax.set_xlim(*X_RANGE)
    ax.set_ylim(*Y_RANGE)
    ax.set_xlabel(r"$r_1 / r_3$", fontsize=14)
    ax.set_ylabel(r"$r_2 / r_3$", fontsize=14)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_title(
        r"angle between $r_1$ and $r_2$ "
        r"(reversed ordering $r_1 \leq r_2 \leq r_3$)",
        fontsize=12,
    )

    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label(r"$\cos\theta_{12}$", fontsize=14)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Wrote {outfile} ({len(cfg)} triangle configurations).")


if __name__ == "__main__":
    make_plot()
