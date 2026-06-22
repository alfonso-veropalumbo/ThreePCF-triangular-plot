"""Make the Fig. 2-style 3PCF triangle-configuration plot.

Run:

    python plot_3pcf_triangular.py

Produces ``zeta_triangular.png``: the synthetic 3PCF zeta(r1, r2, r3) drawn as
a single jagged curve against the triangle-configuration index, reproducing the
structure of Fig. 2 of Hahn et al. 2019 (arXiv:1909.11107).
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from three_pcf_triangular import (
    build_triangle_configurations,
    synthetic_zeta,
)


def make_plot(outfile: str = "zeta_triangular.png") -> None:
    # Radial bins (Mpc/h): coarse enough to keep the plot readable while still
    # yielding ~10^3 triangle configurations, as in the original figure.
    r_bins = np.arange(5.0, 105.0, 5.0)  # 5 .. 100 Mpc/h, 20 bins

    cfg = build_triangle_configurations(r_bins)
    zeta = synthetic_zeta(cfg)
    index = np.arange(len(cfg))

    fig, ax = plt.subplots(figsize=(12, 4.5))

    # Light shading to mark blocks of constant smallest side r1, which exposes
    # the nested-loop ordering (the analogue of the colour blocks in the paper).
    r1_unique = np.unique(cfg.r1)
    for k, r1v in enumerate(r1_unique):
        sel = np.where(cfg.r1 == r1v)[0]
        if sel.size == 0:
            continue
        if k % 2 == 1:
            ax.axvspan(sel[0] - 0.5, sel[-1] + 0.5, color="0.92", lw=0, zorder=0)

    ax.plot(index, zeta, lw=0.9, color="C0", zorder=2)

    ax.set_yscale("log")
    ax.set_xlim(0, len(cfg))
    ax.set_xlabel("triangle configuration index", fontsize=12)
    ax.set_ylabel(r"$\zeta(r_1, r_2, r_3)$", fontsize=13)
    ax.set_title(
        r"3PCF over triangle configurations  "
        r"($r_1 \leq r_2 \leq r_3$,  $5 \leq r \leq 100\ h^{-1}\mathrm{Mpc}$)",
        fontsize=12,
    )

    # Annotate the ordering, echoing the original figure's style.
    ax.text(
        0.015,
        0.06,
        "configurations ordered by looping\n"
        r"$r_1$ (outer) $\to r_2 \to r_3$ (inner),  $r_1 \leq r_2 \leq r_3$",
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85),
    )

    ax.tick_params(direction="in", top=True, right=True)
    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    print(f"Wrote {outfile} with {len(cfg)} triangle configurations.")


if __name__ == "__main__":
    make_plot()
