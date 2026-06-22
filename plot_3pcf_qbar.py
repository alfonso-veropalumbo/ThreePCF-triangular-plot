"""Single-panel scale-averaged 3PCF over triangle shape, normalised by r1.

Run:

    python plot_3pcf_qbar.py

Produces ``qbar_triangular.png``: one panel of the scale-averaged 3PCF
Qbar(x2, x3) over the shape plane x2 = r2/r1, x3 = r3/r1 (the smallest side r1
is the reference -- the "upside-down" configuration-space counterpart of the
largest k1 used in Fourier space).  Qbar is the configuration-space analogue of
Eq. (44):

    Qbar(x2, x3) = 1/(r_u - r_l) * int_{r_l}^{r_u} dr Qhat(r, r x2, r x3).

The synthetic Qhat changes sign at large separation, so the colour scale is
diverging and centred (white) at zero.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from three_pcf_triangular import q_bar

# Shape grid: x2 = r2/r1 >= 1, x3 = r3/r1 >= 1. Valid triangles satisfy
# x2 <= x3 <= x2 + 1 (ordering r1 <= r2 <= r3 and the triangle inequality).
X2_MAX, X3_MAX = 4.0, 5.0
# Scale-average window (smallest side r1), chosen to straddle the zero-crossing
# regime so the sign of the 3PCF is visible.
R_L, R_U = 40.0, 90.0

CMAP = plt.get_cmap("RdBu_r").copy()
CMAP.set_bad("white")  # cells outside the triangle region stay white


def make_plot(outfile: str = "qbar_triangular.png") -> None:
    x2 = np.linspace(1.0, X2_MAX, 240)
    x3 = np.linspace(1.0, X3_MAX, 320)
    X2, X3 = np.meshgrid(x2, x3)

    qbar = q_bar(X2, X3, r_l=R_L, r_u=R_U)

    # Mask everything outside the valid band x2 <= x3 <= x2 + 1.
    valid = (X3 >= X2) & (X3 <= X2 + 1.0)
    qbar = np.ma.masked_array(qbar, mask=~valid)

    vmax = np.abs(qbar).max()
    vmin = qbar.min()
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    mesh = ax.pcolormesh(X2, X3, qbar, cmap=CMAP, norm=norm, shading="auto")

    # Configuration landmarks.
    ax.plot([1, X2_MAX], [1, X2_MAX], color="0.4", lw=1.0, ls="--")          # x3=x2 (isosceles/squeezed edge)
    ax.plot([1, X2_MAX], [2, X2_MAX + 1], color="0.4", lw=1.0, ls=":")       # x3=x2+1 (folded edge)
    ax.plot(1, 1, "k*", ms=10)                                               # equilateral
    ax.annotate("equilateral", (1, 1), textcoords="offset points",
                xytext=(8, 6), fontsize=10)
    ax.text(3.3, 3.18, "squeezed ($r_2 = r_3$)", fontsize=9, color="0.3",
            rotation=45, ha="center", va="top")
    ax.text(2.7, 3.82, "folded ($r_3 = r_1 + r_2$)", fontsize=9, color="0.3",
            rotation=45, ha="center", va="bottom")

    ax.set_xlim(1.0, X2_MAX)
    ax.set_ylim(1.0, X3_MAX)
    ax.set_xlabel(r"$x_2 = r_2 / r_1$", fontsize=14)
    ax.set_ylabel(r"$x_3 = r_3 / r_1$", fontsize=14)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_title(
        rf"scale-averaged 3PCF, normalised by $r_1$ "
        rf"($r_1 \in [{R_L:.0f}, {R_U:.0f}]\ h^{{-1}}\mathrm{{Mpc}}$)",
        fontsize=12,
    )

    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label(r"$\bar{Q}(x_2, x_3)$", fontsize=14)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    frac_neg = np.mean(qbar.compressed() < 0)
    print(f"Wrote {outfile} (fraction of negative shape cells = {frac_neg:.2f}).")


if __name__ == "__main__":
    make_plot()
