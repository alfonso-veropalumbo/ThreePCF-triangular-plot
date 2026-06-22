"""Single-panel scale-averaged 3PCF over triangle shape, normalised by r3.

Run:

    python plot_3pcf_qbar.py

Produces ``qbar_triangular.png``: one panel of the scale-averaged 3PCF
Qbar(x, y) over the shape plane x = r1/r3, y = r2/r3 (the largest side r3 is the
reference, mirroring k1 in Fourier space).  Qbar is the configuration-space
analogue of Eq. (44):

    Qbar(x, y) = 1/(r_u - r_l) * int_{r_l}^{r_u} dr Qhat(x r, y r, r),

i.e. Qhat averaged over the reference (largest) side r3 = r.  The configurations
fill the inverted-triangle region: squeezed at top-left, equilateral at
top-right, folded at bottom-centre.

The colour scale adapts to the data: a diverging scale centred at zero if Qbar
changes sign, otherwise a logarithmic sequential scale.
"""

from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

from three_pcf_triangular import q_bar

# Scale-average window (reference = largest side r3).
R_L, R_U = 40.0, 90.0
# Lower bound on x = r1/r3, keeping the smallest side physical.
X_MIN = 0.05


def make_plot(outfile: str = "qbar_triangular.png") -> None:
    # Shape plane: x = r1/r3, y = r2/r3, both in (0, 1]. Valid triangles satisfy
    # x <= y <= 1 (ordering r1 <= r2 <= r3) and x + y >= 1 (triangle inequality),
    # i.e. the inverted-triangle region. x is bounded below at X_MIN so the
    # smallest side r1 = x r3 stays physical (the squeezed limit r1 -> 0 makes
    # the model two-point function diverge and would swamp the colour scale).
    x = np.linspace(X_MIN, 1.0, 300)
    y = np.linspace(0.5, 1.0, 200)
    X, Y = np.meshgrid(x, y)

    qbar = q_bar(X, Y, r_l=R_L, r_u=R_U)

    valid = (Y >= X) & (X + Y >= 1.0) & (Y <= 1.0)
    qbar = np.ma.masked_array(qbar, mask=~valid)

    # Adapt the colour scale to whatever the data actually is: diverging and
    # centred at zero if Qbar changes sign, otherwise a decade-aligned log scale.
    qmin = float(qbar.min())
    qmax = float(qbar.max())
    if qmin < 0.0 < qmax:
        vlim = max(abs(qmin), abs(qmax))
        norm = TwoSlopeNorm(vmin=-vlim, vcenter=0.0, vmax=vlim)
        cmap = plt.get_cmap("RdBu_r").copy()
    else:
        vmin = 10 ** np.floor(np.log10(qmin))
        vmax = 10 ** np.ceil(np.log10(qmax))
        norm = LogNorm(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("white")  # outside the triangle region
    cbar_label = r"$\bar{Q}(x, y)$"

    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    mesh = ax.pcolormesh(X, Y, qbar, cmap=cmap, norm=norm, shading="auto")

    # Configuration landmarks at the three corners of the inverted triangle.
    for (cx, cy, label, dx, dy, ha) in [
        (0.0, 1.0, "squeezed", 6, -2, "left"),
        (1.0, 1.0, "equilateral", -6, -2, "right"),
        (0.5, 0.5, "folded", 0, -12, "center"),
    ]:
        ax.plot(cx, cy, "k.", ms=6)
        ax.annotate(label, (cx, cy), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, fontsize=10)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.5, 1.0)
    ax.set_xlabel(r"$x = r_1 / r_3$", fontsize=14)
    ax.set_ylabel(r"$y = r_2 / r_3$", fontsize=14)
    ax.tick_params(direction="in", top=True, right=True)
    ax.set_title(
        rf"scale-averaged 3PCF, normalised by $r_3$ "
        rf"($r_3 \in [{R_L:.0f}, {R_U:.0f}]\ h^{{-1}}\mathrm{{Mpc}}$)",
        fontsize=12,
    )

    cb = fig.colorbar(mesh, ax=ax, pad=0.02)
    cb.set_label(cbar_label, fontsize=14)

    fig.tight_layout()
    fig.savefig(outfile, dpi=150)
    frac_neg = np.mean(qbar.compressed() < 0)
    print(f"Wrote {outfile} (fraction of negative shape cells = {frac_neg:.2f}).")


if __name__ == "__main__":
    make_plot()
