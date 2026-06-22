"""
Generate fig-3-2: the triangular nascent-delta function that arises as the
autocovariance of the smoothed white noise y(t) = [N(t+eps) - N(t)] / eps.

For |tau| <= eps the autocovariance (above its constant floor lambda^2) equals

    g_eps(tau) = lambda/eps - lambda |tau| / eps^2 ,

a triangle of height lambda/eps, base [-eps, eps], and area exactly lambda for
every eps > 0.  As eps -> 0 the triangle becomes taller and narrower while
keeping unit-mass-times-lambda, so in the limit it defines the Dirac delta with
mass lambda, i.e. lambda * delta(tau).

Run:  python3 make_fig_3_2.py
"""

import numpy as np
import matplotlib.pyplot as plt

lam = 1.0                       # mass of the limiting delta (lambda)
eps_values = [1.0, 0.5, 0.25]   # shrinking window widths
colors = ["C0", "C1", "C2"]

fig, ax = plt.subplots(figsize=(8.5, 5.0))

for eps, col in zip(eps_values, colors):
    tau = np.linspace(-eps, eps, 400)
    g = lam / eps - lam * np.abs(tau) / eps**2     # triangle on [-eps, eps]
    # extend with zeros outside the base so each curve reads as a function on R
    tau_full = np.concatenate([[-1.3], [-eps], tau, [eps], [1.3]])
    g_full = np.concatenate([[0.0], [0.0], g, [0.0], [0.0]])
    ax.plot(tau_full, g_full, color=col, lw=2.0,
            label=rf"$\epsilon={eps:g}$,  peak $=\lambda/\epsilon={lam/eps:g}$")

# annotate that every triangle has the same area = lambda
ax.fill_between(np.linspace(-eps_values[0], eps_values[0], 200),
                lam / eps_values[0]
                - lam * np.abs(np.linspace(-eps_values[0], eps_values[0], 200))
                / eps_values[0]**2,
                color="C0", alpha=0.12)
ax.annotate(r"area $=\lambda$ for every $\epsilon$",
            xy=(0.35, 0.55), xytext=(0.75, 2.4), fontsize=12, color="0.3",
            arrowprops=dict(arrowstyle="->", color="0.45", lw=1.0))

# arrow suggesting the eps -> 0 limit (taller, narrower spike)
ax.annotate(r"$\epsilon \to 0$", xy=(0.0, 4.4), xytext=(0.0, 3.3),
            ha="center", fontsize=13, color="0.25",
            arrowprops=dict(arrowstyle="->", color="0.25", lw=1.4))

# mark +/- eps for the widest triangle
for s in (-1, 1):
    ax.plot([s * eps_values[0], s * eps_values[0]], [-0.05, 0.05],
            color="black", lw=1.0,
            transform=ax.get_xaxis_transform(), clip_on=False)
ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
ax.set_xticklabels([r"$-\epsilon$", "", "$0$", "", r"$\epsilon$"], fontsize=12)

ax.set_xlabel(r"$\tau = t_1 - t_2$", fontsize=13)
ax.set_ylabel(r"$\dfrac{\lambda}{\epsilon} - \dfrac{\lambda|\tau|}{\epsilon^{2}}$",
              fontsize=13, rotation=0, labelpad=30)
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(0.0, 4.8)

ax.legend(fontsize=11, loc="upper right", frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig("fig-3-2_delta_triangle.png", dpi=160, bbox_inches="tight")
print("wrote fig-3-2_delta_triangle.png")
