"""
Clean redraws of the eight hand-drawn figures used in the streamlined
chapter 18 (Marcet, "Temporal Aggregation of Economic Time Series").

The kernel/projection figures are *computed*, not sketched: given a continuous
moving-average kernel a(u), we form the L^2 projection h = eta(a | A) of a onto
the closed span of its integer shifts {a(u-k), k>=1} (Proposition 1), the
residual c = a - h (the weighting kernel of the sampled innovation), and the
discrete MAR coefficients A_k via Proposition 3,

    A_k = [\int a(u+k) c(u) du] / [\int c(u)^2 du].

Run:  python3 make_fig_18.py
"""

import numpy as np
import matplotlib.pyplot as plt

BLUE, RED, GREY = "C0", "C3", "0.55"

# ---------------------------------------------------------------- grid + tools
du = 0.005
u = np.arange(0.0, 14.0, du)


def shift(vals, k):
    """Return a(u-k) on the grid (zero where u-k < 0)."""
    s = int(round(k / du))
    out = np.zeros_like(vals)
    if s < len(vals):
        out[s:] = vals[: len(vals) - s]
    return out


def project_onto_shifts(a, K=10):
    """h = eta(a | span{a(u-k), k=1..K}) in L^2; returns h, residual c, coeffs."""
    B = np.column_stack([shift(a, k) for k in range(1, K + 1)])
    G = (B.T @ B) * du
    rhs = (B.T @ a) * du
    lam = np.linalg.solve(G + 1e-10 * np.eye(K), rhs)
    h = B @ lam
    return h, a - h, lam


def discrete_mar(a, c, kmax):
    """A_k = [\int a(u+k)c du] / [\int c^2 du], k = 0..kmax."""
    denom = np.sum(c * c) * du
    A = []
    for k in range(kmax + 1):
        num = np.sum(shift(a, -k) * c) * du if False else np.sum(a[int(round(k/du)):] * c[:len(a) - int(round(k/du))]) * du
        A.append(num / denom)
    return np.array(A)


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ================================================================ Figure 1: c
# A kernel close to a continuous AR(1), so that c is small on [1, infinity).
a1 = np.exp(-0.7 * u) + 0.28 * np.exp(-2.2 * u)
h1, c1, _ = project_onto_shifts(a1)
ratio = np.sum(c1[u >= 1] ** 2) / np.sum(c1[u < 1] ** 2)
print(f"fig-ta-c:  integral c^2 on [1,inf) / [0,1) = {ratio:.3f}")

fig, ax = plt.subplots(figsize=(8.5, 4.6))
m = u <= 4.6
ax.plot(u[m], c1[m], color=BLUE, lw=2.0)
ax.axhline(0, color="0.6", lw=0.8)
for x in range(1, 5):
    ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
ax.annotate(r"$c = a$ on $[0,1)$", xy=(0.45, c1[u <= 0.45][-1]),
            xytext=(1.4, 0.85 * c1.max()), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0))
ax.annotate("small on $[1,\\infty)$", xy=(2.3, c1[u <= 2.3][-1]),
            xytext=(2.7, 0.45 * c1.max()), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0))
ax.set_xlim(-0.15, 4.7)
ax.set_xticks(range(0, 5))
ax.set_xlabel("$u$", fontsize=13)
ax.set_ylabel("$c_{ij}(u)$", fontsize=13)
clean(ax)
fig.tight_layout()
fig.savefig("fig-18-1.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-1.png")

# ================================================================ Figure 2: bias
# Continuous MAR a vs discrete MAR coefficients A_k (step levels).
a2 = np.exp(-0.55 * u) + 0.45 * np.exp(-0.18 * u)
a2 = a2 / a2[0]                       # normalize a(0) = 1
h2, c2, _ = project_onto_shifts(a2, K=12)
A = discrete_mar(a2, c2, kmax=5)
V = np.sum(c2 * c2) * du
print(f"fig-ta-bias: A_0={A[0]:.3f}, a(0)=1.0 ; sum A_k^2 * V = {np.sum(A**2)*V:.3f}, "
      f"int a^2 = {np.sum(a2**2)*du:.3f}")

fig, ax = plt.subplots(figsize=(8.5, 4.8))
m = u <= 6
ax.plot(u[m], a2[m], color=BLUE, lw=2.0, label=r"continuous MAR $a(u)$")
for k in range(6):
    ax.hlines(A[k], k, k + 1, color=RED, lw=2.4,
              label=(r"discrete MAR coefficients $A_k$" if k == 0 else None))
    if k > 0:
        ax.plot([k, k], [A[k - 1], A[k]], color=RED, lw=1.0, ls=(0, (2, 2)))
ax.annotate(r"discrete levels lie above $a(u)$", xy=(0.6, A[0]),
            xytext=(1.5, A[0] + 0.07), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="0.4", lw=1.0))
ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(-0.15, 6.1)
ax.set_xticks(range(0, 7))
ax.set_xlabel(r"$u$ (continuous) / $k$ (discrete)", fontsize=12)
ax.set_ylabel("kernel", fontsize=12)
ax.legend(frameon=False, fontsize=11, loc="upper right")
clean(ax)
fig.tight_layout()
fig.savefig("fig-18-2.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-2.png")

# ================================================================ Figure 3: osc
# |a(u)| = e^{-lambda u}, sign flips at u = k + eta within each unit interval.
lam, eta = 0.55, 0.45
fig, ax = plt.subplots(figsize=(8.5, 4.6))
uu = np.arange(0.0, 3.0, du)
env = np.exp(-lam * uu)
sign = np.where((uu % 1.0) < eta, 1.0, -1.0)
a3 = sign * env
# plot each continuous piece separately so sign jumps are not bridged
edges = np.where(np.diff(sign) != 0)[0]
bounds = np.concatenate(([0], edges + 1, [len(uu)]))
for i in range(len(bounds) - 1):
    sl = slice(bounds[i], bounds[i + 1])
    ax.plot(uu[sl], a3[sl], color=BLUE, lw=2.0)
ax.plot(uu, env, color=GREY, lw=1.0, ls=(0, (4, 3)))
ax.plot(uu, -env, color=GREY, lw=1.0, ls=(0, (4, 3)))
ax.text(2.45, env[uu <= 2.45][-1] + 0.03, r"$+e^{-\lambda u}$", color=GREY, fontsize=11)
ax.text(2.45, -env[uu <= 2.45][-1] - 0.07, r"$-e^{-\lambda u}$", color=GREY, fontsize=11)
ax.axhline(0, color="0.6", lw=0.8)
for x in [eta, 1, 1 + eta, 2, 2 + eta]:
    ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
ax.set_xticks([eta, 1, 1 + eta, 2, 2 + eta])
ax.set_xticklabels([r"$\eta$", "1", r"$1+\eta$", "2", r"$2+\eta$"], fontsize=12)
ax.set_xlim(-0.1, 3.0)
ax.set_xlabel("$u$", fontsize=13)
ax.set_ylabel("$a(u)$", fontsize=13)
clean(ax)
fig.tight_layout()
fig.savefig("fig-18-3.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-3.png")

# ================================================================ Figure 4: construct
# a whose shape on [1,2) is NOT a scaled copy of its shape on [0,1);
# the one-lag guess h = mu*a(u-1) therefore misses it -> large residual.
uu = np.arange(0.0, 2.8, du)
a4 = np.zeros_like(uu)
m0 = uu < 1
a4[m0] = 0.45 + 0.55 * np.sin(np.pi * uu[m0])          # broad hump on [0,1)
m1 = (uu >= 1) & (uu < 2)
a4[m1] = 0.45 - 0.25 * (uu[m1] - 1)                    # downward ramp on [1,2), a(1)=0.45
m2 = uu >= 2
a4[m2] = 0.20 * np.exp(-(uu[m2] - 2))                  # continuous tail (a(2)=0.20)
mu = 0.55                                              # one-lag projection coefficient
h4 = np.full_like(uu, np.nan)
mm = uu >= 1
h4[mm] = mu * (0.45 + 0.55 * np.sin(np.pi * (uu[mm] - 1)))
h4[uu >= 2] = np.nan

fig, ax = plt.subplots(figsize=(8.5, 4.6))
ax.plot(uu, a4, color=BLUE, lw=2.0, label=r"$a(u)$")
ax.plot(uu, h4, color=RED, lw=2.0, ls=(0, (5, 3)),
        label=r"guess $h(u)=\mu\,a(u-1)$ on $[1,2)$")
# shade the residual on [1,2)
reg = (uu >= 1) & (uu <= 2)
ax.fill_between(uu[reg], a4[reg], np.nan_to_num(h4[reg]), color=RED, alpha=0.12)
ax.text(1.5, 0.13, r"residual $\int_1^2 c^2$", color="0.35", fontsize=11, ha="center")
ax.axhline(0, color="0.6", lw=0.8)
for x in [1, 2]:
    ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
ax.set_xticks([0, 1, 2])
ax.set_xlim(-0.1, 2.8)
ax.set_ylim(0, 1.15)
ax.set_xlabel("$u$", fontsize=13)
ax.set_ylabel("$a,\\ h$", fontsize=13)
ax.legend(frameon=False, fontsize=11, loc="upper right")
clean(ax)
fig.tight_layout()
fig.savefig("fig-18-4.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-4.png")

# ================================================================ Figure 5: avg
# a(u) = u e^{1-u}: a(0)=0, peak at u=1.  h = mu a(u-1) (one-lag projection).
a5 = u * np.exp(1.0 - u)
num = np.sum(a5[u >= 1] * shift(a5, 1)[u >= 1]) * du
den = np.sum(shift(a5, 1)[u >= 1] ** 2) * du
mu5 = num / den
h5 = np.full_like(u, np.nan)
h5[u >= 1] = mu5 * shift(a5, 1)[u >= 1]
print(f"fig-ta-avg: one-lag mu = {mu5:.3f}")

fig, ax = plt.subplots(figsize=(8.5, 4.6))
m = u <= 4
ax.plot(u[m], a5[m], color=BLUE, lw=2.0, label=r"$a(u)$,  $a(0)=0$")
ax.plot(u[m], h5[m], color=RED, lw=2.0, ls=(0, (5, 3)),
        label=r"guess $h(u)=\mu\,a(u-1)$")
reg = (u >= 1) & (u <= 2)
ax.fill_between(u[reg], a5[reg], np.nan_to_num(h5[reg]), color=RED, alpha=0.12)
ax.text(1.5, 0.18, r"large residual on $[1,2)$", color="0.35", fontsize=11, ha="center")
ax.axhline(0, color="0.6", lw=0.8)
for x in [1, 2]:
    ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
ax.set_xticks([0, 1, 2, 3, 4])
ax.set_xlim(-0.1, 4.0)
ax.set_ylim(0, 1.1)
ax.set_xlabel("$u$", fontsize=13)
ax.set_ylabel("$a,\\ h$", fontsize=13)
ax.legend(frameon=False, fontsize=11, loc="upper right")
clean(ax)
fig.tight_layout()
fig.savefig("fig-18-5.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-5.png")

# ================================================================ Figure 7: Granger 2x2
def hump(uu, c, w, amp, base=0.0):
    return base + amp * np.exp(-((uu - c) / w) ** 2)

uu = np.arange(0.0, 2.8, du)
fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.2))

# a11: hump on [0,1), zero after
a11 = np.where(uu < 1, hump(uu, 0.28, 0.34, 0.75, 0.18), 0.0)
# a12: hump with a small secondary wiggle on [0,1), zero after
a12 = np.where(uu < 1, hump(uu, 0.22, 0.26, 0.55, 0.15) + hump(uu, 0.55, 0.10, 0.12), 0.0)
# a21: hump on [0,1), then a smooth declining tail (nonzero on [1, infinity))
a21 = hump(uu, 0.24, 0.30, 0.45, 0.0) + 0.30 * np.exp(-0.5 * uu) + 0.05
# a22: monotone decay throughout
a22 = 0.78 * np.exp(-0.5 * uu) + 0.05

for ax, a, lab in zip(axes.ravel(), [a11, a12, a21, a22],
                      [r"$a_{11}$", r"$a_{12}$", r"$a_{21}$", r"$a_{22}$"]):
    if lab in (r"$a_{11}$", r"$a_{12}$"):
        # draw compact support: solid on [0,1), then flat zero
        ax.plot(uu[uu < 1], a[uu < 1], color=BLUE, lw=2.0)
        ax.plot(uu[uu >= 1], a[uu >= 1], color=BLUE, lw=2.0)
    else:
        ax.plot(uu, a, color=BLUE, lw=2.0)
    ax.axhline(0, color="0.6", lw=0.8)
    for x in [1, 2]:
        ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
    ax.set_xticks([1, 2])
    ax.set_xlim(-0.05, 2.8)
    ax.set_ylim(-0.05, 1.0)
    ax.set_ylabel(lab, fontsize=13, rotation=0, labelpad=14)
    ax.set_xlabel("$u$", fontsize=12)
    clean(ax)
fig.tight_layout()
fig.savefig("fig-18-7.png", dpi=160, bbox_inches="tight")
print("wrote fig-18-7.png")

# ================================================================ Figures 11 & 12
# Sampled vs averaged: approximating the first kernel by combinations of its
# own integer shifts.  a1(0) != 0 -> good fit;  tilde a1(0) = 0 -> poor fit.
a_s = np.exp(-0.7 * u) + 0.5 * np.exp(-0.25 * u)        # sampled kernel, a(0) != 0
a_s = a_s / a_s[0]
# unit average -> tilde a (continuous, vanishes at 0)
kernel_len = int(round(1.0 / du))
csum = np.concatenate(([0.0], np.cumsum(a_s) * du))
a_avg = (csum[1:] - np.concatenate((np.zeros(kernel_len), csum[1:-kernel_len]))) if False else None
# tilde a(u) = \int_0^1 a(u-s) ds = F(u) - F(u-1)
F = np.cumsum(a_s) * du
a_avg = F - shift(F, 1)

gs, _, _ = project_onto_shifts(a_s, K=10)
ga, _, _ = project_onto_shifts(a_avg, K=10)
err_s = np.sum((a_s[u >= 1] - gs[u >= 1]) ** 2) * du
err_a = np.sum((a_avg[u >= 1] - ga[u >= 1]) ** 2) * du
print(f"fig-ta-11/12: fit error sampled={err_s:.4f}, averaged={err_a:.4f}")

for fname, a_t, g_t, ylab, title in [
    ("fig-18-11.png", a_s, gs, r"$a_1,\ g$", "sampled"),
    ("fig-18-12.png", a_avg, ga, r"$\tilde a_1,\ \hat g$", "averaged"),
]:
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    m = u <= 4.5
    ax.plot(u[m], a_t[m], color=BLUE, lw=1.7,
            label=(r"$a_1$" if title == "sampled" else r"$\tilde a_1$"))
    gplot = np.where(u >= 1, g_t, np.nan)
    ax.plot(u[m], gplot[m], color=RED, lw=2.6,
            label=(r"$g$" if title == "sampled" else r"$\hat g$"))
    ax.axhline(0, color="0.6", lw=0.8)
    for x in [1, 2, 3]:
        ax.axvline(x, color="0.8", lw=0.8, ls=(0, (3, 3)))
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlim(-0.1, 4.5)
    ax.set_xlabel("$u$", fontsize=13)
    ax.set_ylabel(ylab, fontsize=13)
    ax.legend(frameon=False, fontsize=12, loc="upper right")
    clean(ax)
    fig.tight_layout()
    fig.savefig(fname, dpi=160, bbox_inches="tight")
    print(f"wrote {fname}")
