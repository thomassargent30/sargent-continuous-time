"""
Reconstruct Figures 5 and 6 of Chapter 16 (Hansen-Sargent, "Two Difficulties
in Interpreting Vector Autoregressions"), from the parameters of the numerical
example in Table 1.

Continuous time univariate process:

    (D^3 + .6 D^2 + .4 D + .2) z(t) = w(t),      psi(D) = 1.

Roots of theta(s) = s^3 + .6 s^2 + .4 s + .2:

    lambda_1 = -.5424,   lambda_{2,3} = -.0288 +/- .6066 i.

Continuous-time MA kernel        p(tau) = sum_j delta_j e^{lambda_j tau}.
Discrete-time MA kernel          C_k    = MA coefficients of c(L)/d(L).
"Aggregation" kernel             f(tau) = sum_j V_j p(tau - j),  V(L)=d(L)/c(L).

The discrete-time AR/MA polynomials produced by the spectral factorization
(reported in Table 1) are

    d(L) = 1 - 2.1779 L + 1.8722 L^2 - .5485 L^3,
    c(L) = 1 +  .4800 L +  .0192 L^2.

Run:  python3 make_fig_16.py
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- continuous time
# roots of theta(s) (given in Table 1)
lam = np.array([-0.5424,
                -0.0288 + 0.6066j,
                -0.0288 - 0.6066j])


def partial_fraction_residues(roots):
    """delta_j = 1 / prod_{k!=j}(lambda_j - lambda_k)  for P(s)=1/theta(s)."""
    d = np.empty(len(roots), dtype=complex)
    for j in range(len(roots)):
        prod = 1.0 + 0j
        for k in range(len(roots)):
            if k != j:
                prod *= (roots[j] - roots[k])
        d[j] = 1.0 / prod
    return d


delta = partial_fraction_residues(lam)


def p_of(tau):
    """Continuous-time moving-average kernel p(tau) (real)."""
    tau = np.asarray(tau, dtype=float)
    out = np.zeros_like(tau, dtype=complex)
    for dj, lj in zip(delta, lam):
        out += dj * np.exp(lj * tau)
    return out.real


# ---------------------------------------------------------------- discrete time
# AR and MA polynomials from the spectral factorization (Table 1)
d_poly = np.array([1.0, -2.1779, 1.8722, -0.5485])   # d_0..d_3
c_poly = np.array([1.0,  0.4800, 0.0192])            # c_0..c_2


def ma_coeffs(num, den, n):
    """MA coefficients C_0..C_{n} of the rational filter num(L)/den(L)."""
    C = np.zeros(n + 1)
    for k in range(n + 1):
        val = num[k] if k < len(num) else 0.0
        for i in range(1, min(k, len(den) - 1) + 1):
            val -= den[i] * C[k - i]
        C[k] = val            # den[0] == 1
    return C


KMAX = 20
C = ma_coeffs(c_poly, d_poly, KMAX)                  # discrete-time kernel C_k

# V(L) = d(L)/c(L) = C(L)^{-1};  needed for f(tau) = sum_j V_j p(tau-j)
V = ma_coeffs(d_poly, c_poly, 40)


def f_of(tau):
    """Aggregation kernel f(tau) = sum_{j>=0} V_j p(tau-j), p(s)=0 for s<0."""
    tau = np.asarray(tau, dtype=float)
    out = np.zeros_like(tau, dtype=float)
    for j, Vj in enumerate(V):
        s = tau - j
        out += Vj * np.where(s >= 0, p_of(np.maximum(s, 0.0)), 0.0)
    return out


# ---------------------------------------------------------------- sanity check
print("delta_j (residues of 1/theta):")
for j, dj in enumerate(delta, 1):
    print(f"  j={j}: {dj.real:+.4f} {dj.imag:+.4f}i")
print("\n  tau     p(tau)        f(tau)        C_k")
ints = list(range(0, 21))
for k in ints:
    cval = f"{C[k]:+.6f}"
    print(f"  {k:2d}   {p_of(k):+.6f}   {f_of(k):+.6f}   {cval}")

# ---------------------------------------------------------------- Figure 5
tau_fine = np.linspace(0, KMAX, 1000)
p_fine = p_of(tau_fine)
k = np.arange(0, KMAX + 1)
p_int = p_of(k)

fig, ax = plt.subplots(figsize=(8.5, 5.0))

ax.plot(tau_fine, p_fine, color="C0", lw=2.0, zorder=2,
        label=r"$p(\tau)$  (continuous time)")
ax.plot(k, p_int, "o", color="C0", ms=5, zorder=3)

ax.plot(k, C, "--", color="C3", lw=1.6, zorder=2,
        label=r"$C_k$  (discrete time)")
ax.plot(k, C, "s", color="C3", ms=4.5, zorder=3)

ax.axhline(0, color="0.6", lw=0.8)
ax.set_xlim(-0.3, KMAX + 0.3)
ax.set_xticks(range(0, KMAX + 1, 2))
ax.set_xlabel(r"$j$", fontsize=13)
ax.set_ylabel(r"kernel", fontsize=12)
ax.legend(frameon=False, fontsize=12, loc="upper right")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig("fig-16-5_ma_kernels.png", dpi=160, bbox_inches="tight")
print("\nwrote fig-16-5_ma_kernels.png")

# ---------------------------------------------------------------- Figure 6
tau6 = np.linspace(0, 2, 600)
f6 = f_of(tau6)
tau6_pts = np.round(np.arange(0, 2.0001, 0.1), 3)
f6_pts = f_of(tau6_pts)

fig, ax = plt.subplots(figsize=(8.0, 4.6))
ax.plot(tau6, f6, color="C0", lw=2.0, zorder=2)
ax.plot(tau6_pts, f6_pts, "o", color="C0", ms=4.5, zorder=3)
ax.axvline(1.0, color="0.6", lw=0.9, ls=(0, (4, 3)))
ax.text(1.02, 0.02, r"$\tau = 1$", fontsize=11, color="0.4")
ax.set_xlim(0, 2)
ax.set_ylim(0, 0.6)
ax.set_xlabel(r"$\tau$", fontsize=13)
ax.set_ylabel(r"$f(\tau)$", fontsize=13)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig("fig-16-6_f_function.png", dpi=160, bbox_inches="tight")
print("wrote fig-16-6_f_function.png")
