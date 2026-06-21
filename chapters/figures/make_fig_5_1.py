"""
Generate fig-5-1: illustration of the solution concept for a stochastic
differential equation driven by a Poisson counting process,

    dx = f(x,t) dt + g(x,t) dN(t).

The picture shows that between arrival times t_i the path x(t) follows the
ordinary differential equation dx/dt = f(x,t) (smooth segments), while at each
arrival time t_i the path jumps by the amount g(x(t_i-), t_i).  The solution is
continuous from the left: at t_i the value x(t_i) equals the left limit
(filled dot), and the post-jump value is the right limit (open dot).

Run:  python3 make_fig_5_1.py
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- model
# A simple smooth flow between arrivals so the ODE part is clearly curved.
# Here f(x,t) = a*(c - x): x relaxes toward a "target" level c between jumps.
a, c = 1.1, 1.0


def flow(x0, t0, t1, n=200):
    """Solve dx/dt = a*(c - x) on [t0, t1] starting from x0."""
    t = np.linspace(t0, t1, n)
    x = c + (x0 - c) * np.exp(-a * (t - t0))
    return t, x


# arrival times of the Poisson counter and the jump sizes g(x(t_i-), t_i)
arrivals = [1.6, 3.1, 4.7]
jumps = [1.5, 1.1, 1.3]   # amount g added at each arrival

t_end = 6.0
x0 = 0.6   # x(0)

# ---------------------------------------------------------------- build path
segments = []          # list of (t_array, x_array)
left_limits = []       # x(t_i-)   (filled dots, = x(t_i) by left-continuity)
right_values = []      # x(t_i+)   (open dots, post-jump start)

cur_x = x0
breakpoints = [0.0] + arrivals + [t_end]
for k in range(len(arrivals) + 1):
    t0, t1 = breakpoints[k], breakpoints[k + 1]
    t, x = flow(cur_x, t0, t1)
    segments.append((t, x))
    if k < len(arrivals):
        x_minus = x[-1]                 # left limit at the arrival time
        left_limits.append((t1, x_minus))
        cur_x = x_minus + jumps[k]       # apply the jump g
        right_values.append((t1, cur_x))

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(8.5, 5.0))

for t, x in segments:
    ax.plot(t, x, color="C0", lw=2.0, solid_capstyle="round")

# vertical dashed jump lines + dots showing left-continuity
for (ti, xm), (_, xp) in zip(left_limits, right_values):
    ax.plot([ti, ti], [xm, xp], color="C0", lw=1.2, ls=(0, (4, 3)))
    ax.plot(ti, xm, "o", color="C0", ms=7, zorder=5)                 # filled: x(t_i)=x(t_i-)
    ax.plot(ti, xp, "o", mfc="white", mec="C0", mew=1.8, ms=7, zorder=5)  # open: x(t_i+)

# annotate the ODE flow on the first segment
t0seg, x0seg = segments[0]
mid = len(t0seg) // 2
ax.annotate(r"$\frac{dx}{dt}=f(x,t)$",
            xy=(t0seg[mid], x0seg[mid]),
            xytext=(0.35, 2.0),
            fontsize=13,
            arrowprops=dict(arrowstyle="->", color="0.35", lw=1.0))

# annotate the jump g at the first arrival
ti, xm = left_limits[0]
xp = right_values[0][1]
ax.annotate(r"jump $=g\left(x(t_1^-),\,t_1\right)$",
            xy=(ti, 0.5 * (xm + xp)),
            xytext=(ti + 0.25, xm - 0.9),
            fontsize=12,
            arrowprops=dict(arrowstyle="->", color="0.35", lw=1.0))
# small brace-like indicator of jump magnitude
ax.annotate("", xy=(ti - 0.12, xp), xytext=(ti - 0.12, xm),
            arrowprops=dict(arrowstyle="<->", color="0.45", lw=1.0))
ax.text(ti - 0.22, 0.5 * (xm + xp), r"$g$", ha="right", va="center",
        fontsize=12, color="0.3")

# x(0) marker
ax.plot(0, x0, "o", color="C0", ms=7, zorder=5)
ax.annotate(r"$x(0)$", xy=(0, x0), xytext=(-0.05, x0 - 0.55),
            fontsize=12, ha="left")

# arrival-time ticks on the horizontal axis
for i, ti in enumerate(arrivals, start=1):
    ax.plot([ti, ti], [-0.05, 0.05], color="black", lw=1.0,
            transform=ax.get_xaxis_transform(), clip_on=False)

ax.set_xticks(arrivals)
ax.set_xticklabels([r"$t_1$", r"$t_2$", r"$t_3$"], fontsize=13)
ax.set_yticks([])

ax.set_xlabel(r"$t$", fontsize=13)
ax.set_ylabel(r"$x(t)$", fontsize=13, rotation=0, labelpad=15)
ax.set_xlim(-0.3, t_end + 0.1)
ax.set_ylim(-0.2, 4.4)

# clean spines: keep left and bottom as axes
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig("fig-5-1_solution_concept.png", dpi=160, bbox_inches="tight")
print("wrote fig-5-1_solution_concept.png")
