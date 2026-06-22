"""
Generate fig-3-1: a realization (sample path) of a Poisson counting process N(t).

N(t) counts the number of arrivals in [0, t].  It is a right-continuous staircase
that starts at N(0) = 0 and jumps up by one unit at each random arrival time
t_1, t_2, ... .  At an arrival time t_i the process takes its post-jump value
(filled dot); the pre-jump left limit N(t_i^-) is shown as an open dot.  The gaps
between successive arrival times (the interarrival times) are independent and
exponentially distributed with mean 1/lambda.

Run:  python3 make_fig_3_1.py
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- arrival times
# Fixed (illustrative) arrival times so the labelled ticks stay readable.
arrivals = [0.7, 1.5, 2.3, 3.6, 4.4, 5.5]
t_end = 6.2
n = len(arrivals)

breakpoints = [0.0] + arrivals + [t_end]

# ---------------------------------------------------------------- plot
fig, ax = plt.subplots(figsize=(8.5, 5.0))

# horizontal segments: N(t) = k on [t_k, t_{k+1})
for k in range(n + 1):
    t0, t1 = breakpoints[k], breakpoints[k + 1]
    ax.plot([t0, t1], [k, k], color="C0", lw=2.0, solid_capstyle="round")

# vertical dashed jump lines + dots showing right-continuity
for k, ti in enumerate(arrivals, start=1):
    ax.plot([ti, ti], [k - 1, k], color="C0", lw=1.2, ls=(0, (4, 3)))
    ax.plot(ti, k, "o", color="C0", ms=7, zorder=5)                       # value N(t_i)
    ax.plot(ti, k - 1, "o", mfc="white", mec="C0", mew=1.8, ms=7, zorder=5)  # left limit

# N(0) = 0 marker
ax.plot(0, 0, "o", color="C0", ms=7, zorder=5)

# annotate a unit jump
ax.annotate(r"$+1$", xy=(arrivals[2], 2.5),
            xytext=(arrivals[2] + 0.25, 2.4), fontsize=12, color="0.3",
            arrowprops=dict(arrowstyle="->", color="0.45", lw=1.0))

# indicate the first interarrival time t_2 - t_1 with a double arrow
y_br = -0.55
ax.annotate("", xy=(arrivals[0], y_br), xytext=(arrivals[1], y_br),
            arrowprops=dict(arrowstyle="<->", color="0.45", lw=1.0),
            annotation_clip=False)
ax.text(0.5 * (arrivals[0] + arrivals[1]), y_br - 0.28,
        r"$t_2 - t_1$", ha="center", va="top", fontsize=11, color="0.3")

# arrival-time ticks on the horizontal axis
ax.set_xticks(arrivals)
ax.set_xticklabels([rf"$t_{{{i}}}$" for i in range(1, n + 1)], fontsize=13)
ax.set_yticks(range(n + 1))

ax.set_xlabel(r"$t$", fontsize=13)
ax.set_ylabel(r"$N(t)$", fontsize=13, rotation=0, labelpad=18)
ax.set_xlim(-0.3, t_end + 0.1)
ax.set_ylim(-1.0, n + 0.5)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout()
fig.savefig("fig-3-1_poisson_sample_path.png", dpi=160, bbox_inches="tight")
print("wrote fig-3-1_poisson_sample_path.png")
