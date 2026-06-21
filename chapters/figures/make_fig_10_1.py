"""
Generate fig-10-1: the band-pass filter ("window") B_{cd}(w) used in the
Cramer-representation chapter,

    B_{cd}(w) = 1  for  w in [c, d]  or  w in [-d, -c],
              = 0  otherwise,                                with d > c > 0.

The filter is defined on the whole frequency axis w in (-infinity, +infinity)
and is symmetric about w = 0: it passes the two symmetric frequency bands
[c, d] and [-d, -c] and stops everything else.

Run:  python3 make_fig_10_1.py
"""

import numpy as np
import matplotlib.pyplot as plt

# band edges  (0 < c < d)
c, d = 1.0, 2.0
w_max = 3.2

fig, ax = plt.subplots(figsize=(8.5, 3.6))


def draw_band(left, right):
    """Draw one unit-height pass band on [left, right] with left/right risers."""
    # the flat top at height 1
    ax.plot([left, right], [1, 1], color="C0", lw=2.2, solid_capstyle="butt")
    # vertical risers up to the top
    ax.plot([left, left], [0, 1], color="C0", lw=2.2)
    ax.plot([right, right], [0, 1], color="C0", lw=2.2)


# zero level across the whole axis (the stop band)
ax.plot([-w_max, w_max], [0, 0], color="C0", lw=2.2, solid_capstyle="butt", zorder=1)

# the two symmetric pass bands
draw_band(c, d)
draw_band(-d, -c)

# guide lines from the band edges up to height 1
for x in (c, d, -c, -d):
    ax.plot([x, x], [0, 1], color="0.7", lw=0.8, ls=(0, (3, 3)), zorder=0)
ax.plot([-w_max, w_max], [1, 1], color="0.85", lw=0.8, ls=(0, (3, 3)), zorder=0)

# axis frequency-tick labels
ax.set_xticks([-d, -c, c, d])
ax.set_xticklabels([r"$-d$", r"$-c$", r"$c$", r"$d$"], fontsize=13)
ax.set_yticks([1])
ax.set_yticklabels(["1"], fontsize=12)

# arrowheads on the frequency axis to stress w in (-inf, +inf)
ax.annotate("", xy=(w_max, 0), xytext=(-w_max, 0),
            arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
ax.text(w_max, -0.12, r"$w$", fontsize=13, ha="center", va="top")

# vertical axis with arrowhead
ax.annotate("", xy=(0, 1.35), xytext=(0, 0),
            arrowprops=dict(arrowstyle="->", color="black", lw=1.0))
ax.text(0.12, 1.33, r"$B_{cd}(w)$", fontsize=13, ha="left", va="top")

ax.set_xlim(-w_max - 0.15, w_max + 0.35)
ax.set_ylim(-0.05, 1.45)

# strip the default frame; we drew our own axes
for s in ax.spines.values():
    s.set_visible(False)

fig.tight_layout()
fig.savefig("fig-10-1_bandpass_window.png", dpi=160, bbox_inches="tight")
print("wrote fig-10-1_bandpass_window.png")
