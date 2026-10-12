"""Shared figure style for the reports: chart surface, ink, hairlines and series colours."""

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ("#2a78d6", "#eb6834")
BLUE_LIGHT, BLUE_DARK = "#86b6ef", "#1c5cab"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)
    ax.title.set_color(INK)


def hide_unused(axes, n_used):
    """Hide empty panels and show x tick labels on the panels above them."""
    ncols = axes.shape[1]
    for i, ax in enumerate(axes.flat):
        if i >= n_used:
            ax.set_visible(False)
        elif i + ncols >= n_used:
            ax.tick_params(labelbottom=True)
