import matplotlib.pyplot as plt

# validated default categorical order (dataviz reference palette, light mode)
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def style():
    plt.rcParams.update({
        "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
        "axes.edgecolor": "#b5b4ae", "axes.labelcolor": "#52514e", "xtick.color": "#52514e",
        "ytick.color": "#52514e", "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.7,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
        "axes.titlesize": 11, "text.color": "#0b0b0b",
    })
