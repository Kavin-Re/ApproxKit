"""Matplotlib helpers with one consistent, restrained style."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Okabe-Ito colour-blind-safe palette; one fixed colour + marker per family.
FAMILY_STYLE = {
    "exact_add": ("#000000", "s"), "exact_mul": ("#000000", "s"),
    "loa": ("#0072B2", "o"), "trunc_add": ("#D55E00", "v"), "eta": ("#009E73", "^"),
    "gear": ("#CC79A7", "D"), "ama": ("#E69F00", "P"),
    "trunc_mul": ("#D55E00", "v"), "bam": ("#0072B2", "o"), "drum": ("#009E73", "^"),
    "mitchell": ("#CC79A7", "D"), "booth_approx": ("#E69F00", "P"),
}
COST_LABEL = "cost proxy: gate equivalents (GE, NAND2=1)  [NOT synthesized area]"
ERR_LABEL = {"MRED": "MRED (mean relative error distance)", "NMED": "NMED (normalised mean error distance)",
             "ER": "error rate", "MED": "MED"}


def _style():
    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 160, "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
        "axes.titlesize": 10, "axes.titleweight": "bold",
        "legend.frameon": False, "legend.fontsize": 8,
    })


def pareto_plot(rows, front, path, title, cost_key="area_ge", err_key="MRED"):
    _style()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    fams = []
    for r in rows:
        if r["family"] not in fams:
            fams.append(r["family"])
    for fam in fams:
        col, mk = FAMILY_STYLE.get(fam, ("#777777", "o"))
        pts = [r for r in rows if r["family"] == fam]
        ax.scatter([p[cost_key] for p in pts], [p[err_key] for p in pts], s=22, c=col,
                   marker=mk, alpha=0.75, linewidths=0, label=fam)
    ax.step([f[cost_key] for f in front], [f[err_key] for f in front], where="post",
            color="#333333", lw=1.2, label="Pareto front")
    ax.scatter([f[cost_key] for f in front], [f[err_key] for f in front], s=60,
               facecolors="none", edgecolors="#333333", linewidths=1.0)
    pos = [r[err_key] for r in rows if r[err_key] > 0]
    ax.set_yscale("symlog", linthresh=min(pos) / 2 if pos else 1e-6)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(COST_LABEL)
    ax.set_ylabel(ERR_LABEL.get(err_key, err_key) + "  (symlog; 0 = exact)")
    ax.set_title(title)
    ax.legend(ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def histogram_plot(metrics_list, path, title, bins=61):
    """Small multiples of the signed error distribution (probability, log y)."""
    _style()
    k = len(metrics_list)
    cols = min(3, k)
    rws = int(np.ceil(k / cols))
    fig, axs = plt.subplots(rws, cols, figsize=(3.3 * cols, 2.5 * rws), squeeze=False)
    for ax in axs.ravel():
        ax.set_visible(False)
    for ax, m in zip(axs.ravel(), metrics_list):
        ax.set_visible(True)
        c, e = m.histogram(bins=bins, normalize=True)
        ax.bar(0.5 * (e[1:] + e[:-1]), c, width=(e[1] - e[0]), color="#0072B2")
        ax.set_yscale("log")
        ax.set_title(m.unit, fontsize=8)
        ax.set_xlabel("error = approx - exact")
        ax.set_ylabel("probability")
    fig.suptitle(title, fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
