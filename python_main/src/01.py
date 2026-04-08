# making plots for /Users/canderson/Documents/school/local-melike-lab/melike-lab/data_assimilation_main/projects/T1DM/T1D_moving_window_smoother_two_betas.m
# translated to python from T1D_moving_window_smoother_two_betas_plotting.m

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path as P
import subprocess as sp
import scipy.io


# ── Setup ─────────────────────────────────────────────────────────────────────
sns.set_theme(style="whitegrid", context="notebook")

main_dir = P("/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports")
pat = "SM002"

# Mount OneDrive if not already mounted
odrive_path = P.home() / "odrive" / "home"
if not odrive_path.exists():
    sp.Popen(
        ["rclone", "mount", "odrive:", str(P.home() / "odrive"),
         "--vfs-cache-mode", "full", "--daemon"]
    )

# Find all moving_window directories
dirs = [d.name for d in main_dir.iterdir() if d.is_dir()]
window_sets = [d for d in dirs if "moving_window" in d]

f = "optimization_summary.mat"

# ── Load data ─────────────────────────────────────────────────────────────────
summaries = []

for window_set in window_sets:
    print(window_set)
    parent_data_dir = main_dir / window_set / pat

    data_names = sorted(
        [d.name for d in parent_data_dir.iterdir()
         if d.is_dir() and not d.name.startswith(".")]
    )

    rows = []
    for i, a_data_dir in enumerate(data_names):
        parts = a_data_dir.split("_")
        interval_start = float(parts[1])
        interval_end   = float(parts[2])

        mat = scipy.io.loadmat(
            str(parent_data_dir / a_data_dir / f),
            simplify_cells=True
        )
        optimal = mat["output"]["best_params"]["Optimal"]

        rows.append({
            "interval_start": interval_start,
            "interval_end":   interval_end,
            "Gb":             optimal[0],
            "gamma":          optimal[1],
            "sigma":          optimal[2],
            "a":              optimal[7],
            "b":              optimal[8],
            "beta_day":       optimal[9],
            "beta_night":     optimal[10],
            "f_val":          mat["output"]["fval_per_meas"],
            "rmse":           float(np.sqrt(mat["output"]["mse"])),
        })

        pct = round((i + 1) / len(data_names) * 100)
        print(f"{pct}%")

    summ = pd.DataFrame(rows).sort_values("interval_start").reset_index(drop=True)
    summaries.append(summ)


def mins_to_timestr(minutes):
    """Convert minutes-since-midnight array to HH:MM strings for tick labels."""
    def _fmt(m):
        h = int(m // 60) % 24
        mn = int(m % 60)
        return f"{h:02d}:{mn:02d}"
    if np.ndim(minutes) == 0:
        return _fmt(float(minutes))
    return [_fmt(float(m)) for m in minutes]


# ── Plot per window set ───────────────────────────────────────────────────────
for window_set, a_summ in zip(window_sets, summaries):
    plot_dir = main_dir / "plots" / window_set
    plot_dir.mkdir(parents=True, exist_ok=True)

    # RMSE vs interval start
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(data=a_summ, x="interval_start", y="rmse", ax=ax)
    ax.set_xlabel("Interval start")
    ax.set_ylabel("RMSE")
    ax.set_title(window_set)
    ticks = ax.get_xticks()
    ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(plot_dir / "rmse_against_interval_start.png", dpi=150)
    plt.close(fig)

    # Gamma vs interval start, colour = RMSE
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.lineplot(data=a_summ, x="interval_start", y="gamma", ax=ax, zorder=1)
    sc = ax.scatter(
        a_summ["interval_start"], a_summ["gamma"],
        c=a_summ["rmse"], s=100, cmap="viridis", zorder=2
    )
    fig.colorbar(sc, ax=ax, label="RMSE")
    ax.set_xlabel("Interval Start")
    ax.set_ylabel("Gamma")
    ax.set_title(window_set)
    ticks = ax.get_xticks()
    ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(plot_dir / "gamma_against_interval_start_rmse_filled.png", dpi=150)
    plt.close(fig)


# ── Gamma at each time of day across all window sets ─────────────────────────
cols = ["interval_start", "interval_end", "gamma"]
gammas = pd.concat([s[cols] for s in summaries], ignore_index=True)
gammas["window"] = gammas["interval_end"] - gammas["interval_start"]

fig, ax = plt.subplots(figsize=(9, 5))

for _, row in gammas.iterrows():
    ax.plot(
        [row["interval_start"], row["interval_end"]],
        [row["gamma"],          row["gamma"]],
        color="steelblue", alpha=0.4, linewidth=1
    )

norm = plt.Normalize(gammas["window"].min(), gammas["window"].max())
cmap = plt.cm.viridis

# Filled circles at start, open circles at end
sc_start = ax.scatter(
    gammas["interval_start"], gammas["gamma"],
    c=gammas["window"], s=100, cmap=cmap, norm=norm, zorder=3
)
ax.scatter(
    gammas["interval_end"], gammas["gamma"],
    c=gammas["window"], s=100, cmap=cmap, norm=norm,
    facecolors="none", linewidths=1.5, zorder=3
)

cb = fig.colorbar(sc_start, ax=ax, label="Interval Size")
cb_ticks = cb.get_ticks()
cb.set_ticklabels(mins_to_timestr(cb_ticks))

ax.set_xlabel("Time of Day")
ax.set_ylabel("Gamma")
ax.set_title("Gamma Estimate For Different Intervals")
ticks = ax.get_xticks()
ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")

plot_root = main_dir / "plots"
plot_root.mkdir(parents=True, exist_ok=True)
fig.tight_layout()
fig.savefig(plot_root / "gammas_at_diff_start_times.png", dpi=150)
plt.close(fig)
