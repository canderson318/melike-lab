import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib.tools import *
import re
from scipy.io import matlab

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\
    
WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)


# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
output_real_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
os.listdir(output_real_path)
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)

# set output directories
out_dir = local_output_path/ "04"
out_dir.mkdir(exist_ok = True)


#\\\\
#\\\\
# ── Load data ─────────────────────────────────────────────────────────────────
#\\\\
#\\\\

# load summaries
pats = [x for x in os.listdir(local_output_path/"param_summary") if x.startswith("SM")]
summary = pd.concat([pd.read_csv(local_output_path/"param_summary"/pat/"summary.csv") for pat in pats])

# make shorter label
summary['label'] = [re.sub(r".*moving_window_of_|_9_am_p.*", "", x) for x in summary.window_name]

# get starttimes
project_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001")
patient_info = pd.read_csv(project_dir / "Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")

# fix timing
start_dt_time = {pat: pd.Timestamp(patient_info.start_time[patient_info.patient_id == pat].iloc[0]) for pat in pats}

# add start_time to all timing
summary["datetime"] = (
    summary['pat'].map(start_dt_time) # make series of start_times for each pat
    + pd.to_timedelta(summary['interval_start'], unit  = "m") # add interval_start to each real start time
)

# get day of episode
summary['delta_day'] = (
    summary
    .groupby("pat")
    ["datetime"]
    .transform(lambda x: (x.dt.normalize() - x.dt.normalize().min()).dt.days )
)


summary['tod'] = summary['datetime'].dt.time
summary['hrod'] = summary.datetime.dt.hour
div = 8
summary['rnd_hrod'] = (summary.hrod / div).round(1) * div
summary['minod'] = summary.datetime.dt.minute + summary.datetime.dt.hour*60

# reorder
summary.sort_values(by = ["pat","window_name", "datetime"],inplace = True)

# subset for specific window    
window_nm = "moving_window_of_360mins_by_180mins_9_am_pm"
window_set = summary.label[summary.window_name == window_nm].iloc[0]
interval_width, interval_stride = [int(x) for x in re.findall(r'\d+',window_set)]
SUB = summary.loc[summary.label == window_set].copy()

eps = SUB.gamma[SUB.gamma != 0].min()
SUB['log_gamma'] = np.log(SUB.gamma + eps)

# load settings
smoother_dir = project_dir/f"Tidepool_Exports/{window_nm}/{summary.pat.values[0]}"
first_window = os.listdir(smoother_dir)[0]
ml_obj = matlab.loadmat(smoother_dir/first_window/"my_settings.mat")
settings = ml_obj['my_settings']
par_names = np.array([str(x[0]) for x in settings['smoother'].item()["theta_est_names"][0][0][0]])
lwr, uppr  = settings["smoother"].item()['lower_bounds'][0][0].ravel() , settings["smoother"].item()['upper_bounds'][0][0].ravel()
par_bounds = {str(x):(lwr[i],uppr[i]) for i,x in enumerate(par_names)}

#\\\\
#\\\\
# ––– Plot param x hour boxplots styled by pat
#\\\\
#\\\\

cols = list(par_bounds.keys()) + ['rmse','log_gamma']

LINE = SUB.copy()
for col in cols:
    new_var = col+"_median"
    LINE[new_var] = (
        LINE.groupby(["pat", "rnd_hrod"])[col].transform("median",)
    )

LINE = LINE.filter(regex = r"_median|pat|rnd_hrod")

fgsz=(10,30)

fig, ax = plt.subplots(len(cols),1, figsize=fgsz)
axes = ax.flatten()
for i,col in enumerate(cols):
    plt_SUB = SUB[~SUB[col].isna()]
    plt_LINE = LINE[~LINE[col+"_median"].isna()]
    
    if col in par_bounds.keys():
        lwr,uppr = par_bounds[col]
        rng = uppr-lwr
        axes[i].axhspan(lwr, uppr, color="grey", alpha=0.1)
        axes[i].axhline(y = lwr, color="grey", alpha=0.5, linestyle = "--")
        axes[i].axhline(y = uppr, color="grey", alpha=0.5, linestyle = "--")
        axes[i].set_ylim((lwr - .1*rng, uppr +.1*rng))
    
    sns.boxplot(plt_SUB, x = "rnd_hrod", y = col, hue = "pat", ax = axes[i])
    sns.lineplot(plt_LINE, x = "rnd_hrod", y = col+"_median", hue = "pat", ax = axes[i])
    axes[i].set_xlabel("Interval start TOD")
    axes[i].set_ylabel(col.upper())
    axes[i].set_title(f"{col.upper()}")
    handles, labels = axes[i].get_legend_handles_labels()
    axes[i].legend(handles, labels, loc='upper right', title=' Pat', bbox_to_anchor=(1.11, 1.11))
    axes[i].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: str(x)+":00"))
    plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right")
fig.tight_layout(rect = [0,0,1, 0.97])
fig.suptitle("Parameters throughout the day", fontsize = 14, y = .999)
plt.savefig(out_dir/'param_against_HROD_boxplots.pdf', bbox_inches = "tight")

del LINE

#\\\\
#\\\\
# ––– Plot param x day styled by patient
#\\\\
#\\\\


LINE = SUB.copy()
for col in cols:
    new_var = col+"_median"
    LINE[new_var] = (
        LINE.groupby(["pat", "delta_day"])[col].transform("median")
    )
LINE = LINE.filter(regex = r"_median|pat|delta_day")

fig, ax = plt.subplots(len(cols),1, figsize = fgsz)
axes = ax.flatten()
for i,col in enumerate(cols):
    plt_SUB = SUB[~SUB[col].isna()]
    plt_LINE = LINE[~LINE[col+"_median"].isna()]

    if col in par_bounds.keys():
        lwr,uppr = par_bounds[col]
        rng = uppr-lwr
        axes[i].axhspan(lwr, uppr, color="grey", alpha=0.1)
        axes[i].axhline(y = lwr, color="grey", alpha=0.5, linestyle = "--")
        axes[i].axhline(y = uppr, color="grey", alpha=0.5, linestyle = "--")
        axes[i].set_ylim((lwr - .1*rng, uppr +.1*rng))

    sns.boxplot(plt_SUB, x = "delta_day", y = col, hue = "pat", ax = axes[i])
    sns.lineplot(plt_LINE, x = "delta_day", y = col+"_median", hue = "pat", ax = axes[i])
    axes[i].set_xlabel("Day")
    axes[i].set_ylabel(col.upper())
    axes[i].set_title(f"{col.upper()}")
    axes[i].legend(loc='upper right', title=' Pat', bbox_to_anchor=(1.11, 1.11))
fig.tight_layout(rect = [0,0,1, 0.97])
fig.suptitle("Parameters across study days", fontsize = 14, y = 1)
plt.savefig(out_dir/'param_against_deltaday_boxplots.pdf')

del LINE

#\\\\
#\\\\
# ––– Correlation heatmap 
#\\\\
#\\\\

corr_m = SUB.loc[:,cols].corr(method = 'spearman')
corr = pd.DataFrame(corr_m, columns = cols.values, index = cols.values)

fig,ax = plt.subplots(figsize = (10,10))
sns.heatmap(corr, ax = ax, cmap = 'RdBu_r', vmin = -1, vmax = 1)
plt.suptitle(f"Parameter Correlations (spearman)")
plt.tight_layout()
plt.savefig(out_dir/'param_corr_heamap.pdf')

