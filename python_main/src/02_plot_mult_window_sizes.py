# making plots for T1D_moving_window_smoother_two_betas.m
# plot multiple window settings params (window size and stride).
# translated to python from T1D_moving_window_smoother_two_betas_plotting.m

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.utils.tools import *

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\
    

os.chdir("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")    

# patient to analyze
pat = setPatient("SM002")

# # aggregate summaries
# sp.run(
#       ['zsh', '-i', '-c', 'srcmatlab T1D_moving_window_smoother_two_betas_param_summary.m'],
#       executable='/bin/zsh'                                                                                        
#   )
            
project_dir = Path(f"/Users/canderson/odrive/home/melike-rotation/project001/")
parent_dir = project_dir/f"outputs/param_summary/{pat}"
out_dir = parent_dir/ "01/"
out_dir.mkdir(exist_ok = True)

# Mount OneDrive if not already mounted
mount_odrive()

#\\\\
# ── Load data ─────────────────────────────────────────────────────────────────
#\\\\
    

# get starttime
patient_info = pd.read_csv(project_dir/"Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")
start_time = timestr_to_mins(patient_info.start_time[patient_info.patient_id == pat].iloc[0])

# load summaries
summary = pd.read_csv(parent_dir/"summary.csv")

# add start_time to all timing
summary.update(summary.filter(regex = "^interval") + 3)

#\\\\
# ── Plot per window set ───────────────────────────────────────────────────────
#\\\\
summary['label'] = [x.strip("moving_window_of_|_9_am_p") for x in summary.window_name]
    
# RMSE vs interval start
fig, ax = plt.subplots(figsize=(8, 4))
sns.lineplot(data=summary, x="interval_start", y="rmse", ax=ax, hue = 'label')
ax.set_xlabel("Interval start")
ax.set_ylabel("RMSE")
ax.set_title("RMSE At Different Start Times")
ticks = ax.get_xticks()
ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
ax.legend(loc = 'upper left', fontsize = "small",  bbox_to_anchor=(1, 1))
fig.tight_layout()
plt.show()

# GAMMA vs interval start
fig, ax = plt.subplots(figsize=(8, 4))
sns.lineplot(data=summary, x="interval_start", y="gamma", ax=ax, hue = 'label' )
ax.set_xlabel("Interval start")
ax.set_ylabel("Gamma")
ax.set_title("Gamma At Different Start Times")
ticks = ax.get_xticks()
ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
ax.legend(loc = 'upper left', fontsize = "small",  bbox_to_anchor=(1, 1))
fig.tight_layout()
plt.show()


#\\\\
# ── Gamma at each time of day across all window sets ─────────────────────────
#\\\\

fig, ax = plt.subplots(figsize = (10,10))
sns.lineplot(summary.assign(gamma = summary.gamma + .001), x = "minod", y = "gamma", hue = 'delta_day', palette = 'viridis')
ticks = ax.get_xticks()
ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
plt.yscale("log")
plt.show()


plt.scatter(summary.Gb, summary.gamma)
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Gb")
plt.ylabel("Gamma")