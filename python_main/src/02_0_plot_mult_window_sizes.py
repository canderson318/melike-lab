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
from src.lib.tools import *

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\
     

    
WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)

# patient to analyze
# pat = setPatient("SM002")
pat = open("PAT.txt", 'rt').read()
print(pat)

# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
real_output_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
local_output_path = WD/'outputs'

sym_link(real_output_path, local_output_path)

for pat in os.listdir(local_output_path/"param_summary"):

    # set output directory for this script 
    out_dir = local_output_path/ "02" / pat
    out_dir.mkdir(exist_ok = True,parents = True)

    # aggregate summaries
    if False:
        runPatient(command = 'srcmatlab src/lib/T1D_moving_window_smoother_two_betas_param_summary.m')


    #\\\\
    # ── Load data ─────────────────────────────────────────────────────────────────
    #\\\\
        
    # load summaries
    summary = pd.read_csv(local_output_path/"param_summary"/pat/"summary.csv")

    # get starttime
    project_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001")
    patient_info = pd.read_csv(project_dir/"Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")
    start_time = timestr_to_mins(patient_info.start_time[patient_info.patient_id == pat].iloc[0])


    # add start_time to all timing
    summary.update(summary.filter(regex = "^interval") + start_time)

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

    # fig, ax = plt.subplots(figsize = (10,10))
    # sns.lineplot(summary.assign(gamma = summary.gamma + .001), x = "minod", y = "gamma", hue = 'delta_day', palette = 'viridis')
    # ticks = ax.get_xticks()
    # ax.set_xticklabels(mins_to_timestr(ticks), rotation=45, ha="right")
    # plt.yscale("log")
    # plt.show()


    # plt.scatter(summary.Gb, summary.gamma)
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.xlabel("Gb")
    # plt.ylabel("Gamma")