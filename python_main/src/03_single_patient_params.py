import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.utils.tools import *
import re
import time

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
output_real_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)


# aggregate summaries
if False:
    runPatient(command = 'srcmatlab src/T1D_moving_window_smoother_two_betas_param_summary.m')


for pat in os.listdir(local_output_path/"param_summary"):
    # set output directory
    out_dir = local_output_path/ "03"/pat
    out_dir.mkdir(exist_ok = True,parents = True)

    #\\\\
    #\\\\
    # ── Load data ─────────────────────────────────────────────────────────────────
    #\\\\
    #\\\\

    # load summaries
    summary = pd.read_csv(local_output_path/"param_summary"/pat/"summary.csv")

    # get starttime
    project_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001")
    patient_info = pd.read_csv(project_dir/"Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")

    # make shorter label
    summary['label'] = [re.sub(r".*moving_window_of_|_9_am_p.*", "", x) for x in summary.window_name]

    # fix timing
    start_dt_time = pd.Timestamp(patient_info.start_time[patient_info.patient_id == pat].iloc[0])

    # add start_time to all timing
    time_delta = [pd.Timestamp(start_dt_time + pd.Timedelta(minutes = t)) for t in summary.interval_start]
    summary['datetime'] = time_delta
    summary['delta_day'] = summary.datetime.dt.day - summary.datetime.dt.day.min()
    summary['tod'] = summary['datetime'].dt.time
    summary['hrod'] = summary.datetime.dt.hour
    summary['minod'] = summary.datetime.dt.minute + summary.datetime.dt.hour*60

    # reorder
    summary.sort_values(by = ["window_name", "datetime"],inplace = True)

    #\\\\
    #\\\\
    # ── Plot Gamma over time
    #\\\\
    #\\\\

    # subset for specific window    
    window_set = "360mins_by_180mins"
    interval_width, interval_stride = [int(x) for x in re.findall(r'\d+',window_set)]
    SUB = summary.loc[summary.label == window_set].copy()

    # var vs interval start
    cols = pd.Index(['Gb', 'gamma', 'sigma', 'a', 'b', 'beta_day', 'beta_night','rmse'])

    fig, ax = plt.subplots(len(cols),1, figsize=(10, 20))
    axes = ax.flatten()
    for i,col in enumerate(cols):
        # axes[i].hlines(y=0, xmin=0, xmax=interval_width, color='blue', linewidth = 10, alpha = .2, label = "Interval Width")
        # axes[i].hlines(y=0, xmin=0, xmax=interval_stride, color='green', linewidth = 10, alpha = .2, label = "Interval Stride")
        sns.lineplot(SUB, x = "minod", y = col, ax=axes[i], hue = "delta_day", palette = 'viridis')
        if col == 'gamma':
            axes[i].set_yscale("linear")
        axes[i].set_xlabel("Interval start TOD")
        axes[i].set_ylabel(col.upper())
        axes[i].set_title(f"{col.upper()}")
        axes[i].legend(loc = 'upper right',title = 'Day Delta', bbox_to_anchor = (1.11,1.11))
        axes[i].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: mins_to_timestr(x)))
        plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout(rect = [0,0,1, 0.97])
    fig.suptitle(f"{pat}\nParameters over time of day colored by day", fontsize = 14, y = 1)
    plt.savefig(out_dir/'param_against_interval_start.pdf')


    # # plot all params minmax normalized 
    # fig, ax = plt.subplots( figsize=(20, 10))
    # ax.hlines(y=0, xmin=0, xmax=interval_width, color='blue', linewidth = 10, alpha = .2, label = "Interval Width")
    # ax.hlines(y=0, xmin=0, xmax=interval_stride, color='green', linewidth = 10, alpha = .2, label = "Interval Stride")

    # for i,col in enumerate(cols):
    #     sns.lineplot(data=SUB.assign(Y = SUB[col].sub(SUB[col].min()) / (SUB[col].max() - SUB[col].min()) ), x="interval_start", y="Y", ax=ax, label = col)

    # ax.set_xlabel("Interval start")
    # ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: mins_to_timestr(x)))
    # plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    # ax.legend(loc = 'upper right')
    # fig.tight_layout()
    # plt.show()

    # plot each var against every other 
    plt.figure(figsize = (10,10))
    sns.pairplot(SUB.loc[:,cols])
    plt.suptitle(f"{pat}\n",y = .99)
    plt.savefig(out_dir/'param_pairplot.pdf')

    # corr heatmap
    corr_m = SUB.loc[:,cols].corr()
    corr = pd.DataFrame(corr_m, columns = cols.values, index = cols.values)

    fig,ax = plt.subplots(figsize = (10,10))
    sns.heatmap(corr, ax = ax, cmap = 'RdBu_r',)
    plt.suptitle(f"{pat}\n Parameter Correlations")
    plt.tight_layout()
    plt.savefig(out_dir/'param_corr_heamap.pdf')
