
import numpy as np
import seaborn as sns
import matplotlib
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib.tools import *
import re
import time
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
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)

# pat = "SM002"
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
    
    # window_set = "90mins_by_90mins"
    for window_nm in summary.window_name.unique():
        
        window_set = summary.label[summary.window_name == window_nm].iloc[0]
        
        # load settings
        smoother_dir = project_dir/f"Tidepool_Exports/{window_nm}/{pat}"
        first_window = os.listdir(smoother_dir)[0]
        ml_obj = matlab.loadmat(smoother_dir/first_window/"my_settings.mat")
        settings = ml_obj['my_settings']
        par_names = np.array([str(x[0]) for x in settings['smoother'].item()["theta_est_names"][0][0][0]])
        lwr, uppr  = settings["smoother"].item()['lower_bounds'][0][0].ravel() , settings["smoother"].item()['upper_bounds'][0][0].ravel()
        par_bounds = {str(x):(lwr[i],uppr[i]) for i,x in enumerate(par_names)}
        
        # subset for specific window    
        interval_width, interval_stride = [int(x) for x in re.findall(r'\d+',window_set)]
        SUB = summary.loc[summary.label == window_set].copy()

        cols = list(par_bounds.keys()) + ["rmse"]
        
        # var vs interval start plot
        fig, ax = plt.subplots(len(cols),1, figsize=(10, 20))
        axes = ax.flatten()
        for i,col in enumerate(cols):
            
            if col in par_bounds.keys():
                lwr,uppr = par_bounds[col]
                rng = uppr-lwr
                axes[i].axhspan(lwr, uppr, color="grey", alpha=0.3)
                axes[i].axhline(y = lwr, color="grey", alpha=0.5, linestyle = "--")
                axes[i].axhline(y = uppr, color="grey", alpha=0.5, linestyle = "--")
                axes[i].set_ylim((lwr - .1*rng, uppr +.1*rng))

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
        fig.suptitle(f"Parameters over time of day colored by day\n{pat}, {window_set}", fontsize = 14, y = 1)
        plt.savefig(out_dir/f'{window_set}_param_against_interval_start.pdf')

        # plot each var against every other 
        plt.figure(figsize = (10,10))
        sns.pairplot(SUB.loc[:,cols])
        plt.suptitle(f"{pat}, {window_set}\n",y = .99)
        plt.savefig(out_dir/f'{window_set}_param_pairplot.pdf')
        
        # corr heatmap
        corr_m = SUB.loc[:,cols].corr()
        corr = pd.DataFrame(corr_m, columns = cols, index = cols)

        fig,ax = plt.subplots(figsize = (10,10))
        sns.heatmap(corr, ax = ax, cmap = 'RdBu_r',)
        plt.suptitle(f"Parameter Correlations\n{pat}, {window_set}")
        plt.tight_layout()
        plt.savefig(out_dir/f'{window_set}_param_corr_heamap.pdf')
        
        
        plt.close("all")

