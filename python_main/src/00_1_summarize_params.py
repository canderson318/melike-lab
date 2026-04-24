
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

# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
real_output_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
local_output_path = WD/'outputs'
sym_link(real_output_path, local_output_path)

# where summaries go
out_dir = WD / "param_summary"
out_dir.mkdir(exist_ok = True)

#\\\
# ––– Summarize parameters for all patients —–––––––––––––––––––––––––––––––––
#\\\
#   For each window set, e.g. moving_window_of_360mins_by_180mins, find all patients with that data and loop over their windows
#   appending to a full dataframe 
    

# in_dir = Path('/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports')
top_in_dir = local_output_path/"00_0"
in_dirs = [top_in_dir / x for x in os.listdir(top_in_dir)]

# in_dir = in_dirs[0]; 
for in_dir in in_dirs:
    window_sets = [x for x in os.listdir(in_dir) if x.find("moving_window_of") !=-1]
    summs = []

    # pat = patients[0]; window_set = "moving_window_of_360mins_by_180mins"
    for i,window_set in enumerate(window_sets):
        print(window_set)
        patients = os.listdir(in_dir/window_set)
        for j,pat in enumerate(patients):    
            print(f"\r\t{j/len(patients)*100:.0f}% patients completed")
            a_summ = param_summary(parent_dir = in_dir, window_set =  window_set, pat = pat)
            summs.append(a_summ)
            
    summs_filtered = [x for x in summs if x is not None]
    summary = pd.concat(summs_filtered, axis = 0, ignore_index = True)

    # –––– Save to csv ––––
    # COLUMNS:   interval_start interval_end Gb gamma sigma a b beta beta_d beta_n rmse fval fval_per_meas pat window_name
    summary.to_csv(in_dir/"param_summary.csv", index = False)