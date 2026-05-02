##
#+ Run a smoother script on a patient specified in `settings.json` which is read into the matlab script.
#+
##
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib.tools import *
import json
from concurrent.futures import ThreadPoolExecutor

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)

mount_odrive()
output_real_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
local_output_path = WD/'outputs'
sym_link(output_real_path, local_output_path)

out_dir = local_output_path / "00_0/"
out_dir.mkdir(exist_ok=True)

#\\\
# ––– Run Smoother on Patient —–––––––––––––––––––––––––––––––––
#\\\


settings = {
    # "output_dir": "/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports/",
    "smoother":{
        "bounds":{
            # "Gb":(0,1000), # original
            "Gb":(0,700),

            # "gamma": (.0001,.5), # original
            # "gamma": (.012,.035),
            "gamma": (.01,.05),

            # "sigma":(0,100), # original
            "sigma":(0,60),

            "a":(.001,.1), # original

            "b":(.001,.1), # original

            # "beta": (10,120) # original
            # "beta": (15,200)
            "beta": (10,200)

            },
        "window":{
            # "size":360,
            "size":60*8,

            "stride":180,

            # "num_windows": "all"
            # "num_windows": 8,
            "max_time": 48*60
        }
        }
    }

bounds_string = "_".join([f"{x}{y[0]}_{y[1]}" for x,y in  settings['smoother']["bounds"].items()])
script_output_dir = out_dir / bounds_string
script_output_dir.mkdir(exist_ok = True)

# save settings to dir
json.dump(settings, open(script_output_dir / "settings.json", 'w'), indent = 4)

settings.update({"output_dir": str( script_output_dir ) })


patients = ['SM001', 'SM012', 'SM020', 'SM022']

cmd = "srcmatlab src/lib/T1D_moving_window_smoother.m"

# Run command on each patient
for pat in patients:
    settings.update({"pat":pat})
    # source matlab script
    print(f"\n\n****** Running Smoother on {pat} ****** ")
    runPatient(command = cmd, settings = settings)


# \\\\
# \\\\
# –––– Summarize params
# \\\\
# \\\\


# in_dir = Path('/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports')
top_in_dir = local_output_path/"00_0"
in_dirs = [top_in_dir / x for x in os.listdir(top_in_dir)]

# in_dir = in_dirs[1];
for in_dir in in_dirs:
    window_sets = [x for x in os.listdir(in_dir) if x.find("moving_window_of") !=-1]
    summs = []
    print(in_dir)
    # pat = patients[0]; window_set = "moving_window_of_360mins_by_180mins"
    for i,window_set in enumerate(window_sets):
        print("\t\t",window_set)
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