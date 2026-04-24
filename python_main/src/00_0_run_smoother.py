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

cmd = "srcmatlab src/lib/T1D_moving_window_smoother.m"

settings = {
    # "output_dir": "/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports/",
    "smoother":{
        "bounds":{
            "Gb":(0,1000),
            
            # "gamma": (.0001,.5), # original
            "gamma": (.012,.035),
            
            "sigma":(0,100),
            
            "a":(.01,.1),
            
            "b":(.01,.1),
            
            # "beta": (10,120) # original
            "beta": (15,130)
            
            },
        "window":{
            "size":360,
            
            "stride":180,
            
            # "num_windows": "all"
            "num_windows": 8,
        }
        }
    }

bounds_string = "_".join([f"{x}{y[0]}_{y[1]}" for x,y in  settings['smoother']["bounds"].items()])
script_output_dir = out_dir / bounds_string
script_output_dir.mkdir(exist_ok = True)

# save settings to dir
json.dump(settings, open(script_output_dir / "settings.json", 'w'), indent = 4)

settings.update({"output_dir": str( script_output_dir ) })

patients = ["SM022", "SM012"]
pat = "SM012"
for pat in patients:
    settings.update({"pat":pat})
    # source matlab script
    print(f"\n\n****** Running Smoother on {pat} ****** ")
    runPatient(command = cmd)
