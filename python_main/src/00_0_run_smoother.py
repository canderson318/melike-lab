##
#+ Run a smoother script on a patient
#++ Saves patient specified to ./PAT.txt and reads that from teh matlab script, T1D_moving_window_smoother_two_betas.m
##
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib import tools
import json

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

os.chdir("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")    

#\\\
# ––– Run Smoother on Patient —–––––––––––––––––––––––––––––––––
#\\\

# patients = ["SM022","SM020","SM012"]
patients = ["SM001"]

# command = "srcmatlab src/lib/T1D_moving_window_smoother_two_betas.m"
cmd = "srcmatlab src/lib/T1D_moving_window_smoother.m"

settings = {
    "smoother":{
        "bounds":{
            "Gb":(0,1000),
            # "gamma": (.0001,.5), # orig
            "gamma": (.012,.035),
            "sigma":(0,100),
            "a":(.01,.1),
            "b":(.01,.1),
            # "beta": (10,120) # orig
            "beta": (15,130) 
            },
        "window":{
            "size":360,
            "stride":180
        }
        }
    }

for pat in patients:
    # patient to analyze
    settings.update({"pat":pat})
    tools.makeSettings(settings) # feeds into matlab scripts
    # source matlab script
    print(f"\n\n****** Running Smoother on {pat} ****** ")
    tools.runPatient(command = cmd)