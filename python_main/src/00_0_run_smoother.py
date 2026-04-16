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

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

os.chdir("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")    

#\\\
# ––– Run Smoother on Patient —–––––––––––––––––––––––––––––––––
#\\\

# patients = ["SM022","SM020","SM012"]
patients = ["SM012"]

# command = "srcmatlab src/lib/T1D_moving_window_smoother_two_betas.m"
cmd = "srcmatlab src/lib/T1D_moving_window_smoother.m"

# check if script exists
scrpt = cmd.split(" ")[1]
if not Path(scrpt).exists():
    raise FileExistsError(f"`{scrpt}` not found")

for pat in patients:
    # patient to analyze
    tools.setPatient(pat) # feeds into matlab scripts
    double_check = open('PAT.txt', 'rt').read()
    if pat != double_check:
        raise ValueError("ERROR: current patient and logged patient do not match!")
    
    # source matlab script
    print(f"\n\n****** Running Smoother on {double_check} ****** ")
    tools.runPatient(command = cmd)