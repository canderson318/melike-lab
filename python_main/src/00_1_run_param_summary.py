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
from src.lib.tools import *

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)

#\\\
# ––– Summarize parameters for all patients —–––––––––––––––––––––––––––––––––
#\\\

# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
real_output_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
local_output_path = WD/'outputs'
# patients = ["SM001","SM002","SM012", "SM020", "SM022"]
patients = ["SM001"]

# aggregate summaries
cmd = 'srcmatlab src/lib/T1D_moving_window_smoother_param_summary.m'

for pat in patients:
    print(f"\n\n****** Compiling Param Summary for {pat} ****** ")
    runPatient(command = cmd, settings = {"pat": pat})
        
