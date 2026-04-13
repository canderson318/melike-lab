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
from src.utils import tools

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

os.chdir("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")    

#\\\
# ––– Run Smoother on Patient —–––––––––––––––––––––––––––––––––
#\\\
# patient to analyze
pat = tools.setPatient("SM001")
print(f"Running Smoother on {open('PAT.txt', 'rt').read()}")

# source matlab script
tools.runPatient(command = "srcmatlab src/T1D_moving_window_smoother_two_betas.m")