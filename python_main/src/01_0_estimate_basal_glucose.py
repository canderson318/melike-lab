##
#+ Estimate Basal glucose by taking the 40th percentile of BG across all participants.
##
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess 
from src.lib.tools import *
import re
import scipy as sp

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\
    
os.chdir("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")    

data_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports/data")
dirs = os.listdir(data_dir)

frames = []                                                                                                                                                                       
for dir in dirs:                                                                                                                                                                
    dat = pd.read_csv(data_dir / dir / "bg.csv")                                                                                                                                  
    dat["patient_id"] = dir                                                                                                                                                       
    frames.append(dat)
                                                                                                                                                                                
data = pd.concat(frames, axis=0, ignore_index=True)                                                                                                                             

data["DT"] = pd.to_timedelta(data["time"], unit="m")
data["TOD"] = (pd.Timestamp("1970-01-01") + data["DT"]).dt.time
base = pd.Timestamp("1970-01-01") + data["DT"]
data["MINOD"] = base.dt.hour * 60 + base.dt.minute
data["HROD"] = base.dt.hour


sub = data[(data.HROD<5) & (data.HROD > 2)]


print(f"** 05:00 to 02:00\nMedian = {sub.bg.median():.3f}, STDEV = {sub.bg.std():.3f}, N = {sub.patient_id.nunique()}")

percs = []

for pat in sub.patient_id.unique():
    # gb = sub.bg[(sub.patient_id == pat) &( sub.time < 24*60)]
    gb = sub.bg[sub.patient_id == pat]
    if gb.size !=0:
        perc40 = np.percentile(gb, 40)
        percs.append(perc40)
        
percs = np.array(percs)

# Paper DOI: 10.3389/fnut.2023.1203899
print(f"** 40th Percentile ** \nAverage = {percs.mean():.3f}, STDEV = {percs.std():.3f}, N = {sub.patient_id.nunique()}")