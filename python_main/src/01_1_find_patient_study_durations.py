import pandas as pd
import os
from pathlib import Path
from src.lib.tools import *

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\
WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)

# load data
data_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports/data")
dirs = os.listdir(data_dir)

frames = []                                                                                                                                                                       
for dir in dirs:                                                                                                                                                                
    dat = pd.read_csv(data_dir / dir / "bg.csv")                                                                                                                                  
    dat["patient_id"] = dir                                                                                                                                                       
    frames.append(dat)
                                                                                                                                                                                
data = pd.concat(frames, axis=0, ignore_index=True)

 
#\\\\
#\\\\
# ––– Find patients with similar duration as SM001
#\\\\
#\\\\
(
    abs(( data.groupby("patient_id")['time'].agg('max') / 60/24) - 15 ) 
).sort_values()

# SM022     0.347465
# SM020     0.371806
# SM001     0.393449
# SM012     0.403252