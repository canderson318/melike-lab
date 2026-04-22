import pandas as pd
import numpy as np
import os
from pathlib import Path
import subprocess as sp
from typing import Optional, Tuple, Dict, Any
import json
from pprint import pprint

def mount_odrive(force = False):
    '''Mount onedrive to home directory if not there.'''
    odrive_mountpoint = Path.home() / "odrive"
    odrive_path = odrive_mountpoint / "home"
    if not odrive_path.exists() or force:
        res = sp.run(
            ["rclone", "mount", "odrive:", str(odrive_mountpoint),
            "--vfs-cache-mode", "full", "--daemon"],
            capture_output=True
        )
        print(f"Odrive now mounted at {odrive_path}: ", res)

def sym_link(real_path, sym_path):
    if not isinstance(sym_path, Path):
        sym_path = Path(sym_path)
                                                                                                        
    if sym_path.exists() or sym_path.is_symlink():                                                       
        sym_path.unlink() if sym_path.is_symlink() else sym_path.rmdir()                                 
                                                                                                        
    return os.symlink(real_path, sym_path)   
    
    

def mins_to_timestr(minutes):
    """Convert minutes-since-midnight array to HH:MM strings for tick labels."""
    def _fmt(m):
        h = int(m // 60) % 24
        mn = int(m % 60)
        return f"{h:02d}:{mn:02d}"
    if np.ndim(minutes) == 0:
        return _fmt(float(minutes))
    return [_fmt(float(m)) for m in minutes]

def timestr_to_mins(string):
    """Undo previous function"""
    H,M = string.split(" ")[1].split(":")
    minutes = int(H)*60 + int(M)
    return minutes


def makeSettings(settings: dict):
    """Write patient to text file for reading in matlab scripts."""
    with open("./settings.json", 'w') as f:
        json.dump(settings, f)
    return settings['pat']
    
def runPatient(command: str,settings: Optional[dict] = None):
    """Run matlab script on patient"""
    if settings is not None:
        makeSettings(settings)
    pprint(json.load(open("settings.json", 'r')))
    frags = command.split(" ")
    path = [x  for x in frags if "/" in x or x.endswith((".m",".py"))][0]
    if not Path(path).exists():
        return FileNotFoundError("Error: File not found in `command`")
    return sp.run(['zsh', '-i', '-c', command],executable='/bin/zsh')

def loadPatientData(pat,data_dir = "/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports/data/"):
    path = Path(data_dir)/ pat
    files = ["bg.csv", "insulin.csv", "nutrition.csv"]
    def proc(path):
        df = pd.read_csv(path)
        df = df.astype(float)
        non_time_cols = [c for c in df.columns if c !="time"]
        mask = (~df.loc[:,non_time_cols].isna()).any(axis = 1) & ~df['time'].isna()
        df = df.loc[mask, :] # where any obs not missing
        return df
    return (proc(path/f) for f in files)
    