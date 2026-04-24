import pandas as pd
import numpy as np
import os
from pathlib import Path
import subprocess as sp
from typing import Optional, Tuple, Dict, Any
import json
from pprint import pprint
import scipy.io as sio
from concurrent.futures import ThreadPoolExecutor

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

def makeSettings(settings: Dict[str,Any]):
    """Write smoother settings to json for reading in matlab scripts."""
    with open("./settings.json", 'w') as f:
        json.dump(settings, f,indent = 4)
    pprint(settings)
    
def runPatient(command: str,settings: Optional[Dict] = None):
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
    """
    returns: bg,insulin,nutrition
    """
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
    

def param_summary(parent_dir , window_set, pat):
    """
    Multi-threaded I/O to load matlab param struct:
        optim_summary.fval
        optim_summary.fval_per_meas
        optim_summary.mse
        optim_summary.mse
        optim_summary.opt_params -> table with rows: Gb, gamma, sigma, a, b, beta_d, beta_n
    """
    DIR = Path(f"{parent_dir}/{window_set}/{pat}/")
    if not DIR.exists():
        print(f"\tDirectory does not exist:\t{DIR}")
        return 
    
    def _load_window(args):
        DIR, window, pat = args
        interval_start, interval_end  = window.split("_")[1:]
        mat_path = DIR/window/"optimization_summary.mat"
        
        if not mat_path.exists():
            print(f"\tFile not found:\t{mat_path}")
            return 
        
        mat = sio.loadmat( mat_path)
        fval_per_meas = mat['output']["fval_per_meas"][0][0][0][0]
        fval = mat['output']["fval"][0][0][0][0]
        mse = mat['output']['mse'][0][0][0][0]
        rmse = np.sqrt(mse)
        pars = mat['output'][0][0]['opt_params'].ravel()
        if len(pars) == 7:
            Gb, gamma, sigma, a, b, beta_d, beta_n = pars
            beta = None
        elif len(pars) == 6:
            Gb, gamma, sigma, a, b, beta = pars
            beta_d, beta_n = None, None
        else:
            raise ValueError("Not able to unpack `opt_params`. Check what parameters (two betas or one?).")
        return {"interval_start": interval_start,"interval_end": interval_end,"Gb": Gb,"gamma": gamma,
            "sigma": sigma,"a": a,"b": b,"beta": beta, "beta_d": beta_d, "beta_n": beta_n, 
            "rmse":rmse,"fval":fval,"fval_per_meas":fval_per_meas ,
            "pat": pat, "window_name": window_set}
        
    windows = os.listdir(DIR)
    args = [(DIR, w, pat) for w in windows]
    with ThreadPoolExecutor() as pool:
        rows = list(pool.map(_load_window, args))

    # filter out nones where file not found
    rows = [x for x in rows if x is not None]
    
    return pd.DataFrame(rows).sort_values("interval_start").reset_index(drop=True)
