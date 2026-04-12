import numpy as np
import os
from pathlib import Path
import subprocess as sp
from typing import Optional, Tuple, Dict, Any

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


def setPatient(pat):
    """Write patient to text file for reading in matlab scripts."""
    with open("./PAT.txt", 'wt', encoding = "UTF+8",) as f:
        f.write(pat)
    return pat
    
def runPatient(command: str,pat: Optional[str]=None):
    """Run window smoother on patient"""
    if pat is not None:
        setPatient(pat)
    sp.run(['zsh', '-i', '-c', command],executable='/bin/zsh')
    