# python_main

Python plotting scripts for the T1D moving-window ODE parameter estimation project. Analyzes how fitted model parameters vary across time of day.

## Data

- **Source:** Tidepool CGM exports (patient `SM002`)
- **Location:** `~/odrive/home/melike-rotation/project001/Tidepool_Exports/` (rclone mount)
- **Input:** `optimization_summary.mat` per time window, output from the upstream MATLAB pipeline

## Scripts

### `src/00_run_smoother.py`
Sets active patient (`PAT.txt`) and runs the upstream MATLAB smoother (`T1D_moving_window_smoother_two_betas.m`).

### `src/01_estimate_basal_glucose.py`
Estimates basal glucose across all patients as the 40th percentile of BG readings between 02:00–05:00.

### `src/02_plot_mult_window_sizes.py`
Single-patient plots across multiple moving-window configurations. Outputs to `outputs/02/<pat>/`:
- RMSE vs. interval start time (per window set)
- Gamma vs. interval start time (per window set)

### `src/03_single_patient_params.py`
Single-patient parameter analysis for a fixed window set (`360mins_by_180mins`). Outputs to `outputs/03/<pat>/`:
- All params vs. time of day, colored by study day
- Pairplot of params
- Param correlation heatmap

### `src/04_multi_patient_params.py`
Multi-patient parameter analysis. Loads summaries for all patients and outputs to `outputs/04/`:
- All params vs. time of day (colored by study day, linetype by patient)
- Boxplots of params by study day with median overlay per patient

## Related

Upstream pipeline: `melike-lab/data_assimilation_main/projects/T1DM/T1D_moving_window_smoother_two_betas.m`
now sourcing from located in src/*.m. Path updated with in these scripts points to `melike-lab/data_assimilation_main/projects/T1DM`.

## Dependencies

### Python
```
python3 -m venv .venv
source .venv/bin/activate

pip install pandas numpy matplotlib seaborn

```

`pip list`>>

```Package   Version
------------------- ------------
appnope   0.1.4
asttokens 3.0.1
comm 0.2.3
contourpy 1.3.0
cycler    0.12.1
debugpy   1.8.20
decorator 5.2.1
et_xmlfile2.0.0
exceptiongroup 1.3.1
executing 2.2.1
fonttools 4.60.2
importlib_metadata  8.7.1
importlib_resources 6.5.2
ipykernel 6.31.0
ipython   8.18.1
jedi 0.19.2
jupyter_client 8.6.3
jupyter_core   5.8.1
kiwisolver1.4.7
matplotlib3.9.4
matplotlib-inline   0.2.1
nest-asyncio   1.6.0
numpy2.0.2
openpyxl  3.1.5
packaging 26.0
pandas    2.3.3
parso0.8.6
pexpect   4.9.0
pillow    11.3.0
pip  26.0.1
platformdirs   4.4.0
prompt_toolkit 3.0.52
psutil    7.2.2
ptyprocess0.7.0
pure_eval 0.2.3
Pygments  2.20.0
pyparsing 3.3.2
python-dateutil2.9.0.post0
pytz 2026.1.post1
pyzmq27.1.0
scipy1.13.1
seaborn   0.13.2
setuptools58.0.4
six  1.17.0
stack-data0.6.3
tornado   6.5.5
traitlets 5.14.3
typing_extensions   4.15.0
tzdata    2026.1
wcwidth   0.6.0
zipp 3.23.0
```

### Matlab
Version: 25.2.0.3150157 (R2025b) Update 4

Use zsh/bash rc functions to enable commandline matlab without gui grossness. 

```
# Add to ~/.zshrc for matlab commandline integration

# add to path
export PATH="/Applications/MATLAB_R2025b.app/bin:$PATH"
# MatLab cli repl
#+ Opens terminal session.
matlabrepl(){
    matlab -nodisplay 
}
# Matlab source function
# ++ Adds matlab script parent dir to PATH and sources it without opening figure GUI.
# ++ Usage: `srcmatlab path/to/script.m` or `srcmatlab path/to/script`
srcmatlab(){
    local abs
    abs=$(realpath "$1")
    local dir
    dir=$(dirname "$abs")
    local base    
    base=$(basename "${abs%.m}")
    matlab -nodesktop -nosplash -noFigureWindows -batch "addpath('$dir'); $base"   
}
```
