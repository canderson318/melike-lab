# python_main

Python plotting scripts for the T1D moving-window ODE parameter estimation project. Analyzes how fitted model parameters vary across time of day.

## Data

- **Source:** Tidepool CGM exports (patient `SM002`)
- **Location:** `~/odrive/home/melike-rotation/project001/Tidepool_Exports/` (rclone mount)
- **Input:** `optimization_summary.mat` per time window, output from the upstream MATLAB pipeline

## Scripts

### `src/01.py`
Loads optimization results across all moving-window directories and generates plots saved to `Tidepool_Exports/plots/`:

- RMSE vs. interval start time (per window set)
- Gamma vs. interval start time, colored by RMSE (per window set)
- Gamma estimates across all window sizes and start times (combined)

Translated from `src/T1D_moving_window_smoother_two_betas_plotting.m`.

## Dependencies

```
numpy, pandas, scipy, matplotlib, seaborn
```

## Related

Upstream pipeline: `melike-lab/data_assimilation_main/projects/T1DM/T1D_moving_window_smoother_two_betas.m`
