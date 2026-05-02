import numpy as np
import seaborn as sns
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import plotly.express as px
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib.tools import *
import re
import json
import re
from sklearn.decomposition import PCA
import scipy
import statsmodels.formula.api as smf
import itertools
#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)


# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
output_real_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
os.listdir(output_real_path)
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)

# set output directories
out_dir = local_output_path/ "04_1"
out_dir.mkdir(exist_ok = True)



#\\\\
#\\\\
# ── Load data ─────────────────────────────────────────────────────────────────
#\\\\
#\\\\

# get starttimes
project_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001")
patient_info = pd.read_csv(project_dir / "Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")

# fix timing
pats = patient_info.patient_id.dropna().unique()
start_dt_time = {pat: pd.Timestamp(patient_info.start_time[patient_info.patient_id == pat].iloc[0]) for pat in pats}

# load smoother results and settings
def parse_bounds(s, prepend = ''):
    parts = re.findall(r'([a-zA-Z]+)([\d.]+)_([\d.]+)', s)
    return {f"{prepend}{name}": (float(lo), float(hi)) for name, lo, hi in parts}

in_dirs = [ local_output_path / "00_0" / D for D in os.listdir(local_output_path/"00_0")]

summs = []
# in_dir = in_dirs[0]
for in_dir in in_dirs:
    # load settings
    smoother_settings = json.load(open(in_dir / "settings.json", 'r'))

    # load summary
    summary = pd.read_csv(in_dir/"param_summary.csv")

    summary['interval_size'] = summary['interval_end'] - summary['interval_start']
    # make shorter label
    summary['label'] = [re.sub(r".*moving_window_of_|_9_am_p.*", "", x) for x in summary.window_name]

    # add start_time to all timing
    summary["datetime"] = (
        summary
        ['pat']
        .map(start_dt_time) # make series of start_times for each pat
        + pd.to_timedelta(summary['interval_start'], unit  = "m") # add interval_start to each real start time
    )

    # get day of episode
    summary['delta_day'] = (
        summary
        .groupby("pat")
        ["datetime"]
        .transform(lambda x: (x.dt.normalize() - x.dt.normalize().min()).dt.days )
    )
    summary['tod'] = summary['datetime'].dt.time
    summary['hrod'] = summary.datetime.dt.hour
    div = 8
    summary['rnd_hrod'] = (summary.hrod / div).round(1) * div
    summary['minod'] = summary.datetime.dt.minute + summary.datetime.dt.hour*60
    # reorder
    summary.sort_values(by = ["pat","window_name", "datetime"],inplace = True)
    # add settings
    bounds_dict = parse_bounds(in_dir.name, "bounds_")
    new_cols = pd.DataFrame([bounds_dict ] * summary.shape[0])
    summary = pd.concat([summary, new_cols], axis = 1)
    summary['par_string'] = in_dir.name
    summs.append(summary)

# summary.columns
raw_SUMM = pd.concat(summs, ignore_index = True)


# \\\\
# \\\\
# ––– Plot param x hour
# \\\\
# \\\\

bounds_cols = [x for x in raw_SUMM.columns if x.find("bounds_") != -1]
cols = [x.replace("bounds_", "") for x in bounds_cols] + ['rmse']

# in_dir = in_dirs[0]; window_name  = raw_SUMM.window_name.unique()[0]
for in_dir in in_dirs:
    for window_name in raw_SUMM.window_name.unique():
        plot_dir = out_dir/window_name
        plot_dir.mkdir(exist_ok = True)

        # subset for specific window
        SUB = raw_SUMM.loc[(raw_SUMM.window_name == window_name )& (raw_SUMM.par_string == in_dir.name)].copy()

        if SUB.empty:
            continue

        # \\\\
        # ––– Plot param x hour
        # \\\\
        fgsz=(20,20)
        fig, ax = plt.subplots(round(len(cols)/2),2, figsize=fgsz)
        axes = ax.flatten()
        for i,col in enumerate(cols):
            bound_col = f'bounds_{col}'
            if  bound_col in bounds_cols:
                lwr,uppr = SUB[bound_col].unique()[0]
                rng = uppr-lwr
                axes[i].axhspan(lwr, uppr, color="grey", alpha=0.1)
                axes[i].axhline(y = lwr, color="grey", alpha=0.5, linestyle = "--")
                axes[i].axhline(y = uppr, color="grey", alpha=0.5, linestyle = "--")
                axes[i].set_ylim((lwr - .1*rng, uppr +.1*rng))
            plt_SUB = SUB[~SUB[col].isna()]
            sns.lineplot(SUB, x = 'hrod', y = col, ax = axes[i], hue = 'pat')
            axes[i].set_xlabel("Interval start TOD")
            axes[i].set_ylabel(col.upper())
            axes[i].set_title(f"{col.upper()}")
            handles, labels = axes[i].get_legend_handles_labels()
            axes[i].legend(handles, labels, loc='upper right', title=' Pat', bbox_to_anchor=(1.11, 1.11))
            axes[i].xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: str(x)+":00"))
            plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right")
        fig.tight_layout(rect = [0,0,1, 0.97])
        fig.suptitle(f"Parameters throughout the day\n{in_dir.name}", fontsize = 14, y = .999)
        plt.savefig(plot_dir / f'{in_dir.name}.pdf', bbox_inches = "tight")


# fig = px.scatter_3d(SUMM[SUMM.window_name == window_nm], x='gamma', y='beta', z='hrod', opacity=0.4, color= 'par_string')
# fig.write_html("out.html")

features = ['Gb', 'gamma', 'sigma', 'a', 'b', 'beta']
X  = raw_SUMM.copy()[features]
X = X.sub(X.mean(axis = 0), axis = 1).div(X.std(axis = 0), axis =1)

n_components = 5
pca = PCA(n_components=n_components)
pca.fit(X)
components = pd.DataFrame(pca.components_.T,index=features,columns=[f'PC{i+1}' for i in range(n_components)])
pcs = pd.DataFrame(pca.transform(X), columns = [f"PC{i+1}" for i in range(n_components)])

SUMM = raw_SUMM.copy()
SUMM = SUMM.drop(columns=pcs.columns, errors='ignore')
SUMM = pd.concat([SUMM, pcs], axis = 1)

categories = SUMM.par_string.astype('category')
n = len(categories.cat.categories)
cmap = mcolors.ListedColormap(plt.cm.tab10.colors[:n])
colors = categories.cat.codes

# SUMM['day_code'] = np.select(
#       [
#         (SUMM.hrod >= 12) & (SUMM.hrod <= 18),
#         (SUMM.hrod > 18),
#         (SUMM.hrod < 6),
#     ],
#       [0, 1,2]
#   )# afternoon, night, morning
# SUMM['day_category'] = np.array(["afternoon", "night", "morning"])[SUMM.day_code]

mask = [True] * SUMM.shape[0]
# mask = (SUMM.hrod >=6) & (SUMM.hrod <=18)
hue_var = 'fval'
sns.pairplot(SUMM.loc[mask,:], vars=pcs.columns[:3], hue=hue_var,diag_kind = 'kde', palette = 'viridis')

_, ax = plt.subplots()
sns.kdeplot(SUMM, x='PC1', hue=hue_var, fill=True, alpha=0.4, ax=ax)
legend = ax.get_legend()
handles, labels = legend.legend_handles, [t.get_text() for t in legend.get_texts()]
legend.remove()
ax.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='lower right', borderaxespad=0)
plt.tight_layout()


corr = np.array([scipy.stats.spearmanr(SUMM.loc[mask,hue_var].values, pcs.loc[mask].iloc[:,i])[0] for i in range(n_components)]).ravel()
plt.figure()
hm_args = dict(cmap = 'RdBu_r', vmin=-1, vmax=1, center=0)
sns.heatmap(corr.reshape(1, -1),**hm_args)


df = SUMM.copy()
df["interval_size_hr"] = (df["interval_size"] / 60).round()
x_var = 'Gb'
y_var = 'rmse'
model = smf.ols(f"{y_var} ~ {x_var}*par_string", data=df).fit()
model.summary()

colors = plt.cm.tab10.colors
x_range = np.linspace(df[x_var].min(), df[x_var].max(), 200)

fig, ax = plt.subplots()
for i, par in enumerate(df.par_string.unique()):
    mask = df.par_string == par
    ax.scatter(df.loc[mask, x_var], df.loc[mask, y_var], alpha=0.4, color=colors[i])
    pred_df = pd.DataFrame({x_var: x_range, "par_string": [par]*len(x_range)})
    pred_summary = model.get_prediction(pred_df).summary_frame(alpha=0.05)
    ax.fill_between(x_range, pred_summary['mean_ci_lower'],  pred_summary['mean_ci_upper'], color = colors[i], alpha = .1)
    ax.plot(x_range, pred_summary['mean'], color=colors[i], label=par)

ax.legend(bbox_to_anchor = (1.01,1))
ax.set_xlabel(x_var.upper())
ax.set_ylabel(y_var.upper())

df['xx'] = [str(x) + "--" + y for x,y in zip( df.interval_size_hr, df.par_string)]
_, ax = plt.subplots()
sns.boxplot(df.sort_values("xx"), y = 'rmse',hue = "xx", ax = ax)
ax.legend(bbox_to_anchor = (1,1))


#\\\\
#\\\\
# ––– Predict rmse from params
#\\\\
#\\\\


vars = ['Gb','gamma','sigma','a','b','beta']
model = smf.ols(f"rmse ~ { ' + '.join(vars)} + par_string + interval_size_hr", data=df ).fit()
model.summary()

pred_summary = model.get_prediction(df).summary_frame(alpha = .05)

order = pred_summary['mean'].argsort().values
_, ax = plt.subplots()
ax.plot(pred_summary['mean'].iloc[order].values, df.rmse.iloc[order].values, 'o', alpha=0.3)
ax.plot(pred_summary['mean'].iloc[order].values, pred_summary['mean'].iloc[order].values, 'k--')
ax.fill_between(  pred_summary['mean'].iloc[order].values,  pred_summary['mean_ci_lower'].iloc[order].values,  pred_summary['mean_ci_upper'].iloc[order].values, alpha=0.3)
ax.set_xlabel("Predicted RMSE")
ax.set_ylabel("Actual RMSE")
plt.show()

resid =  df.rmse - pred_summary['mean']
x = np.arange(len(df.rmse))
sd = resid.std()

fig, ax = plt.subplots()
ax.fill_between(x = x, y1 = [sd]*len(x), y2= [-sd]*len(x), alpha = .3, color = 'orange')
ax.axhline(y = 0, linestyle = "--")
ax.scatter(x = x, y =resid)


# # \\\\
# # \\\\
# # –––– Generate posteriors of params
# # \\\\
# # \\\\

# import patsy
# import emcee

# # design matrix
# pars = [ "Gb", "gamma",  "sigma", "a", "b", "beta"]
# y, X = patsy.dmatrices(f"rmse ~ {' + '.join(pars)}",
#                         data=df, return_type='dataframe')
# X, y = X.values, y.values.ravel()

# def log_lik(theta):
#     coefs, log_sigma = theta[:-1], theta[-1]
#     sigma = np.exp(log_sigma)
#     resid = y - X @ coefs
#     return -0.5 * np.sum((resid/sigma)**2) - len(y)*log_sigma  # log normal

# ndim = X.shape[1] + 1  # coefs + log_sigma
# nwalkers = 32
# p0 = np.random.randn(nwalkers, ndim) * 0.1

# sampler = emcee.EnsembleSampler(nwalkers, ndim, log_lik)
# sampler.run_mcmc(initial_state=p0, nsteps =  10000)
# samples = sampler.get_chain(discard=5000, flat=True)


# labels = ['intercept'] + pars + ['log_sigma']
# fig, ax = plt.subplots(samples.shape[1], 2, figsize=(20, 20))
# coefs = []
# for i in range(samples.shape[1]):
#     Z = samples[:, i]
#     x = np.linspace(Z.min(), Z.max(), 200)
#     kde = scipy.stats.gaussian_kde(Z)
#     y = kde(x)
#     y = y/y.max()
#     map = x[np.argmax(y)]
#     P = .95
#     lwr_p, uppr_p = (1 - P)/2, P + (1 - P)/2
#     lwr, uppr = np.quantile(Z, q = [lwr_p, uppr_p])
#     coefs.append({"par": labels[i], "MAP":map, f"CI_{lwr_p:.3f}": lwr, f"CI_{uppr_p:.3f}":uppr})

#     counts, _, _ = ax[i,1].hist(Z, bins=50, alpha=0.5, color = "black")
#     mask = (x >= lwr) & (x <= uppr)
#     ax[i,0].plot(x, y*counts.max())
#     ax[i,0].fill_between(x[mask], y[mask]*counts.max(), alpha=0.3, color="blue")
#     ax[i,0].set_xlabel(labels[i])
#     ax[i,0].axvline(np.median(Z), linestyle = '--', color = 'green', label = "Median")
#     ax[i,0].axvline(map, linestyle = '-',linewidth = 5, color = 'orange', label = 'MAP')
#     ax[i,0].legend()
#     ax[i,1].plot(Z)
#     ax[i,1].set_ylabel(labels[i])
#     ax[i,1].set_ylim((Z.min(), Z.max()))
# plt.tight_layout()
# plt.show()

# MCMC_summ = pd.DataFrame(coefs)

# MAP = MCMC_summ['MAP'].values
# pred = X @ MAP[:-1].T

# plt.scatter(y,pred); plt.xlabel("Actual"); plt.ylabel("Predicted")
# plt.axline((0, 0), slope=1, color='k', linestyle='--')

