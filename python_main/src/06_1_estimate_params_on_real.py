import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.lib.tools import *
from src.lib.GlucoseSim import *
import pickle as pkl
from scipy import optimize
from scipy.stats import norm
from numdifftools import Hessian
from scipy.linalg import inv

#\\\\
# ── Setup ─────────────────────────────────────────────────────────────────────
#\\\\

WD = Path("/Users/canderson/Documents/school/local-melike-lab/melike-lab/python_main")
os.chdir(WD)

# Mount OneDrive if not already mounted
mount_odrive()

# link outputs it to this directory
project_dir = Path("/Users/canderson/odrive/home/melike-rotation/project001")
output_real_path = project_dir/'outputs/'
os.listdir(output_real_path)
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)

# set output directories
out_dir = local_output_path/ "06_1"
out_dir.mkdir(exist_ok = True)

#\\\\\
#\\\\\
# –––– Load Patient Data
#\\\\\
#\\\\\
patient_info = pd.read_csv(project_dir / "Tidepool_Exports/Tandem_Tidepool_Deidentified.csv")

pat = "SM002"
pat_info = {k:v.get(1) for k,v in dict(patient_info[patient_info.patient_id == pat]).items()}

bg, ins, nut = loadPatientData(pat)

# match matlab t1dm case
nut = nut[nut.carbs.notna() & (nut.carbs != 0)].copy()
nut.carbs *= 1000  # g → mg

ins_basal = ins[ins.Basal.notna()][['time', 'Basal']].copy()
ins_bolus = ins[ins.Bolus.notna() & (ins.Bolus != 0)][['time', 'Bolus']].copy()
bg = bg[bg.bg.notna()].copy()

# unified time index — union of all event times
all_times = sorted(set(bg.time) | set(ins_basal.time) | set(ins_bolus.time) | set(nut.time))
dat = pd.DataFrame({'time': all_times})

# bg: NaN where no measurement (don't interpolate)
dat = dat.merge(bg[['time', 'bg']], on='time', how='left')
# basal: piecewise constant — forward-fill between rate changes
dat = dat.merge(ins_basal, on='time', how='left')
dat['Basal'] = dat['Basal'].ffill().fillna(0)
# bolus: event-based — 0 between events
dat = dat.merge(ins_bolus, on='time', how='left')
dat['Bolus'] = dat['Bolus'].fillna(0)
# carbs: event-based — 0 between meals
dat = dat.merge(nut[['time', 'carbs']], on='time', how='left')
dat['carbs'] = dat['carbs'].fillna(0)

dat['time_hr'] = dat['time'] / 60
dat['time_day'] = dat['time_hr'] / 24

par_summ = param_summary(Path("outputs/00_0/Gb0_1000_gamma0.012_0.035_sigma0_100_a0.01_0.1_b0.01_0.1_beta15_130/"),pat="SM001",window_set = "moving_window_of_360mins_by_180mins")
window = (1001,1361)
pars = par_summ[(par_summ.interval_start == window[0]) &  (par_summ.interval_end == window[1])]

SUB = dat[(dat.time>=window[0]) & (dat.time<= window[1])]

fig, ax = plt.subplots(3,1,figsize = (15,10))
sns.lineplot(SUB, x="time_day", y="bg", ax=ax[0])
ax[0].set_ylabel("Glucose", color="steelblue")
ax[0].tick_params(axis='y', labelcolor="steelblue")
ax[0].set_xlim((SUB.time_day.min(), SUB.time_day.max()))
ax[0].set_xlabel("")

ax[1].step(SUB['time_day'], SUB['carbs'], color='orange', label='Carbs',alpha = .8)
ax[1].set_ylabel('Carbs', color="orange")
ax[1].tick_params(axis='y', labelcolor='orange')
ax[1].set_xlim((SUB.time_day.min(), SUB.time_day.max()))
ax[1].set_xlabel("")

ax[2].step(SUB['time_day'], SUB['Basal'], color='green', label='Basal', alpha = .8)
ax[2].set_ylabel('Basal Insulin', color="green")
ax[2].tick_params(axis='y', labelcolor='green')

ax[2].scatter(SUB['time_day'], SUB['Bolus'], color='green', label='Bolus', alpha = 1,  s= 50)
ax[2].set_ylabel("Bolus Insulin", color = "green")
ax[2].tick_params(axis = "y", labelcolor = "green")

ax[2].set_xlabel("")
ax[2].set_xlim((SUB.time_day.min(), SUB.time_day.max()))
ax[2].legend()

[ax[i].xaxis.set_major_locator(plt.MultipleLocator(.5)) for i in range(3)]
plt.suptitle(f"{pat} Glucose evolution", y=1, fontsize=14)
plt.tight_layout()
plt.xlabel("Day")
plt.savefig(out_dir/f"{pat}_glucose_evolution.png")


# \\\\
# \\\\
# –––– Compare simulation to matlab estimates
# \\\\
# \\\\

BV = estimate_blood_volume_dL(age = pat_info.get("age"), weight_kg = pat_info.get('weight_kg'), height = pat_info.get("height_cm"), sex = pat_info.get("sex")) * 30
par_dict = {key:val.item() for key,val in dict(pars).items()}
sim = GlucoseSim(
    t = SUB.time.values, Gb = par_dict.get("Gb"),sigma = par_dict.get("sigma"),gamma = par_dict.get("gamma"),
    Gmeal = SUB['carbs'].values/BV,b_meal=par_dict.get("b"), a_meal = par_dict.get("a"), 
    Ibolus = SUB['Bolus'].values/BV, beta=par_dict.get("beta"), a_ins = 0.01, b_ins = 0.02)
sim.run(num_sim=10)

fig, ax = plt.subplots(figsize = (10,10))
ax.plot(SUB['time'], SUB['bg'])

df = pd.DataFrame(sim.results)
mean_G = df.mean(axis=1)
std_G  = df.std(axis=1)
if ax is None:
    fig, ax = plt.subplots()
color = kwargs.pop('color', 'steelblue')
ax.plot(sim.t,mean_G, label='mean', color=color, **kwargs)
ax.fill_between(sim.t, mean_G - std_G, mean_G + std_G,
                alpha=0.3, color=color, label='±1 SD')
ax.legend(loc = "upper left", framealpha = 1,bbox_to_anchor=(1.15, 1))


#\\\\
#\\\\
# ––– Optimize
#\\\\
#\\\\
def neg_log_likelihood(params, t, G_obs, Gmeal, Ibolus):
    Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta = params
    # optimize b as function a
    #++ b_meal = a + np.exp(log_d_meal)
    b_meal = a_meal + np.exp(log_d_meal)
    b_ins = a_ins + np.exp(log_d_ins)

    # check gamma not near anything
    EPS = 1e-5
    if any(abs(gamma - x) < EPS for x in [a_meal, b_meal, a_ins, b_ins]):
        return np.inf

    ll = 0
    # integrate functions
    for k in range(len(t) - 1):
        hk = t[k+1] - t[k]
        # integrate meals up to k
        mk = GlucoseSim.integrate_meal(k,Gmeal, a_meal, b_meal, t, gamma)
        # integrate insulin function
        ik = GlucoseSim.integrate_bolusInsulin(k,Ibolus, a_ins, b_ins, t, gamma)
        # estimated glucose given params, meals, and insulin
        mu = GlucoseSim.nextG_raw(Gb, gamma, hk, G_obs[k], mk, beta, ik)
        # estimated fluctuation around mu
        sigma_k = GlucoseSim.fluctuation(sigma, gamma, hk)
        # log liklihood of mu (next prediction of next G) given compared to observed with estimated sigma
        ll_k = norm.logpdf(G_obs[k+1], loc=mu, scale=sigma_k)
        if not np.isfinite(ll_k):
            print(f"k={k}, mu={mu:.3f}, sigma_k={sigma_k:.3f}, G_obs={G_obs[k+1]:.3f}")
            return np.inf
        ll += ll_k

    return -ll

start,stop = 5*60*24 + (60*4), 5*60*24+ (60*10)
SUB = dat[(dat.time > start) & (dat.time < stop)]
t = SUB.time.values
G_obs = SUB.bg.values
Ibolus = SUB.Bolus.values
Gmeal = SUB.carbs.values

#------------------Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta
start_points   = [100,     10,   .06,    .01,         -1,    .01,        -1,   100]

#------------------Gb,   sigma,   gamma,  a_meal, log_d_meal,    a_ins,  log_d_ins,    beta
# bounds =    [(50,250),(0.1,50),(1e-4,1),(1e-4,1),   (-100,5), (1e-4,1),   (-100,5), (0,500)]
bounds =    [(50,250), (0.1,100),   (0,1),   (0,1),   (-100,0), ((0,1)),   (-100,0), (0,500)]

# optimize.show_options(solver = 'minimize', method = 'Nelder-Mead')
optimize.show_options(solver = 'minimize', method = 'CG')

result = optimize.minimize(
    neg_log_likelihood,
    x0=start_points,
    args=(t, G_obs, Gmeal, Ibolus),
    # bounds=bounds,
    # method='CG', options = {'maxiter':2_000}
    method='Nelder-Mead', options = {'maxiter':2_000, 'adaptive': True}
    # method='L-BFGS-B', options = {'maxiter':2_000,'maxfun':10_000}
)

map_params = result.x

Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta = map_params
b_meal = a_meal + np.exp(log_d_meal)
b_ins = a_ins + np.exp(log_d_ins)

new_sim = GlucoseSim(Gb = Gb, sigma = sigma, gamma = gamma, 
                     a_meal = a_meal, b_meal = b_meal,
                     a_ins =  a_ins, b_ins = b_ins, beta = beta,
                     Gmeal = Gmeal, Ibolus = Ibolus , t  = t)

new_sim.run(num_sim = 30)
new_sim.t.max()
t.max()

fig, ax = plt.subplots(figsize = (15,10))
sns.lineplot(SUB, x = "time", y = "bg", color = 'steelblue', ax = ax, label = "Actual")
new_sim.plot(ax = ax, color = "orange")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Glucose", color = "steelblue")
plt.savefig(out_dir/f"estimate.pdf")

#\\\\
#\\\\
#–––– Estimate confidence around MAP
#\\\\
#\\\\


# Hessian at MAP:
#++ each parameter's second derivative w.r.t. every other parameter
#+++ small along diagonal: low curvature, more bowl like, less identifiable
#+++ large along diagonal: high curvature, ll sensitive to parameter, well identified
#+++ large off diag: parameters are correlated in ll space, changing one, the other will compensate
#+++ small off diag: parameters are uncorrelated in ll


H = Hessian(lambda p: neg_log_likelihood(p, t, G_obs, Gmeal, Ibolus),
            method = "central", order = 2, step = 1e-3)(
                map_params
                )
H = pd.DataFrame(H, index = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"], columns = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"])

# pkl.dump(H, open (out_dir / "hessian.pkl", 'wb'))
# H = pkl.load(open (out_dir / "hessian.pkl", 'rb'))

# posterior covariance and std errors
cov = inv(H)
std_errors = np.sqrt(np.diag(cov))
err  = pd.DataFrame({"SE": std_errors, "MAP":map_params}, index = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"])
err


