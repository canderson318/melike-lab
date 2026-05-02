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
import scipy
from scipy import optimize
from scipy.stats import norm
from numdifftools import Hessian
from multiprocessing import Pool
from functools import partial
import emcee
import schwimmbad
import tqdm
from pprint import pprint

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
out_dir = local_output_path/ "06_0"
out_dir.mkdir(exist_ok = True)

# load simulation
sim = pkl.load(open(local_output_path/"05_0"/"sim.pkl", 'rb'))

# extract simulated patient
t_true, G_obs_true, Gmeal_true, Ibolus_true, I_true = sim.t, sim.results.mean(axis = 1), sim.Gmeal, sim.Ibolus, sim.Ibasal

plt.plot(t_true, G_obs_true)
plt.step(t_true, Ibolus_true)
plt.step(t_true, Gmeal_true)

# subset
# inds = range(len(sim.t))
inds = range(0,100)
t, G_obs, Gmeal, Ibolus, I  = t_true[inds], G_obs_true[inds], Gmeal_true[inds], Ibolus_true[inds], I_true[inds]



# \\\\
# \\\\
# –––– Estimate Parameters From Simulation
# \\\\
# \\\\

def log_likelihood(params, t, G_obs, Gmeal, Ibolus, I, bounds):

    par_dict = dict(zip(['Gb', 'sigma', 'gamma', 'a_meal', 'log_d_meal', 'a_ins', 'log_d_ins', "beta"], params))
    Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta = params

    # optimize b as function a
    #++ b_meal = a + np.exp(log_d_meal)
    b_meal = a_meal + np.exp(log_d_meal)
    b_ins = a_ins + np.exp(log_d_ins)

    par_dict.update({"b_meal":b_meal, "b_ins": b_ins})

    # ensure params correct and within bounds
    check_pars = ["Gb", "sigma", "gamma", "a_meal",'b_meal', "a_ins","b_ins", "beta" ]
    EPS = 1e-3
    # if any(abs(gamma - x) < EPS for x in [a_meal, b_meal, a_ins, b_ins] ) or \
    #     any((par_dict[x] < bounds[x][0]) | (par_dict[x] > bounds[x][1]) for x in check_pars):
    #     return - np.inf

    if any(abs(gamma - x) < EPS for x in [a_meal, b_meal, a_ins, b_ins] ):
        # print(f"a/b for meal/insulin function within {EPS:.2e} of gamma.")
        return - np.inf

    for par in check_pars:
        if par_dict[par] < bounds[par][0] or par_dict[par] > bounds[par][1]:
            # print(f"{par} = {par_dict[par]:.2f} outside of bounds ({bounds[par][0]:.3f}, {bounds[par][1]:.3f})")
            return -np.inf

    ll = 0
    # integrate meal function
    for k in range(len(t) - 1):
        hk = t[k+1] - t[k]
        # integrate meals up to k
        mk = GlucoseSim.integrate_meal(k,Gmeal, a_meal, b_meal, t, gamma)
        # integrate insulin function
        ik = GlucoseSim.integrate_bolusInsulin(k,Ibolus, a_ins, b_ins, t, gamma) + GlucoseSim.integrate_basalInsulin(k, I, t, gamma)
        # estimated glucose given params, meals, and insulin
        mu = GlucoseSim.nextG_raw(Gb, gamma, hk, G_obs[k], mk, beta, ik)
        # estimated fluctuation around mu
        sigma_k = GlucoseSim.fluctuation(sigma, gamma, hk)
        # log liklihood of mu (next prediction of next G) given compared to observed with estimated sigma
        ll_k = norm.logpdf(G_obs[k+1], loc=mu, scale=sigma_k)
        if not np.isfinite(ll_k):
            print(f"k={k}, mu={mu:.3f}, sigma_k={sigma_k:.3f}, G_obs={G_obs[k+1]:.3f}")
            return - np.inf
        ll += ll_k
    return ll

def neg_log_likelihood(params, bounds):
    return - log_likelihood(params, t, G_obs, Gmeal, Ibolus, I, bounds)

start_points = {'Gb': 145,
            'sigma':     10,
            'gamma':     .07,
            'a_meal':    .04,
            'log_d_meal':-4,
            'a_ins':     .01,
            'log_d_ins': -4,
            'beta':      100
            }
bounds =  {'Gb':      (100, 200),
        'sigma':      (5, 15),
        'gamma':      (0.06, 0.08),
        'a_meal':     (0.01, 0.1),
        # 'log_d_meal': (-100, 0),
        'b_meal':  (0.001, 0.1),
        'a_ins':      (0.001, 0.1),
        # 'log_d_ins':  (-100, 0),
        'b_ins':  (0.001, 0.1),
        'beta':       (50, 200)}


nwalkers = int(len(bounds)*2.5)
ndim = len(start_points.values())
nsteps=1000
scales = np.array([10, 1, 0.005, 0.005, 0.5, 0.001, 0.5, 10])  # rough scale per param
p0 = np.array(list(start_points.values())) + np.random.randn(nwalkers, ndim) * scales

# log_likelihood(p0[0,:], t, G_obs, Gmeal, Ibolus, I, bounds)

log_likelihood_wrapped = partial(log_likelihood, t=t, G_obs=G_obs, Gmeal=Gmeal, Ibolus=Ibolus, I=I, bounds = bounds)

with schwimmbad.MultiPool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood_wrapped, pool=pool)
    sampler.run_mcmc(initial_state=p0, nsteps=nsteps, progress=True)
pkl.dump(sampler, open(out_dir/"mcmc_res.pkl", 'wb'))

chains = sampler.get_chain(flat=False)  # (nsteps, nwalkers, ndim)
labels = ['Gb', 'sigma', 'gamma', 'a_meal', 'log_d_meal', 'a_ins', 'log_d_ins', "beta"]

fig, ax = plt.subplots(chains.shape[2], 1, figsize=(10, 20))
for i in range(chains.shape[2]):
    lab = labels[i]
    for j in range(chains.shape[1]):
        Z = chains[:, j, i].copy()
        if labels[i] == "log_d_meal":
            Z = chains[:, j, i-1] +  np.exp(Z)
            lab = "b_meal"
        elif labels[i] == "log_d_ins":
            Z = chains[:, j, i-1] + np.exp(Z)
            lab = "b_ins"
        ax[i].plot(Z, alpha=0.4, linewidth=0.5)
    ax[i].set_ylabel(lab)
plt.tight_layout()
plt.show()

ndiscard = 200
samples = sampler.get_chain(discard = ndiscard, flat = True) # ((nsteps*nwalkers), ndim)

# \\\\
# \\\\
# –––– Plot
# \\\\
# \\\\
fig, ax = plt.subplots(int(samples.shape[1]/2), 2, figsize=(12, 10))
axes = ax.flatten()
coefs = []
for i in range(samples.shape[1]):
    
    Z = samples[:, i].copy()
    lab = labels[i]

    if labels[i] == "log_d_meal":
        Z = samples[:, i-1] +  np.exp(Z)
        lab = "b_meal"
    elif labels[i] == "log_d_ins":
        Z = samples[:, i-1] + np.exp(Z)
        lab = "b_ins"

    x = np.linspace(Z.min(), Z.max(), 500)
    kde = scipy.stats.gaussian_kde(Z)
    y = kde(x)
    y = y/y.max()
    map = x[np.argmax(y)]
    med = np.median(Z)
    P = .95
    lwr_p, uppr_p = (1 - P)/2, P + (1 - P)/2
    lwr, uppr = np.quantile(Z, q = [lwr_p, uppr_p])
    coefs.append({"par": lab, "median": med,"MAP":map, f"CI_{lwr_p:.3f}": lwr, f"CI_{uppr_p:.3f}":uppr})

    counts, _, _ = axes[i].hist(Z, bins=50, alpha=0.5, color = "black")
    mask = (x >= lwr) & (x <= uppr)
    axes[i].plot(x, y*counts.max())
    axes[i].fill_between(x[mask], y[mask]*counts.max(), alpha=0.3, color="blue")
    axes[i].set_xlabel(lab)
    axes[i].axvline(med, linestyle = '--', color = 'green', label = "Median")
    axes[i].axvline(map, linestyle = '-',linewidth = 5, color = 'orange', label = 'MAP')
    axes[i].legend()
plt.tight_layout()
plt.show()

# \\\\
# \\\\
# –––– Compare Estimates 
# \\\\
# \\\\

MCMC_summ = pd.DataFrame(coefs)
map_params = dict(zip(MCMC_summ.par, MCMC_summ["MAP"]))
# map_params["b_meal"] = map_params["a_meal"] + np.exp(map_params["log_d_meal"])
# map_params["b_ins"] = map_params["a_ins"] + np.exp(map_params["log_d_ins"])
# del map_params["log_d_ins"]; del map_params["log_d_meal"]

actual_pars = {x:sim.__dict__[x] for i, x in enumerate(sim.__dict__) if i <10}


for actual, pred  in zip(actual_pars.items(), map_params.items()):
    k_actual, v_actual = actual
    k_pred, v_pred = pred
    diff = v_actual - v_pred
    print(f"{k_actual}: {diff:.2e}, {100 * diff/ v_actual:.2f}% deviation")

# \\\\
# \\\\
# –––– Simulate again with estimated params
# \\\\
# \\\\
Gb, sigma, gamma, a_meal, a_ins, beta, b_meal, b_ins = map_params.values()

new_sim = GlucoseSim(Gb = Gb, sigma = sigma, gamma = gamma,
                     a_meal = a_meal, b_meal = b_meal,
                     a_ins =  a_ins, b_ins = b_ins, beta = beta,
                     Gmeal = Gmeal, Ibolus = Ibolus , t  = t)
new_sim.run(num_sim = 5)

fig, ax = plt.subplots(figsize = (10,8))
ax.plot(t, G_obs, color = "orange", label = "Actual")
new_sim.plot(ax = ax)
ax.set_ylabel("Blood Glucose")
ax.set_xlabel("Time")
plt.title("Actual Versus Predicted Simulated Glucose Evolution")

# RMSE
np.sqrt(((sim.results[inds,:].mean(axis = 1) - new_sim.results.mean(axis = 1))**2).mean())
# 3.0784249149024587

# MAD
resid  = sim.results[inds,:].mean(axis = 1) - new_sim.results.mean(axis = 1)
np.abs(resid).mean()

plt.plot(range(len(resid)),resid)
plt.axhline(0, color = "black")

x = sim.results[inds,:].mean(axis = 1)
y = new_sim.results.mean(axis = 1)

def lm(x,y):
    MM = np.column_stack([np.ones(len(x)), x])
    return np.linalg.inv(MM.T@MM)@MM.T@y

intrcpt,slope = lm(x,y)
y_hat = intrcpt + slope * x

fig, ax = plt.subplots(figsize = (10,10))
ax.axline(xy1 = (140,140),slope =1, color = "black", linestyle = '--', alpha = .5)
ax.scatter(x,y)
ax.plot(x,y_hat.ravel(), color = "orange")
ax.set_xlabel("Actual")
ax.set_ylabel("Predicted")