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
output_real_path = '/Users/canderson/odrive/home/melike-rotation/project001/outputs/'
os.listdir(output_real_path)
local_output_path = WD/'outputs'

sym_link(output_real_path, local_output_path)

# set output directories
out_dir = local_output_path/ "06"
out_dir.mkdir(exist_ok = True)

# load simulation
sim = pkl.load(open(local_output_path/"05"/"sim.pkl", 'rb'))

# \\\\
# \\\\
# –––– Estimate Parameters From Simulation
# \\\\
# \\\\

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
    # integrate meal function
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

# extract simulated patient
t, G_obs, Gmeal, Ibolus = sim.t, sim.results.mean(axis = 1), sim.Gmeal, sim.Ibolus

# subset 
# N = 100
# inds = np.sort(np.random.choice(range(len(t)), N, replace = False))
inds = range(0,100)
t, G_obs, Gmeal, Ibolus  = t[inds], G_obs[inds], Gmeal[inds], Ibolus[inds] 


# # optimize b as function a
# a_meal_actual, b_meal_actual, a_ins_actual, b_ins_actual = [sim.__dict__[x] for x in ['a_meal', 'b_meal', 'a_ins', 'b_ins']]
# # b_meal = a_meal + np.exp(log_d_meal)
# log_d_meal_actual = np.log(b_meal_actual - a_meal_actual)
# log_d_ins_actual = np.log(b_ins_actual - a_ins_actual)
# np.exp(log_d_meal_actual) + a_meal_actual

#\\\\
#\\\\
# ––– Optimize
#\\\\
#\\\\

plt.plot(t, G_obs)

#                  Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta
# start_points = [140,     5,   .06,     .1,         -1,    .1,        -1,   50],
start_points   = [150,     10,   .07,    .04,         -4,    .01,        -2,   100]

#                Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta
bounds =  [(100,200),(5,15),(.06,.08),(.01,.07),(-100,5), (.001,.02),   (-100,5), (50,200)]

result = optimize.minimize(
    neg_log_likelihood,
    x0=start_points,
    args=(t, G_obs, Gmeal, Ibolus),
    bounds=bounds,
    method='Nelder-Mead', options = {'maxiter':2_000, "adaptive": False}
)

# result = optimize.differential_evolution(
#     neg_log_likelihood, 
#      x0=start_points,
#     args=(t, G_obs, Gmeal, Ibolus),
#     bounds=bounds, 
#     workers = 1,
#     maxiter=10,
# )

map_params = result.x
map_params


Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta = map_params
b_meal = a_meal + np.exp(log_d_meal)
b_ins = a_ins + np.exp(log_d_ins)

new_sim = GlucoseSim(Gb = Gb, sigma = sigma, gamma = gamma, 
                     a_meal = a_meal, b_meal = b_meal,
                     a_ins =  a_ins, b_ins = b_ins, beta = beta,
                     Gmeal = Gmeal, Ibolus = Ibolus , t  = t)
new_sim.run(num_sim = 5)

fig, ax = plt.subplots(figsize = (10,8))
ax.plot(t, G_obs, color = "orange", label = "Actual")
new_sim.plot(ax = ax)
ax.set_ylabel("Blood Glucose")


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

MM = np.column_stack([np.ones(len(x)), x])

intrcpt,slope = np.linalg.inv(MM.T@MM)@MM.T@y
y_hat = intrcpt + slope * x


# from sklearn.linear_model import LinearRegression
# LR = LinearRegression()
# LR.fit(x.reshape((-1,1)),y.reshape((-1,1)))
# intrcpt,slope = LR.intercept_,LR.coef_
# y_hat = intrcpt + slope * x


fig, ax = plt.subplots(figsize = (10,10))
ax.axline(xy1 = (140,140),slope =1, color = "black", linestyle = '--', alpha = .5)
ax.scatter(x,y)
ax.plot(x,y_hat.ravel())
ax.set_xlabel("Actual")
ax.set_ylabel("Predicted")