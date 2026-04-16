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
from scipy.optimize import minimize
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

# optimize b as function a
a_meal_actual, b_meal_actual, a_ins_actual, b_ins_actual = [sim.__dict__[x] for x in ['a_meal', 'b_meal', 'a_ins', 'b_ins']]

# b_meal = a_meal + np.exp(log_d_meal)
log_d_meal_actual = np.log(b_meal_actual - a_meal_actual)
log_d_ins_actual = np.log(b_ins_actual - a_ins_actual)
np.exp(log_d_meal_actual) + a_meal_actual

# neg_log_likelihood((140,5,.06,a_meal_actual, log_d_meal_actual, a_ins_actual, log_d_ins_actual, 50), t, G_obs, Gmeal, Ibolus)

#\\\\
#\\\\
# ––– Optimize 
#\\\\
#\\\\

#                  Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta 
# start_points = [140,     5,   .06,     .1,         -1,    .1,        -1,   50],   
start_points   = [200,     10,   .06,    .01,         -1,    .01,        -1,   100]

#                Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta 
bounds =  [(80,300),(0.1,50),(1e-4,1),(1e-4,1),   (-100,5), (1e-4,1),   (-100,5), (0,500)]

# result = minimize(  
#     neg_log_likelihood,
#     x0=start_points,   
#     args=(t, G_obs, Gmeal, Ibolus), 
#     bounds=bounds,
#     # method='L-BFGS-B',
#     method='Nelder-Mead',
#     options = {'factr':1e5,'maxiter':1000,'maxfun':10_000}
# )
# pkl.dump(result, open (out_dir / "optim_sim_res.pkl", 'wb'))
result = pkl.load(open(out_dir / "optim_sim_res.pkl", 'rb'))

map_params = result.x
map_params


Gb, sigma, gamma, a_meal, log_d_meal, a_ins, log_d_ins, beta = map_params
b_meal = a_meal + np.exp(log_d_meal)
b_ins = a_ins + np.exp(log_d_ins)
                                                                                                                   
new_sim = GlucoseSim(Gb = Gb, sigma = sigma, gamma = gamma, a_meal = a_meal, b_meal = b_meal,a_ins =  a_ins, b_ins = b_ins, beta = beta, 
                     Gmeal = Gmeal, Ibolus = Ibolus , t  = t)


fig, ax = plt.subplots(figsize = (10,8))
ax.plot(t, G_obs, color = "orange", label = "Actual")
new_sim.run(num_sim = 5).plot(ax = ax) 

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

# H = Hessian(lambda p: neg_log_likelihood(p, t, G_obs, Gmeal, Ibolus), 
#             method = "central", order = 2, step = 1e-3)(
#                 map_params
#                 )
# H = pd.DataFrame(H, index = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"], columns = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"])
# pkl.dump(H, open (out_dir / "hessian.pkl", 'wb'))
H = pkl.load(open (out_dir / "hessian.pkl", 'rb'))

# posterior covariance and std errors                                                                    
cov = inv(H)
std_errors = np.sqrt(np.diag(cov))
err  = pd.DataFrame({"SE": std_errors, "MAP":map_params}, index = ["Gb", "sigma", "gamma", "a_meal", "log_d_meal", "a_ins", "log_d_ins", "beta"])
err