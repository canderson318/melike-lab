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
out_dir = local_output_path/ "05"
out_dir.mkdir(exist_ok = True)


# \\\\
# \\\\
# –––– Simulate Glucose Evolution over time with different parameters
# \\\\
# \\\\

np.random.seed(1020)
# simulate meals
K = 500
t = np.arange(K)
num_meals = int(10)
meals = np.zeros(K)
meals[np.random.choice(range(K-100), num_meals)] = 500
blood_vol = 5 * 10
Gmeal = np.maximum(meals / blood_vol, .01)

# simulate bolus insulin
num_Ibolus = 4
I = np.zeros(K)
I[np.random.choice(range(K-100), num_Ibolus)] = 100
Ib = np.maximum(I / blood_vol, .01)


# set params
Gb = 140
sigma = 10
gamma = .08
# a     = .5 + gamma # peaks quickly
# b     = .3 + gamma # drops towards baseline quickly
a_meal = .05
b_meal = .07
a_ins = .01
b_ins = .1
beta = 50

# run simulation
sim = GlucoseSim(t = t, Gb=Gb, sigma=sigma, gamma=gamma, a_meal=a_meal, b_meal= b_meal, Gmeal=Gmeal, Ibolus = Ib, a_ins = a_ins, b_ins = b_ins, beta = beta)
sim.run(num_sim=5)

params = {x[0]: x[1] for x in sim.__dict__.items() if x[0] not in ['t', 'Gmeal', "Ibolus", 'results']}
param_str = ", ".join([f"{x[0]}: {x[1]}" for x in params.items()])
param_str = param_str.replace("b_meal", "\nb_meal")

fig, ax1 = plt.subplots()
sim.plot(ax = ax1)
plt.xlabel("Time")
ax1.set_ylabel("Glucose", color = "steelblue")
ax1.tick_params(axis = 'y',  labelcolor = "steelblue")
ax2 = ax1.twinx()
ax2.plot(t, Gmeal, color='orange', alpha=0.6, label='Gmeal')
ax2.set_ylabel('Gmeal (mg/dl)', color = "orange", labelpad = 40)
ax2.tick_params(axis='y', labelcolor='orange', direction='in', pad=-8)
for label in ax2.get_yticklabels():                                                                  
      label.set_horizontalalignment('right')       
ax2.legend(loc = "upper right", bbox_to_anchor=(1.4, .8))
ax3 = ax1.twinx()
ax3.spines['right'].set_position(('outward', 0))
# ax3.set_ylim((None,None))
ax3.plot(t, Ib, color='green', alpha=0.6, label='Ibolus')
ax3.set_ylabel('Ib (mg/dl)', color = "green", labelpad = 25)
ax3.tick_params(axis = 'y',labelcolor = "green")
ax3.legend(loc="upper right", bbox_to_anchor=(1.4, .7))
plt.suptitle(f"Simulated Glucose evolution", y = 1.02, fontsize = 14)
plt.title(param_str, fontsize = 9)
param_path = "_".join([f"{x[0]}{x[1]}" for x in params.items()])
plt.savefig(out_dir/f"simple_sim_{param_path}.pdf")


with  open(out_dir/"sim.pkl", 'wb') as f:
    pkl.dump(sim,f)
    