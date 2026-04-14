import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import os
from pathlib import Path
import subprocess as sp
from src.utils.tools import *
from src.utils.GlucoseSim import *
import re

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






# \\\\
# \\\\
# –––– Simulate Glucose Evolution ––––
# \\\\
# \\\\
# simulate meals
K = 500
t = np.arange(K)
num_meals = int(10)
meals = np.zeros(K)
meals[np.random.choice(range(K-100), num_meals)] = 500
blood_vol = 5 * 10
Gmeal = np.maximum(meals / blood_vol, .01)

# set params
Gb = 140
sigma = 5
gamma = .06
# a     = .5 + gamma # peaks quickly
# b     = .3 + gamma # drops towards baseline quickly
a = .1
b = .15

# run simulation
sim = GlucoseSim(Gb=Gb, sigma=sigma, gamma=gamma, a=a, b=b, t=t, Gmeal=Gmeal)
sim.run(num_sim=5)

fig, ax1 = plt.subplots()
sim.plot(ax=ax1)
plt.xlabel("Time")
ax1.set_ylabel("Glucose", color = "steelblue")
ax1.tick_params(axis = 'y',  labelcolor = "steelblue")
ax2 = ax1.twinx()
ax2.plot(t, Gmeal, color='orange', alpha=0.6, label='Gmeal')
ax2.set_ylabel('Gmeal (mg/dl)', color = "orange")
ax2.tick_params(axis = 'y',labelcolor = "orange")
ax2.legend(loc='upper right')
plt.title(f"Simulated Glucose evolution\nGb = {Gb}, sigma = {sigma}, gamma = {gamma}, a = {a}, b = {b}")
plt.savefig(out_dir/f"simple_sim_Gb{Gb}_sigma{sigma}_a{a}_b{b}_gamma{gamma}.pdf")


# \\\
# \\\
# ––– Plot effect of gamma on sim
# \\\
# \\\
gammas = np.linspace(0.0001,.06, 20)

fig, ax = plt.subplots(5,4, figsize = (20,10))
axes = ax.flatten()
for i,gamm in enumerate(gammas):
    sim = GlucoseSim(Gb=Gb, sigma=sigma, gamma=gamm, a=a, b=b, t=t, Gmeal=Gmeal)
    sim.run(num_sim=30)
    sim.plot(ax = axes[i])
    axes[i].set_title(f"Gamma = {gamm:.3f}")
plt.tight_layout()
plt.savefig(out_dir/"gamma_space_sims.pdf")
plt.show()
