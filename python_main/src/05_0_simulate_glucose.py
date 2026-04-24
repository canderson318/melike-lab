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
out_dir = local_output_path/ "05_0"
out_dir.mkdir(exist_ok = True)


# \\\\
# \\\\
# –––– Simulate Glucose Evolution over time with different parameters
# \\\\
# \\\\

np.random.seed(123)
# simulate meals
K = 500
blood_vol = 5 * 10

t = np.arange(K)
num_meals = int(10)
meals = np.zeros(K)
meals[np.random.choice(range(K-100), num_meals)] = 2000
Gmeal = np.maximum(meals / blood_vol, .01)

# simulate bolus insulin
Ib = np.random.normal(150,3,K) / blood_vol
Ib = np.maximum(Ib,0)
num_Ibolus = 15
Ib[np.random.choice(range(K), K-num_Ibolus, replace = False)] = 0

# simulate basal insulin
I = np.random.normal(60,1,K) / blood_vol
I = np.maximum(I,0)


# set params
Gb = 140
sigma = 20
gamma = .035
a_meal = .01
b_meal = .03
a_ins = .05
b_ins = .06
beta = 30

# run simulation
sim = GlucoseSim(t = t, Gb=Gb, sigma=sigma, gamma=gamma, 
                 a_meal=a_meal, b_meal= b_meal, 
                 Gmeal=Gmeal, 
                 Ibasal = I, Ibolus = Ib, a_ins = a_ins, b_ins = b_ins, beta = beta)

sim.run(num_sim=10)

# sim.plot()

# \\\
# ––– Plot 
# \\\
params = {x[0]: x[1] for x in sim.__dict__.items() if x[0] not in ['t', 'Gmeal', "Ibolus", 'results', "Ibasal"]}
param_str = ", ".join([f"{x[0]}: {x[1]}" for x in params.items()])
param_str = param_str.replace("b_meal", "\nb_meal")

fig, ax = plt.subplots(4,1,figsize = (15,10))
# ax[0].plot(t, sim.results)
sim.plot(ax = ax[0])
ax[0].set_ylabel("Glucose", color="steelblue")
ax[0].tick_params(axis='y', labelcolor="steelblue")
ax[0].legend(bbox_to_anchor=(1.05, 1))

ax[1].bar(t, Gmeal, color='orange', width=1, label='Carbs',alpha = 1)
ax[1].set_ylabel('Carbs', color="orange")
ax[1].tick_params(axis='y', labelcolor='orange')
ax[1].set_xlabel("")

ax[2].plot(t, I, color='green', label='Basal', alpha = .3)
ax[2].set_ylabel('Basal Insulin', color="green")
ax[2].tick_params(axis = "y", labelcolor = "green")
ax[2].legend(bbox_to_anchor=(1.05, 1))

ax[3].bar(t, Ib, color='green', label='Bolus', alpha = .7)
ax[3].set_ylabel('Bolus Insulin', color="green")
ax[3].tick_params(axis = "y", labelcolor = "green")
ax[3].legend(bbox_to_anchor=(1.05, 1))

# [ax[i].xaxis.set_major_locator(plt.MultipleLocator(.5)) for i in range(3)]
plt.suptitle(f"Simulated Glucose evolution", y=1, fontsize=14)
plt.tight_layout()
plt.xlabel("Time")

param_path = "_".join([f"{x[0]}{x[1]}" for x in params.items()])
plt.savefig(out_dir/f"simple_sim_{param_path}.pdf")

with  open(out_dir/"sim.pkl", 'wb') as f:
    pkl.dump(sim,f)
    

#@\\\
#@\\\
      # Play around with params
#@\\\
#@\\\
      
