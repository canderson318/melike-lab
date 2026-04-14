import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


class GlucoseSim:
    """
    Simulate blood glucose (BG) trajectories under the T2DM event-time model
    from Sirlanci et al. 2022 (eq. 8-9).

    The model is a discrete-time Ornstein-Uhlenbeck process:
        g_{k+1} = Gb + exp(-gamma*h_k)*(g_k - Gb) + m_k + sigma*sqrt(1 - exp(-2*gamma*h_k))*XI_k

    where m_k is the cumulative meal contribution over [t_k, t_{k+1}] and xi_k ~ N(0,1).

    Parameters
    ----------
    Gb    : basal glucose (mg/dl) — mean of the unforced process
    sigma : amplitude of BG oscillations (mg/dl)
    gamma : BG decay rate toward Gb (1/min)
    a, b  : meal absorption shape params (1/min); control rise/fall rate of glucose uptake; a < b
    t     : array of event times
    Gmeal : array of meal glucose loads at each event time (mg/dl)
    """

    def __init__(self, Gb, sigma, gamma, a, b, t, Gmeal):
        self.Gb = Gb
        self.sigma = sigma
        self.gamma = gamma
        self.a = a
        self.b = b
        self.t = t
        self.Gmeal = Gmeal
        self.K = len(t)
        self.results = None  # this is set after calling .run()
    
    @staticmethod
    def calculate_c(a, b, t, tk_m):
        dt = t[1] - t[0]
        t_minus_tk_m = (t - tk_m)
        res = np.sum(np.exp(-a*t_minus_tk_m) - np.exp(-b*t_minus_tk_m)) * dt
        return res
    
    @staticmethod
    def meal(Gj, a, b, tj_m, gamma, t_k, t_k1):
        """
        Meal contribution for a single meal (eq. 9).
        Gj: glucose load of meal (mg/dl), tj_m: meal time,
        t_k/t_k1: current and next event times.
        """
        # # grid of time to integrate a,b over to get c
        # #++ upper bound of integral at t+10/a not inf, 10/a close enough to 0
        # t_grid = np.linspace(tj_m, tj_m + 10/a, 10000)
        # c = GlucoseSim.calculate_c(a, b, t_grid, tj_m)
        c = (a*b)/(b-a) 
        
        # dt
        hk = t_k1 - t_k
        term1 = (np.exp(-a*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-a*(t_k - tj_m))) / (gamma - a)
        term2 = (np.exp(-b*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-b*(t_k - tj_m))) / (gamma - b)
        return (Gj / c) * (term1 - term2)
    
    @staticmethod
    def bolusInsulin(Gj, d, tj_m, gamma, t_k, t_k1):
        """
        """
        # # grid of time to integrate a,b over to get c
        # # ++ upper bound of integral at t+10/a not inf, 10/a close enough to 0
        # t_grid = np.linspace(tj_m, tj_m + 10/a, 10000)
        # c = GlucoseSim.calculate_c(a, b, t_grid, tj_m)
        c = 1/((np.exp(-d*(z-x))-1)/d - (np.exp(-c*(z-x))-1)/c);
        
        # dt
        hk = t_k1 - t_k
        term1 = (np.exp(-a*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-a*(t_k - tj_m))) / (gamma - a)
        term2 = (np.exp(-b*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-b*(t_k - tj_m))) / (gamma - b)
        return (Gj / c) * (term1 - term2)
    
    @staticmethod
    def nextGk(Gb,hk, gamma, gk, mk, sigma):
        XIk =  np.random.normal(0,1)
        evolution = Gb + np.exp(-gamma * hk) * (gk - Gb) + mk 
        noise = sigma * np.sqrt(1-np.exp(-2*gamma*hk))* XIk
        return evolution + noise

    def run(self, num_sim):
        res = []
        for _ in range(num_sim):
            G = np.zeros(self.K, dtype=float)
            G[0] = self.Gb
            for k in range(self.K - 1):
                # dt
                hk = self.t[k+1] - self.t[k]
                # meal function: integrate over meals up to current time
                mk = sum(
                    self.meal(self.Gmeal[j], self.a, self.b, self.t[j], self.gamma, self.t[k], self.t[k+1])
                    for j in range(k+1) if self.Gmeal[j] > 0
                    )
                # evolve current glucose as function meal and params
                G[k+1] = self.nextGk(Gb=self.Gb, hk=hk, gamma=self.gamma, gk=G[k], mk=mk, sigma=self.sigma)
            res.append(G.copy())
        self.results = np.array(res).T  # shape (K, num_sim)
        return self

    def plot(self, ax=None):
        if self.results is None:
            raise RuntimeError("Call .run() before .plot()")
        df = pd.DataFrame(self.results)
        mean_G = df.mean(axis=1)
        std_G  = df.std(axis=1)
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(mean_G, color='steelblue', label='mean')
        ax.fill_between(df.index, mean_G - 1*std_G, mean_G + 1*std_G,
                        alpha=0.3, color='steelblue', label='±1 SD')
        ax.legend(loc = "upper left")
        return ax
    