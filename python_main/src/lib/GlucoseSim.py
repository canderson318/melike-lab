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

    def __init__(self, Gb, sigma, gamma, a_meal, b_meal,beta, a_ins, b_ins, t, Gmeal, Ibolus, x = 5, z = 300):
        self.Gb = Gb
        self.sigma = sigma
        self.gamma = gamma
        self.a_meal = a_meal
        self.b_meal = b_meal
        self.a_ins = a_ins
        self.b_ins = b_ins
        self.beta = beta
        self.x= x
        self.z = z
        self.t = t
        self.Gmeal = Gmeal
        self.Ibolus = Ibolus
        self.K = len(t)
        self.results = None  # this is set after calling .run()

    # @staticmethod
    # def calculate_c(a, b, t, tk_m):
    #     dt = t[1] - t[0]
    #     t_minus_tk_m = (t - tk_m)
    #     res = np.sum(np.exp(-a*t_minus_tk_m) - np.exp(-b*t_minus_tk_m)) * dt
    #     return res

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
        # const = GlucoseSim.calculate_const(a, b, t_grid, tj_m)
        const = (a*b)/(b-a)

        # dt
        hk = t_k1 - t_k
        term1 = (np.exp(-a*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-a*(t_k - tj_m))) / (gamma - a)
        term2 = (np.exp(-b*(t_k1 - tj_m)) - np.exp(-gamma*hk) * np.exp(-b*(t_k - tj_m))) / (gamma - b)
        return (Gj * const) * (term1 - term2)

    @staticmethod
    def bolusInsulin(Ij, a, b, tj_ins, gamma, t_k, t_k1, x = 5, z = 300):
        """
        Bolus insulin contribution for a single bolus event.
        Ij    : insulin dose
        c, d  : double-exponential shape params (c < d)
        tj_ins: bolus event time
        x     : onset delay (min) — kernel active starting at tj_ins + x
        z     : end of effect (min) — kernel inactive after tj_ins + z
        gamma : BG decay rate (1/min)
        t_k, t_k1: current and next event times
        """
        onset = tj_ins + x
        end   = tj_ins + z
        hk    = t_k1 - t_k

        # before onset or after end of effect: no contribution
        if not (onset <= t_k and t_k1 <= end):
            return 0.0

        # normalizing const: 1 / ∫_0^{z-x} [exp(-a·τ) - exp(-b·τ)] dτ
        const = 1.0 / ((np.exp(-b*(z-x)) - 1)/b - (np.exp(-a*(z-x)) - 1)/a)

        term1 = (np.exp(-a*(t_k1 - onset)) - np.exp(-gamma*hk) * np.exp(-a*(t_k - onset))) / (gamma - a)
        term2 = (np.exp(-b*(t_k1 - onset)) - np.exp(-gamma*hk) * np.exp(-b*(t_k - onset))) / (gamma - b)
        return const * Ij * (term1 - term2)

    @staticmethod
    def integrate_meal(k, Gmeal, a_meal, b_meal, t, gamma):
        return sum(
            GlucoseSim.meal(Gmeal[j], a_meal, b_meal, t[j], gamma, t[k], t[k+1])
            for j in range(k+1) if Gmeal[j] > 0
        )

    @staticmethod
    def integrate_bolusInsulin(k, Ibolus, a_ins, b_ins, t, gamma, x  = 5, z = 300):
        return sum(
            GlucoseSim.bolusInsulin(Ibolus[j], a_ins, b_ins, t[j], gamma, t[k], t[k+1], x, z)
            for j in range(k+1) if Ibolus[j] > 0
        )

    @staticmethod
    def nextG_raw(Gb, gamma, hk, gk, mk, beta, ik):
        return Gb + np.exp(-gamma * hk) * (gk - Gb) + mk - beta*ik

    @staticmethod
    def fluctuation(sigma, gamma, hk):
        return sigma * np.sqrt(1-np.exp(-2*gamma*hk))

    @staticmethod
    def nextGk(Gb,hk, gamma, gk, mk,ik,beta, sigma):
        """
        Evolve current glucose gk to next as a function of reversion to baseline, meal, and bolus insulin.

        Gb: Basal Glucose.
        hk: time step from current to next timepoint.
        gamma: BG decay rate (1/min)
        gk: current glucose (mg/dl)
        mk: integral of meal function up to next timepoint (mg/dl)
        ik: integral of insulin function up to next timepoint (mg/dl)
        beta: scales effect of insulin on blood glucose (mg/(dl·U))
        sigma: average of random fluctuations
        """
        XIk =  np.random.normal(0,1)
        # evolution = Gb + np.exp(-gamma * hk) * (gk - Gb) + mk - beta*ik
        nextG = GlucoseSim.nextG_raw(Gb, gamma, hk, gk, mk, beta, ik)
        noise = GlucoseSim.fluctuation(sigma, gamma, hk) * XIk
        return nextG + noise

    def run(self, num_sim = 1):
        res = []
        for _ in range(num_sim):
            G = np.zeros(self.K, dtype=float)
            G[0] = self.Gb
            for k in range(self.K - 1):
                # dt
                hk = self.t[k+1] - self.t[k]

                mk_integral = GlucoseSim.integrate_meal(k, self.Gmeal, self.a_meal, self.b_meal, self.t, self.gamma)
                ik_integral = GlucoseSim.integrate_bolusInsulin(k, self.Ibolus, self.a_ins, self.b_ins, self.t, self.gamma, self.x, self.z)

                # evolve current glucose as function meal and params
                G[k+1] = self.nextGk(Gb=self.Gb, hk=hk, gamma=self.gamma, gk=G[k], mk=mk_integral, ik = ik_integral, sigma=self.sigma, beta = self.beta)

            res.append(G.copy())
        self.results = np.array(res).T  # shape (K, num_sim)
        return self

    def plot(self, ax=None, **kwargs):
        if self.results is None:
            raise RuntimeError("Call .run() before .plot()")
        df = pd.DataFrame(self.results)
        mean_G = df.mean(axis=1)
        std_G  = df.std(axis=1)
        if ax is None:
            fig, ax = plt.subplots()
        color = kwargs.pop('color', 'steelblue')
        ax.plot(self.t,mean_G, label='mean', color=color, **kwargs)
        ax.fill_between(self.t, mean_G - std_G, mean_G + std_G,
                        alpha=0.3, color=color, label='±1 SD')
        ax.legend(loc = "upper left", framealpha = 1,bbox_to_anchor=(1.15, 1))
        return ax
