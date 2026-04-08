import numpy as np
import matplotlib.pyplot as plt


n = 1000
t = np.linspace(0,10, n)


# f = lambda t,w: np.cos(w*t)
# x = f(t, 2)
# f = lambda t: (t**2)/ np.sqrt(t)
# x = f(t)
f = lambda t: (t**2)
x = f(t)


plt.axhline(y =0, color  = "grey", linestyle = "--")
plt.plot(t, x, label = "x(t)")

dt = t[1] - t[0]
dxdt = np.zeros(n)
dxdt[1:] = np.diff(x) / dt

plt.plot(t, dxdt, label = "dxdt")


ddxdt= np.zeros(n)
ddxdt[1:] = np.diff(dxdt) / dt


plt.plot(t, ddxdt, label = "ddxdt")
plt.legend()
plt.show()
