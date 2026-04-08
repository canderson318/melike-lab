# OU process simulation
# credit: https://www.quantstart.com/articles/ornstein-uhlenbeck-simulation-with-python/
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression as lm

np.random.seed(1103)
# Parameters for the OU process
theta = 0.03      # Speed of mean reversion
mu = 0     # Long-term mean
sigma = 10      # Volatility
X0 = 1.0         # Initial value
T = 500         # Total time
dt = 20       # Time step
N = int(T / dt)  # Number of time steps

# Pre-allocate array for efficiency
X = np.zeros(N)
X[0] = X0

# Generate the OU process
for t in range(1, N):
    dW = np.sqrt(dt) * np.random.normal(0, 1)
    X[t] = X[t-1] + theta * (mu - X[t-1]) * dt + sigma * dW

# Plot the result
plt.plot(np.linspace(0, T, N), X, color = "orange")
plt.title("")
plt.xlabel("Time")
plt.ylabel("X(t)")
plt.show()



# estimate OU params
diff = np.diff(X)
_X = X[:-1].reshape(-1,1)

reg  = lm(fit_intercept=True)
reg.fit(_X, diff)

theta_hat = - reg.coef_[0] / dt 
mu_hat = reg.intercept_ / theta_hat

y_hat = reg.predict(_X)
sigma_hat = np.std(diff - y_hat) / np.sqrt(dt)

rmse = lambda y, yhat: np.sqrt(np.mean((y-yhat)**2))


# reintegrate to function
X_hat = np.zeros(N)
X_hat[0] = X[0] 

for t in range(1, N):
    # Predict the next difference based on the PREDICTED value
    X_hat[t] = X_hat[t-1] + reg.predict(X_hat[t-1].reshape(-1,1))


print(f"**Parameters**\ntheta = {theta}\nmu = {mu}\nsigma = {sigma}\n")
print(f"**Estimated Parameters**\ntheta = {theta_hat}\nmu = {mu_hat}\nsigma = {sigma_hat}\nRMSE = {rmse(diff,y_hat)}")
print(f"**Integraged**\nRMSE = {rmse(X, X_hat)}")

# plot
# plot
# plot
# plot
# plot



fig, ax = plt.subplots(2,1, figsize = (10,5))
time = np.linspace(0,T,N)
ax[0].plot(time, X, color="orange", label="Actual Function")
ax[0].plot(time, X_hat, color="blue", label="Fitted Function")
ax[0].set_title("Ornstein-Uhlenbeck Process Simulation")
ax[0].set_xlabel("Time")
ax[0].set_ylabel("X(t)")
ax[0].legend()

time = np.arange(len(_X))
ax[1].plot(time, diff, color="orange", label="Actual diff")
ax[1].plot(time, y_hat, color="blue", label="Fitted diff")
ax[1].set_xlabel("Time")
ax[1].set_ylabel("dX(t)")
ax[1].legend()
plt.show()

