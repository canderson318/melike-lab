% =========================================================================
% ODE Modeling Example: Logistic Growth with Noisy Data
% =========================================================================
% This example demonstrates:
% 1. Simulating data from an ODE (logistic growth model)
% 2. Adding realistic noise
% 3. Fitting the ODE to recover parameters
% 4. Comparing ODE to simple linear regression
% =========================================================================

clear ; close all; clc;

%% STEP 1: Define the TRUE ODE model
% We'll use logistic growth: dy/dt = r*y*(1 - y/K)
% where:
%   y = population/disease severity/metabolite level
%   r = growth rate
%   K = carrying capacity (maximum value)

% True parameters (what we'll try to recover)
r_true = 0.5;      % growth rate
K_true = 100;      % carrying capacity
y0_true = 5;       % initial condition

% Define ODE function
logistic_ode = @(t, y, r, K) r * y * (1 - y/K);

%% STEP 2: Generate "true" data by solving the ODE
% Time points where we'll "observe" data
t_obs = linspace(0, 15, 20)';  % 20 observations over 15 time units

% Solve ODE with true parameters
[t_true, y_true] = ode45(@(t,y) logistic_ode(t, y, r_true, K_true), t_obs, y0_true);

%% STEP 3: Add realistic noise to create "observed" data
% Add Gaussian noise (like measurement error)
noise_level = 10;  % standard deviation of noise
y_obs = y_true + noise_level * randn(size(y_true));

% Ensure no negative values (realistic for many biological measurements)
y_obs = max(y_obs, 0.1);

fprintf('=== TRUE PARAMETERS ===\n');
fprintf('Growth rate (r): %.3f\n', r_true);
fprintf('Carrying capacity (K): %.3f\n', K_true);
fprintf('Initial value (y0): %.3f\n\n', y0_true);

%% STEP 4: Naive approach - Linear Regression
% Fit: y = beta0 + beta1*t
X = [ones(length(t_obs), 1), t_obs];
beta_linear = X \ y_obs;  % least squares

% Predict with linear model
t_fine = linspace(0, 20, 200)';
y_pred_linear = beta_linear(1) + beta_linear(2) * t_fine;

fprintf('=== LINEAR REGRESSION RESULTS ===\n');
fprintf('Intercept: %.3f\n', beta_linear(1));
fprintf('Slope: %.3f\n', beta_linear(2));
fprintf('(Notice: predicts y(20) = %.1f - could go negative or exceed limits!)\n\n', ...
        beta_linear(1) + beta_linear(2)*20);

%% =========================================================================
% Objective function for ODE fitting
% =========================================================================
function error = ode_objective(params, t_obs, y_obs)
    % Extract parameters
    r = params(1);
    K = params(2);
    y0 = params(3);
    
    % Check for invalid parameters
    if r <= 0 || K <= 0 || y0 <= 0
        error = 1e10;  % penalty for invalid parameters
        return;
    end
    
    % Solve ODE with current parameters
    try
        [~, y_pred] = ode45(@(t,y) r * y * (1 - y/K), t_obs, y0);
        
        % Calculate sum of squared errors
        error = sum((y_obs - y_pred).^2);
    catch
        % If ODE solver fails, return large error
        error = 1e10;
    end
end



%% STEP 5: Fit ODE model to recover parameters
% We need to estimate r, K, and y0 from the noisy data

% Define objective function: sum of squared errors
objective = @(params) ode_objective(params, t_obs, y_obs);

% Initial guess for parameters (deliberately wrong to test recovery)
params_init = [0.3, 80, 3];  % [r, K, y0]

% Optimize to find best parameters
options = optimset('Display', 'iter', 'MaxFunEvals', 10e3);
params_fit = fminsearch(objective, params_init, options);

r_fit = params_fit(1);
K_fit = params_fit(2);
y0_fit = params_fit(3);

fprintf('\n=== ODE FIT RESULTS ===\n');
fprintf('Estimated growth rate (r): %.3f (true: %.3f)\n', r_fit, r_true);
fprintf('Estimated capacity (K): %.3f (true: %.3f)\n', K_fit, K_true);
fprintf('Estimated initial (y0): %.3f (true: %.3f)\n\n', y0_fit, y0_true);

% Solve ODE with fitted parameters for visualization
[~, y_fit] = ode45(@(t,y) logistic_ode(t, y, r_fit, K_fit), ...
                   t_fine, y0_fit);

%% STEP 6: Visualize everything
figure('Position', [100 100 1200 500]);

% Panel 1: Data and both fits
subplot(2,1,1);
hold on;
plot(t_obs, y_obs, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Observed data (noisy)');
plot(t_true, y_true, 'b-', 'LineWidth', 2, 'DisplayName', 'True ODE trajectory');
plot(t_fine, y_fit, 'r--', 'LineWidth', 2,'DisplayName', 'Fitted ODE');
plot(t_fine, y_pred_linear, 'g:', 'LineWidth', 2, 'DisplayName', 'Linear regression');
hold off;
xlabel('Time', 'FontSize', 12);
ylabel('Value (y)', 'FontSize', 12);
title('ODE vs Linear Regression Fit', 'FontSize', 14);
legend('Location', 'southeast');
grid on;
xlim([0 20]);
ylim([0 120]);

% Panel 2: Residuals comparison
subplot(2,1,2);
% Calculate predictions at observation times for both models
[~, y_ode_at_obs] = ode45(@(t,y) logistic_ode(t, y, r_fit, K_fit), ...
                          t_obs, y0_fit);
y_linear_at_obs = beta_linear(1) + beta_linear(2) * t_obs;

residuals_ode = y_obs - y_ode_at_obs;
residuals_linear = y_obs - y_linear_at_obs;

hold on;
plot(t_obs, residuals_ode, 'ro-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('ODE (RMSE=%.2f)', rms(residuals_ode)));
plot(t_obs, residuals_linear, 'go-', 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Linear (RMSE=%.2f)', rms(residuals_linear)));
plot([0 15], [0 0], 'w--', 'LineWidth', 1);
hold off;
xlabel('Time', 'FontSize', 12);
ylabel('Residuals (observed - predicted)', 'FontSize', 12);
title('Model Residuals Comparison', 'FontSize', 14);
legend('Location', 'best');
grid on;

%% STEP 7: Show the danger of extrapolation
fprintf('=== EXTRAPOLATION COMPARISON (t=30) ===\n');
[~, y_ode_extrap] = ode45(@(t,y) logistic_ode(t, y, r_fit, K_fit), ...
                          [0 30], y0_fit);
y_linear_extrap = beta_linear(1) + beta_linear(2) * 30;

fprintf('ODE prediction at t=30: %.2f (plateaus at K=%.1f)\n', ...
        y_ode_extrap(end), K_fit);
fprintf('Linear prediction at t=30: %.2f (keeps growing!)\n', ...
        y_linear_extrap);
fprintf('\nThe ODE correctly captures saturation, while linear model fails.\n');

