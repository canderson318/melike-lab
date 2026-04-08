%%
clc;
clear;


set(0, 'DefaultFigureColor', 'white')
set(0, 'DefaultAxesColor', 'white')


addpath(genpath("/Users/canderson/Documents/school/local-melike-lab/melike-lab/practice-scripts/"))



%% read Data
dat = readtable("/Users/canderson/odrive/home/melike-rotation/project001/data/patXXX/bg.csv");
head(dat)

%% Clean Data
dat = sortrows(dat, "time");
[~,idx] = unique(dat,  "first");
dat = dat(idx, :);

bg = dat.bg;
t = dat.time;

bg = bg(1:1000);
t = t(1:1000);

%% Model Baseline and Reversion Strength
% dG/dt =-𝛉(G-µ)dt
% integrate wrt t
% G(t) = µ + (G_0 - µ)*e^(-𝛉t)  
% get next from current
% G_(i+1) = µ + (G_i - µ)*e^(-𝛉∆t)
% factor
% G_(i+1) = G_i*e^(-𝛉∆t) + µ(1-e^(-𝛉∆t))
% let b = e^(-𝛉∆t), then let a = µ(1-b)
% G_(i+1) = a + b*G_i
% Use linear regression to solve for a and b
% G_(i+1) = a+b*G_i + ε, where ε ~ ꫜ(0, σ^2)

bg_curr = bg(1:end-1);
bg_next = bg(2:end);

%  make training data
p_train = .7;
n_train = floor(length(t) * p_train);

bg_train = bg(1:n_train);
t_train = t(1:n_train);

bg_train_curr = bg_train(1:end-1);
bg_train_next = bg_train(2:end);

% least squares approx. of a and b
model_matrix = [ ones(length(bg_train_curr),1), bg_train_curr];

coeffs = model_matrix \ bg_train_next;
a = coeffs(1);
b = coeffs(2);

% calculate 𝛉 from b
% b = e^(-𝛉∆t)
% 𝛉 = -ln|b| / ∆t
delta_t = mean(diff(t));
theta = -log(abs(b)) / delta_t;

% calculate µ from a
% a = µ(1-b)
% µ = a / (1-b)
mu = a / (1-b);


% G(t) = µ + (G_0 - µ)*e^(-𝛉t)
% G_(i+1) = G_i*e^(-𝛉∆t) + µ(1-e^(-𝛉∆t))
bg_train_next_pred = bg_train_curr * exp(-theta*delta_t) + mu* (1 - exp(-theta*delta_t));

resid = bg_train_next_pred - bg_train(2:end); 
sigma = std(resid);

bg_next_pred = bg_curr * exp(-theta*delta_t) + mu* (1 - exp(-theta*delta_t));

plot(t, bg, "g-", 'DisplayName', "Actual", "LineWidth", 5)
hold on
plot(t(2:end), bg_next_pred, 'b--', 'DisplayName', "Predicted Testing",'LineWidth',3)
plot(t_train(2:end), bg_train_next_pred, 'r--', 'DisplayName', "Predicted Training",'LineWidth',3)
legend()
hold off

%% ODE
% dG/dt =-𝛉(G-µ)dt
ode_fun = @(t,G) -theta * ( G - mu);
G0 = bg(1);  % start from first observed value
[t_pred, G_pred] = ode45(ode_fun, t, G0);

plot(t, bg, "g-", 'DisplayName', "Actual", "LineWidth", 3)
hold on
plot(t_pred, G_pred, "r--", "DisplayName", "Predicted", "LineWidth", 2)
title("Simple ODE")
subtitle("dG/dt = -𝛉(G-µ)dt")
legend()
hold off

%% SDE

nsims = 100;
sims = zeros(length(t), nsims);

for i = 1:nsims
    sims(:,i) = SDE(theta, sigma, mu, bg(1), t);
end

sds = sqrt(var(sims, 0, 2));
mu  = mean(sims, 2);

n = size(sims, 2);
ci = 1.96 * sds / sqrt(n);

fig = figure;

plot(t, mu, 'r-', 'DisplayName', 'Average Prediction')
hold on
plot(t, bg, 'g-', 'LineWidth', 2, 'DisplayName', 'True')
fill([t; flipud(t)], [mu+ci; flipud(mu-ci)], 'b', ...
    'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '95% CI')
title("SDE")
subtitle("dG/dt = -𝛉(G-µ)dt + dW/dt")
legend()
hold off

save_dir = '/Users/canderson/Documents/school/local-melike-lab/melike-lab/practice-scripts/bg_dynamics_out/';

if ~exist(save_dir)
    mkdir(save_dir);
end

figurename_png = '01.png';

print(fig, fullfile(save_dir, figurename_png), '-dpng', '-r300')

%% Function Definitions
function G_SDE = SDE(theta, sigma, mu, y0, t)
n = length(t);
dt = mean(diff(t));
G_SDE = zeros(n,1);
G_SDE(1) = y0;  % semicolon here
for i = 1:n-1
    dW = sqrt(dt) * randn(1);
    G_SDE(i+1) = G_SDE(i) + theta * (mu - G_SDE(i)) * dt + sigma * dW;
end
end