function [phi_prop,chainTime,outerTime,outMSE,outRcorr] = f_likelihood_OU_T1D_v2(likelihood_input)

% read input structure
theta_est_names = likelihood_input.theta_est_names;
proposal = likelihood_input.proposal;
preprocessed_data = likelihood_input.preprocessed_data;
%data = likelihood_input.data;
% foo_model_parameters = likelihood_input.foo_model_parameters;
foo_model_parameters = likelihood_input.settingsModel.model_parameters;
%driverConverter = likelihood_input.driverConverter;
dim_meas = likelihood_input.dim_meas;
%settingsSmoother = likelihood_input.settingsSmoother;
% settingsModel = likelihood_input.settingsModel;
% forecast_times = likelihood_input.forecast_times;
% state_inits = likelihood_input.state_inits;
% ODEoptions = likelihood_input.ODEoptions;
% empty_driver_data = likelihood_input.empty_driver_data;
% measuredStates_model = likelihood_input.measuredStates_model;
% measurementConverter = likelihood_input.measurementConverter;
% t_end_driver = likelihood_input.t_end_driver;
outer_timer = likelihood_input.outer_timer;
%sig = likelihood_input.sig;

%% changing values of the parameters to be estimated to the proposal values
for i=1:numel(theta_est_names)
    p_name = char(theta_est_names(i));
    foo_var = proposal(i);
    foo_model_parameters.(p_name) = foo_var;
end

%% extracting parameter values
delta = foo_model_parameters.delta;
% % % % delta = 0; %CHANGE (1,1)-ELEMENT OF S BELOW
c_gamma = foo_model_parameters.gamma;
c_sigma = foo_model_parameters.sigma;
c_beta_d = foo_model_parameters.beta_d;
c_beta_n = foo_model_parameters.beta_n;
Gb = foo_model_parameters.Gb;
T = foo_model_parameters.T;
dL_blood_volume = foo_model_parameters.vg;
g0 = foo_model_parameters.g0;
g0_var = foo_model_parameters.g0_var;
a = foo_model_parameters.a;
b = foo_model_parameters.b;

%% form the big table containing nutrition, insulin, and bg measurements all together

nutrition = preprocessed_data.drivers.nutrition_oral;
basal_insulin = preprocessed_data.drivers.insulin_basal;
bolus_insulin = preprocessed_data.drivers.insulin_bolus;

t = reshape(preprocessed_data.raw_measurements.time, [length(preprocessed_data.raw_measurements.time),1]);
bg = preprocessed_data.raw_measurements.bg;
bg = reshape(bg,[length(bg),1]);
measurement = array2table([t,bg],'VariableNames',{'time','bg'});
measurement.nutrition = zeros(length(t),1);
measurement.basal_insulin = zeros(length(t),1);
measurement.bolus_insulin = zeros(length(t),1);
measurement.meas_time = ones(length(t),1);

meal_time = nutrition(:,1);
nut = table();
nut.time = meal_time;
nut.bg = zeros(length(meal_time),1);
nut.nutrition = nutrition(:,2);
nut.basal_insulin = zeros(length(meal_time),1);
nut.bolus_insulin = zeros(length(meal_time),1);
nut.meas_time = zeros(length(meal_time),1);

basal_ins_time = basal_insulin(:,1);
basal_ins = table();
basal_ins.time = basal_ins_time;
basal_ins.bg = zeros(length(basal_ins_time),1);
basal_ins.nutrition = zeros(length(basal_ins_time),1);
basal_ins.basal_insulin = basal_insulin(:,2);
basal_ins.bolus_insulin = zeros(length(basal_ins_time),1);
basal_ins.meas_time = zeros(length(basal_ins_time),1);

x = 5; %starts acting in 5 min
z = 300; %effective around 5 hours (neglecting the first 5 min)
bolus_ins_time = sort(unique([bolus_insulin(:,1); bolus_insulin(:,1)+x; bolus_insulin(:,1)+z]));
bolus_ins = table();
bolus_ins.time = bolus_ins_time;
bolus_ins.bg = zeros(length(bolus_ins_time),1);
bolus_ins.nutrition = zeros(length(bolus_ins_time),1);
bolus_ins.basal_insulin = zeros(length(bolus_ins_time),1);
bolus_ins.bolus_insulin = zeros(length(bolus_ins_time),1);
bolus_ins.bolus_insulin(ismember(bolus_ins.time,bolus_insulin(:,1))) = bolus_insulin(:,2);
bolus_ins.meas_time = zeros(length(bolus_ins_time),1);

big_nut = [nut;basal_ins;bolus_ins;measurement];
big_nut = sortrows(big_nut,1);

groups = findgroups(big_nut.time);
my_f = @(time,bg,nut,basal_ins,bolus_ins,is_meas) [time(1),sum(bg),max(nut),max(basal_ins),max(bolus_ins),max(is_meas)];
G = splitapply(my_f,big_nut,groups);

G(:,7) = DoubleExponentialMeals4(G(:,1),nutrition,c_gamma,a,b)/dL_blood_volume; %processed_nutrition (might check this but also integrated nutrition)

%% actual timing for daytime and nighttime beta

start_time = likelihood_input.actual_start_time;

% Compute the actual times for each row
real_times = start_time + minutes(G(:,1));

% Extract the hour of day for each timestamp
hour_of_day = hour(real_times) + minute(real_times)/60;

% Logical mask for day vs night period
which_time = likelihood_input.which_time;
is_day = (hour_of_day >= which_time) & (hour_of_day < which_time+12);

hk = G(2:end,1)-G(1:end-1,1);
hk(end+1) = hk(end); %the last one is only a placeholder
y = find(G(:,1) <= basal_insulin(1,1));
b = y(end);
for i = b:size(G,1)
    y = find(basal_insulin(:,1) <= G(i,1));
    ind = y(end);
    G(i,4) = basal_insulin(ind,2)/60; %convert from U/hr to U/min.     
end
integrated_basal_insulin = BasalInsulin(G(:,1),basal_insulin,c_gamma);
G(is_day,8) = c_beta_d * integrated_basal_insulin(is_day); %this is the integrated basal insulin rate
G(~is_day,8) = c_beta_n * integrated_basal_insulin(~is_day); %this is the integrated basal insulin rate

% y = 60; %peaks at 60 min
%Currently the specified peak time by y is not used. The values c = 0.015
%and d = 0.02 results peak value at around 55 min.
integrated_bolus_insulin = DoubleExponentialInsulin2(G(:,1),bolus_insulin,c_gamma,0.015,0.02,x,z);
G(is_day,9) = c_beta_d * integrated_bolus_insulin(is_day);
G(~is_day,9) = c_beta_n * integrated_bolus_insulin(~is_day);

norm_cons = 1./(1-exp(-c_gamma*hk));
G(:,10) = norm_cons.*(G(:,7)-G(:,8)-G(:,9)); %this is n_k = norm_constant*(m_k - i1_k - i2_k)
G(:,10) = Gb+G(:,10);

big_nut = array2table(G,'VariableNames',{'time','bg','nutrition','basal_insulin',...
    'bolus_insulin','meas_time','integrated_nutrition','integrated_basal_insulin',...
    'integrated_bolus_insulin','total_nut_and_ins_rate'});

% %         groups = findgroups(big_nut.time);
% %         my_f = @(time,bg,nut,is_meas) [time(1),sum(bg),sum(nut),max(is_meas)];
% %         G = splitapply(my_f,big_nut,groups);
G = table2array(big_nut);
first_idx = find(G(:,1)==measurement.time(1));
last_idx = find(G(:,1)==measurement.time(end));
G = G(first_idx:last_idx,:);
big_table = array2table(G,'VariableNames',{'time','bg','nutrition_rate',...
    'basal_insulin_rate','bolus_insulin_rate','is_meas','integrated_nutrition',...
    'integrated_basal_insulin','integrated_bolus_insulun','total_nut_and_ins_rate'});

%% check if the bg measurements is nan at first time
% % % % if isnan(big_table.bg(1))
% % % %     g0 = foo_model_parameters.g0;
% % % % else
% % % %     g0 = big_table.bg(1);
% % % % end

chainTimer = tic;

%% constrcting mean
N=size(G,1)-1;

A = -c_gamma * (repmat(G(:,1),1,N+1)-repmat(G(:,1)',N+1,1));
expA = exp(A);

A1 = expA(:,2:end);
A1 = tril(A1,-1);

A2 = expA(:,1:end-1);
A2 = tril(A2,-1);

m_0 = ((A1-A2)*G(1:end-1,10))+(expA(:,1)*g0);

%% constructing covariance matrix

%% the followig part is for piecewise constant sigma
% B1 = sparse(tril(expA(:,2:end),-1));
% C1 = B1;
% 
% B1 = repmat(B1,N+1,1);
% C1 = reshape(C1,(N+1)*N,1);
% C1 = repmat(C1,1,N+1);
% C1 = reshape(C1',(N+1)^2,N);
% 
% B2 = sparse(tril(expA(:,1:end-1),-1));
% C2 = B2;
% 
% B2 = repmat(B2,N+1,1);
% C2 = reshape(C2,(N+1)*N,1);
% C2 = repmat(C2,1,N+1);
% C2 = reshape(C2',(N+1)^2,N);
% 
% X = (B1.*C1)-(B2.*C2);
% 
% Gm = sparse(G(1:end-1,3));
% % % % % C_0 = X*(G(1:end-1,3).^2);
% C_0 = X*(Gm.^2);
% C_0 = (c_alpha^2)*reshape(C_0,N+1,N+1);

%% the following part is for variable sigma computed with DE meal function

% D = tril(exp(-c_gamma*(repmat(G(:,1),1,N+1)-repmat(G(:,1)',N+1,1))));
% 
% int = zeros(N,1);
% C_0 = zeros(N+1,N+1);
% 
% for i = 1:N
% %     tt = (G(i,1):0.001:G(i+1,1));
% %     values = eval_integrand_for_cov_DE_meal(nutrition,G(1,1),G(i+1,1),tt,Gb,c_gamma,a,b);
% %     int1(i) = trapz(tt,values);
%     t_end = G(i+1,1);
%     g = @(t)eval_integrand_for_cov_DE_meal(nutrition,G(1,1),t_end,t,Gb,c_gamma,a,b);
%     int(i) = integral(g,G(i,1),G(i+1,1),'ArrayValued',true);
%     C_0(:,i+1) = exp(-c_gamma*(G(i+1,1)-G(i,1)))*C_0(:,i) + D(:,i+1)*int(i);
% end
% 
% C_0 = 2 * c_gamma * c_alpha^2 * C_0;
% 
% % if norm(C_0-C_0') > 1e-4 || isnan(norm(C_0-C_0'))
% %     warning('C_0 is not a symmetric matrix or it has nan entries')
% % end
% 
% C_0 = tril(C_0)+tril(C_0,-1)';
% 
% if norm(C_0-C_0') > 1e-4 || isnan(norm(C_0-C_0'))
%     warning('C_0 is not a symmetric matrix or it has nan entries')
% end

%% the following part is for sigma estimated directly
E1 = exp(-c_gamma*abs(repmat(G(:,1),1,N+1)-repmat(G(:,1)',N+1,1)));
E2 = exp(-c_gamma*(repmat(G(:,1),1,N+1)+repmat(G(:,1)',N+1,1)-2*G(1,1)));

% C_0 = c_sigma^2 * (E1-E2);
C_0 = (c_sigma^2 * (E1-E2)) + (E2 * g0_var);

%%
idx = find(G(:,6));
if idx(1) == 1
    idx(1) = [];
end
m0 = m_0(idx);
% C_0(1,1) = (0.1*g0)^2;
C0 = C_0(idx,idx);
%delta = 0.1*m0;z
S = C0;%+diag(delta.^2);
% S = C0+delta^2*eye(size(C0,1));
data = big_table.bg(idx);
bb = linsolve(S,data-m0);
log_det_S = sum(log(eig(S)));

phi_prop = ((data-m0).')*bb+log_det_S;

for dd=1:dim_meas
    outMSE(dd) = norm(data-m0);
    %outMSE(dd) = nan;
    outRcorr(dd) = nan;
end

chainTime = toc(chainTimer);
outerTime = toc(outer_timer);

end