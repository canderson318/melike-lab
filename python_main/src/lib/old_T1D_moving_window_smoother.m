%this scriot is to estimate parameters over a moving time interval of
%length (ideally) 1-2 days and then predict bg level until the next event
%time
clear
%% initialize the settings

foo = pwd;
token = '/data_assimilation_main';
code_dir = foo(1:strfind(pwd,token)+length(token)-1);

if ~exist(code_dir,'dir')
    error('%s could not be found. Make sure you run this script downstream of it.\nYou are currently running from: %s',token,foo);
end

main_dir = '~/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Tidepool_Exports';
master_data_dir = fullfile(main_dir,'data');
master_output_dir = fullfile(main_dir,'moving_window_of4hours_by30mins');

addpath(genpath(code_dir))
rmpath(genpath(fullfile(code_dir,'local')))

all_data = readtable(fullfile(main_dir,'Tandem_Tidepool_Deidentified.xlsx'));

setup_style = struct();
setup_style.data = 't1dm'; % used in createModelSettings
setup_style.model = 'OU_v1_T1D'; % used in createModelSettings 
setup_style.filter = 'EnKF_joint_v1';
setup_style.simulator = 'OU_v1';
setup_style.smoother = 'fmincon';

pat = char('SM001');
data_dir = fullfile(master_data_dir,pat);

% SET output directory.
% IMPORTANT: do not let this be inside code repository.
% master_output_dir = fullfile('~/OU_ICU_prediction');

my_settings = struct();

my_settings.model = createModelSettings(setup_style,code_dir);
my_settings.model.f_nonlincon = [];

[DAstate,my_settings.filter] = createFilterSettings(setup_style);
my_settings.simulator = createSimulatorSettings(setup_style,DAstate,my_settings.model);

my_settings.smoother = struct();
my_settings.smoother.patient = pat;

data = readData(data_dir,setup_style);

times = [data.raw_measurements.time;data.raw_drivers.nutrition_oral.time;...
    data.raw_drivers.insulin_basal.time;data.raw_drivers.insulin_bolus.time];
times = unique(sort(times));
count_trend = 0;

DAtype = 'smooth';

data = readData(data_dir,setup_style);
data.raw_drivers.nutrition_oral(data.raw_drivers.nutrition_oral.carbs == 0,:) = [];

length_of_interval = 1440/6;
length_of_moving = 30;

estimation_table = zeros(1,10);

t_0 = data.raw_measurements.time(1); %starting point of estimation interval
tt = t_0+length_of_interval;
[~,closestIndex_1] = min(abs(tt-data.raw_measurements.time));
t_1 = data.raw_measurements.time(closestIndex_1);

counter = 1;

while t_1 <= data.raw_measurements.time(1) + (1440*8)
    run_output_dir = fullfile(master_output_dir,pat,sprintf('smooth_%s_%s',int2str(t_0),int2str(t_1)));

    %find the time and value of the latest measurement before t_0
    indd = find(data.raw_measurements.time <= t_0);
    indx = indd(end);
    bg_val = data.raw_measurements.bg(indx);

    %use this measurement value to pick initial point g_0
    my_settings.model.model_parameters.g0 = bg_val;

    data.raw_drivers.nutrition_oral = data.raw_drivers.nutrition_oral((data.raw_drivers.nutrition_oral.time <= t_1),:);
    data.raw_drivers.insulin_basal = data.raw_drivers.insulin_basal((data.raw_drivers.insulin_basal.time <= t_1),:);
    data.raw_drivers.insulin_bolus = data.raw_drivers.insulin_bolus((data.raw_drivers.insulin_bolus.time <= t_1),:);
    data.raw_measurements = data.raw_measurements(((data.raw_measurements.time >= t_0) & (data.raw_measurements.time <= t_1)),:);

    % % my_settings.simulator.t_sim_interval = 2;

    my_settings.smoother.t_opt_start = t_0;
    my_settings.smoother.t_opt_end = t_1;
    my_settings.smoother.theta_est_names = {'Gb','gamma','sigma','a','b','beta'};
    my_settings.smoother.lower_bounds = [0 0.0001 0 0.01 0.01 10]';
    my_settings.smoother.upper_bounds = [1000 0.5 100 0.1 0.1 120]';
    my_settings.smoother = createSmootherSettings(setup_style,DAstate,my_settings.model,my_settings.smoother);
    % my_settings.model.model_parameters.beta = 0;
    % my_settings.model.model_parameters.a = 0.04;
    % my_settings.model.model_parameters.b = 0.06;
    my_settings.smoother.method = @Multi_start;
    my_settings.smoother.Aineq = [0 0 0 1 -1 0];
    my_settings.smoother.bineq = 0;
    % my_settings.smoother.Aineq = [0 0 0 1 -1 0; 1 0 0 0 0 -3.5; 0 -1110 0 0 0 1];
    % my_settings.smoother.bineq = [0; 115; 10];

    my_settings.model.f_nonlincon = []; %@nonlin_const_DE_meal_params;
    which_ind = find(strcmp(all_data.patient_id,pat));
    age = all_data.age(which_ind);
    sex = all_data.sex(which_ind);
    weight_kg = all_data.weight_kg(which_ind);
    height = all_data.height_cm(which_ind);
    my_settings.model.model_parameters.vg = EstimateBloodVolume_dL(age, sex, weight_kg, height);
    clear which ind age sex weight_kg height
    my_settings.smoother.is_use_parallel = 1;

    my_settings.smoother.num_runs = 20;

    t_prior = 0;
    solver_dir = prep_main(data,my_settings,t_prior,DAtype,DAstate,run_output_dir,code_dir);
    [output,~] = main(data,my_settings,t_prior,DAtype,DAstate,run_output_dir,code_dir,solver_dir);
    
    estimation_table(counter,1) = output.best_params.Optimal(1);%Gb
    estimation_table(counter,2) = output.best_params.Optimal(2);%gamma
    estimation_table(counter,3) = output.best_params.Optimal(3);%sigma
    estimation_table(counter,4) = output.best_params.Optimal(8);%a
    estimation_table(counter,5) = output.best_params.Optimal(9);%b
    estimation_table(counter,6) = output.best_params.Optimal(10);%beta
    estimation_table(counter,7) = output.fval_per_meas;
    estimation_table(counter,8) = sqrt(output.mse);
    estimation_table(counter,9) = t_0;
    estimation_table(counter,10) = t_1;
    estimation_table(counter,11) = t_1-t_0;

    data = readData(data_dir,setup_style);
    data.raw_drivers.nutrition_oral(data.raw_drivers.nutrition_oral.carbs == 0,:) = [];

    t_00 = t_0+length_of_moving;
    [~,closestIndex_2] = min(abs(t_00-data.raw_measurements.time));
    t_0_new = data.raw_measurements.time(closestIndex_2);
    if t_0_new == t_0 %there is no measurement in the interval t_0 to t_0+length_of_moving
        t_0 = data.raw_measurements.time(closestIndex_2+1);
    else
        t_0 = t_0_new;
    end
    clear t_0_new
    tt = t_0+length_of_interval;
    [~,closestIndex_1] = min(abs(tt-data.raw_measurements.time));
    t_1 = data.raw_measurements.time(closestIndex_1);
    if t_1-t_0 < 110 %there is a measurement gap after t_1 and estimation windows getting smaller
        t_0 = data.raw_measurements.time(closestIndex_1+1);
        clear closestIndex_1
        tt = t_0+length_of_interval;
        [~,closestIndex_1] = min(abs(tt-data.raw_measurements.time));
        t_1 = data.raw_measurements.time(closestIndex_1);
    end
    counter = counter+1;
end

estimation_table = array2table(estimation_table,'VariableNames',...
    {'Gb','gamma','sigma','a','b','beta','fval_per_meas','rmse','t0','t1','length'});

save(fullfile(master_output_dir,pat,'estimation_table.mat'),'estimation_table')
writetable(estimation_table,fullfile(master_output_dir,pat,'estimation_table.csv'))

figure('Renderer','painters', 'Position',[10 10 1000 1000]);
subplot(8,1,1)
plot(estimation_table.t0,estimation_table.Gb,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('Gb')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,2)
plot(estimation_table.t0,estimation_table.gamma,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('gamma')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,3)
plot(estimation_table.t0,estimation_table.sigma,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('sigma')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,4)      
plot(estimation_table.t0,estimation_table.a,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('a')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,5)
plot(estimation_table.t0,estimation_table.b,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('b')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,6)
plot(estimation_table.t0,estimation_table.beta,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('beta')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,7)
plot(estimation_table.t0,estimation_table.fval_per_meas,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('fval per measurement')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)

subplot(8,1,8)
plot(estimation_table.t0,estimation_table.rmse,'-*','LineWidth',1.5)
xlim([estimation_table.t0(1),estimation_table.t0(end)])
xticks(estimation_table.t0(1):180:estimation_table.t0(end))
ylabel('rmse')
xt = get(gca, 'XTick');
labels = datestr(xt/60/24, 'HH:MM');
xticklabels(labels)
xlabel('Time (hh:mm)')

%% To understand and temporal patterns in beta estimation

% % % % [ind,meanVals] = ischange(estimation_table.beta,'mean');
% % % % plot(estimation_table.beta);
% % % % hold on;
% % % % plot(ind,estimation_table.beta(ind),'ro')
% % % % 
% % % % time_hr = estimation_table.t0/60;
% % % % time_of_day = mod(time_hr,24);
% % % % 
% % % % figure;
% % % % scatter(time_of_day, estimation_table.beta, 'filled')
% % % % xlabel('Hour of day')
% % % % ylabel('Value')
% % % % title('Time-of-day pattern over 2 days')
% % % % 
% % % % edges = 0:1:24;  % 1-hour bins
% % % % [N, ~, bin] = histcounts(time_of_day, edges);  % bin indices
% % % % % Remove values that fall outside (bin = 0)
% % % % valid = bin > 0;
% % % % bin = bin(valid);
% % % % data_valid = estimation_table.beta(valid);
% % % % 
% % % % % Compute mean per hour
% % % % mean_vals = accumarray(bin, data_valid, [length(edges)-1, 1], @mean);
% % % % 
% % % % % Plot
% % % % figure;
% % % % bar(edges(1:end-1)+0.5, mean_vals)
% % % % xlabel('Hour of day')
% % % % ylabel('Mean value')
% % % % title('Average value per hour of day')tiew
% % % % 
% % % % smoothed = smooth(estimation_table.beta, 180); % e.g., 1-hour window
% % % % plot(time_of_day, smoothed, 'LineWidth', 1.5)
% % % % xlabel('Hour of day')
% % % % ylabel('Smoothed value')