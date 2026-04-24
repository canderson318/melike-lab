
%this scriot is to estimate parameters over a moving time interval of
%length (ideally) 1-2 days and then predict bg level until the next event
%time. Uses one beta estimate not two day/night ones.

%% Script Setup //////////////////////////////////////////////////////////

clear; clc;

% patient to model 
settings = jsondecode(fileread("settings.json"));
pat = settings.pat;

code_dir = '/Users/canderson/Documents/school/local-melike-lab/melike-lab/data_assimilation_main';

main_dir = '/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports';
master_data_dir = fullfile(main_dir,'data');

% Mount Onedrive
if ~exist("~/odrive/home")
    % For read/write access where MATLAB can create directories
    [status, output] = system('rclone mount odrive: ~/odrive --vfs-cache-mode full --daemon');
end

addpath(genpath(code_dir))
rmpath(genpath(fullfile(code_dir,'local')))


% load data
all_data = readtable(fullfile(main_dir,'Tandem_Tidepool_Deidentified.csv'),"NumHeaderLines" ,1);

all_data = all_data(:, 1:10);
names = ["patient_id","upload_date","sex","age","diagnosis_date","race","ethnicity","wt_kg","ht_cm","start_time",];
all_data.Properties.VariableNames = names;

%% Modelling Preferences //////////////////////////////////////////////////////////
setup_style = struct();
setup_style.data = 't1dm'; % used in createModelSettings
setup_style.model = 'OU_v1_T1D'; % used in createModelSettings 
setup_style.filter = 'EnKF_joint_v1';
setup_style.simulator = 'OU_v1';
setup_style.smoother = 'fmincon';


%% Model Settings  //////////////////////////////////////////////////////////

%  whers ther data
pat_data_dir = fullfile(master_data_dir,pat);

% filter and simulator settings
my_settings = struct();
my_settings.model = createModelSettings(setup_style,code_dir);
my_settings.model.f_nonlincon = [];
[DAstate,my_settings.filter] = createFilterSettings(setup_style);
my_settings.simulator = createSimulatorSettings(setup_style,DAstate,my_settings.model);

% smoother settings
t_prior = 0;
DAtype = 'smooth'; % run smoothing
my_settings.smoother = struct();
my_settings.smoother.patient = pat;
my_settings.smoother.theta_est_names = {'Gb','gamma','sigma','a','b','beta'};
if isfield(settings.smoother,"bounds")
    my_settings.smoother.lower_bounds = [settings.smoother.bounds.Gb(1) settings.smoother.bounds.gamma(1) settings.smoother.bounds.sigma(1) settings.smoother.bounds.a(1) settings.smoother.bounds.b(1) settings.smoother.bounds.beta(1)]';
    my_settings.smoother.upper_bounds = [settings.smoother.bounds.Gb(2) settings.smoother.bounds.gamma(2) settings.smoother.bounds.sigma(2) settings.smoother.bounds.a(2) settings.smoother.bounds.b(2) settings.smoother.bounds.beta(2)]';
else 
    % defaults
    my_settings.smoother.lower_bounds = [0 0.0001 0 0.01 0.01 10]';
    my_settings.smoother.upper_bounds = [1000 0.5 100 0.1 0.1 120]';
end

my_settings.smoother = createSmootherSettings(setup_style,DAstate,my_settings.model,my_settings.smoother);
my_settings.smoother.method = @Multi_start;
my_settings.smoother.Aineq = [0 0 0 1 -1 0];
my_settings.smoother.bineq = 0;
my_settings.smoother.is_use_parallel = 1;
my_settings.smoother.num_runs = 20;

 which_ind = find(strcmp(all_data.patient_id, pat));                                                      
if isempty(which_ind)                                                                                    
    error("error: Patient not found.")                                                                                               
end
age = all_data.age(which_ind);
sex = all_data.sex(which_ind);
wt_kg = all_data.wt_kg(which_ind);
ht_cm = all_data.ht_cm(which_ind);
my_settings.model.model_parameters.vg = EstimateBloodVolume_dL(age, sex, wt_kg, ht_cm);
my_settings.model.start_time = all_data.start_time(which_ind);
clear which ind age sex wt_kg ht_cm


% load patient data with setup preferences
orig_data = readData(pat_data_dir,setup_style);

% remove obs where no carbs
orig_data.raw_drivers.nutrition_oral(orig_data.raw_drivers.nutrition_oral.carbs == 0,:) = [];

%% Loop over diff window sizes ------
%+++ for different window params, loop over each of those windows and estimate smooth BG
%+++ save ewach outer loop res to estimation table

interval_length = settings.smoother.window.size;
interval_stride = settings.smoother.window.stride;

fprintf("\n****** Number of windows summary ******\n")
cumm_time = 0;
% how many windows? 
max_t = max(orig_data.raw_measurements.time);
% numbner of iterations 
num_iter = @(tmax,interval,stride) 1  + (( tmax - interval ) / stride);
num = floor(num_iter(max_t, interval_length, interval_stride));
fprintf("\t%.0f windows for %.0fmin window size at a %.0fmin stride\n",num,interval_length, interval_stride)
cumm_time = cumm_time + (interval_length*num); 
fprintf("\t>>> %.0fmin total patient time to model\n", cumm_time)


fprintf('\n**********************\nInterval Length: %.1fhr (%dmin)\tInterval Stride: %.1fhr (%dmin)\n**********************\n', interval_length/60, interval_length, interval_stride/60, interval_stride)

% my_settings.model.which_time = which_time;

% SET output directory;; IMPORTANT: do not let this be inside code repository.
if isfield(settings, 'output_dir')
    parent_output_dir = settings.output_dir    
else
    parent_output_dir = main_dir;
end

disp("Output Directory:")
output_dir = fullfile( parent_output_dir,strcat('moving_window_of_',num2str(interval_length),'mins_by_',num2str(interval_stride),'mins '))


%% Params for moving window --------------------

% patient window initial conditions
%starting point of estimation interval
t_0 = orig_data.raw_measurements.time(1); 

% end of interval
tt = t_0+interval_length;
% where observed time closest to end of interval/day
[~,closestIndex_1] = min(abs(tt-orig_data.raw_measurements.time));
% observeable end of interval
t_1 = orig_data.raw_measurements.time(closestIndex_1);

%% Loop over windows

% when to stop looping
% time_end = orig_data.raw_measurements.time(1) + interval_length + 60*12; % run until hit this time: first time + windowsize + 6hours -> if stride is 1hr this will run on 6 windows
num_windows = settings.smoother.window.num_windows
if ~strcmp(num_windows, "all") &&  ~isnumeric( num_windows )
    error(sprintf('%d is an invalid argument.\nUse "all" for all of the data or an integer for the number of windows to run.', num_windows))
end
% time_end = 12*60 + interval_length % run for 12 hours

time_end = orig_data.raw_measurements.time(end); % run for all data
counter = 1;
data = orig_data;

while t_1 <= time_end && (strcmp(num_windows, "all") | counter <= num_windows)
    
    fprintf('\n\t**********************\n\t\t(%.0f) Window: %.1fhr (%dmin) - %.1fhr (%dmin)\n\t**********************\n',i,t_0/60, t_0, t_1/60, t_1)

    
    run_output_dir = fullfile(output_dir,pat,sprintf('smooth_%s_%s',int2str(t_0),int2str(t_1)));

    %find the time and value of the latest measurement before t_0
    indd = find(data.raw_measurements.time <= t_0);
    indx = indd(end);
    bg_val = data.raw_measurements.bg(indx);

    %use this measurement value to pick initial point g_0
    my_settings.model.model_parameters.g0 = bg_val;

    % subset data for time window 
    data.raw_drivers.nutrition_oral = data.raw_drivers.nutrition_oral((data.raw_drivers.nutrition_oral.time <= t_1),:);
    data.raw_drivers.insulin_basal = data.raw_drivers.insulin_basal((data.raw_drivers.insulin_basal.time <= t_1),:);
    data.raw_drivers.insulin_bolus = data.raw_drivers.insulin_bolus((data.raw_drivers.insulin_bolus.time <= t_1),:);
    data.raw_measurements = data.raw_measurements(((data.raw_measurements.time >= t_0) & (data.raw_measurements.time <= t_1)),:);

    % update settings for this window
    my_settings.smoother.t_opt_start = t_0;
    my_settings.smoother.t_opt_end = t_1;

    % Run smoother 
    solver_dir = prep_main(data,my_settings,t_prior,DAtype,DAstate,run_output_dir,code_dir);
    [output,~] = main(data,my_settings,t_prior,DAtype,DAstate,run_output_dir,code_dir,solver_dir);

    % reset data 
    data = orig_data;

    % set next window_start as current_start + stride
    t_00 = t_0+interval_stride;
    [~,closestIndex_2] = min(abs(t_00-data.raw_measurements.time)); % patient obs near time
    t_0_new = data.raw_measurements.time(closestIndex_2); % new window start
    
    if t_0_new == t_0 %there is no measurement in the interval t_0 to t_0+interval_stride
        t_0 = data.raw_measurements.time(closestIndex_2+1);
    else
        t_0 = t_0_new;
    end
    clear t_0_new

    % set end of window
    tt = t_0+interval_length;
    [~,closestIndex_1] = min(abs(tt-data.raw_measurements.time));
    t_1 = data.raw_measurements.time(closestIndex_1);
    if t_1-t_0 < 110 % window getting small
        if closestIndex_1+1 > length(data.raw_measurements.time)
            warning("Warning: index exceeds array length.")
            % break
            t_1 = inf;
        else 
            t_0 = data.raw_measurements.time(closestIndex_1+1);
            clear closestIndex_1
            tt = t_0+interval_length;
            [~,closestIndex_1] = min(abs(tt-data.raw_measurements.time));
            t_1 = data.raw_measurements.time(closestIndex_1);
        end
    end
    counter = counter+1;

end

