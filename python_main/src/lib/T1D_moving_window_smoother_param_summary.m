% Create parameter summaries for each window 
% save to csv files
%% setup
clear;clc;

% patient 
pat = fileread('PAT.txt');
% pat = "SM001";
fprintf("%s\n",pat)

% add scripts to path
code_dir = '/Users/canderson/Documents/school/local-melike-lab/melike-lab/data_assimilation_main';
addpath(genpath(code_dir))
% addpath(fullfile(code_dir, "projects/T1DM"));  
rmpath(genpath(fullfile(code_dir,'local')))


% Mount Onedrive
if ~exist("~/odrive/home")
    % For read/write access where MATLAB can create directories
    [status, output] = system('rclone mount odrive: ~/odrive --vfs-cache-mode full --daemon');
end

% input dir
in_dir = '/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports';

% output dir
out_main_dir = "/Users/canderson/odrive/home/melike-rotation/project001/outputs/";
if ~exist(out_main_dir)
    mkdir(out_main_dir)
end

% dir for this script
out_dir = fullfile(out_main_dir, "param_summary");
if ~exist(out_dir)
    mkdir(out_dir)
end


dirs = string({dir(in_dir).name});
window_sets = dirs(contains(dirs, "moving_window"));

% which file to read from each dir
f = "optimization_summary.mat";

summaries =  cell(length(window_sets), 1);


% w = 3;
for w=1:length(window_sets)
    window_set = window_sets(w);
    parent_data_dir = fullfile(in_dir,fullfile(window_set,pat));
    
    if ~exist(parent_data_dir)
        continue
    end

    fprintf("\n%s\n", window_set)

    %% Load Data 
    
    % each windows result dir
    data_listing = dir(parent_data_dir);
    data_dirs = data_listing(~startsWith({data_listing.name}, '.'));
    data_names = {data_dirs.name};  % Just the names as cell array

    % what fields to grab
    
    % fields = ["interval_start", "interval_end","Gb", "gamma", "sigma", "a", "b", "beta_day", "beta_night", "beta","f_val","rmse"];
    % types = repelem("double",length(fields));
    % summ = table('Size', [length(data_names), length(fields)], 'VariableTypes',types, 'VariableNames', fields);
    
    x_summary_path = fullfile(parent_data_dir, char(data_names(1)), f);
    x_optimization_summary = load(x_summary_path);
    
    % fixed fields always present
    base_fields = ["interval_start", "interval_end", "f_val", "rmse"];
    
    % param fields — intersect with actual RowNames
    param_candidates = ["Gb","gamma","sigma","delta","g0","vg","T","a","b","beta","beta_d","beta_n","g0_var"];
    row_names = x_optimization_summary.output.best_params.Properties.RowNames;
    param_fields = intersect(param_candidates, string(row_names), 'stable');
    
    % combine and build table
    fields = [base_fields, param_fields.'];
    types = repelem("double", length(fields));
    summ = table('Size', [length(data_names), length(fields)], ...
                'VariableTypes', types, ...
                'VariableNames', fields);

    % i = 1;
    for i=1:length(data_names)

        a_data_dir = char(data_names(i));

        splote = strsplit(a_data_dir, "_");
        [interval_start interval_end ] = deal(str2num(splote{2}),str2num(splote{3}));
        
        summary_path = fullfile(parent_data_dir, a_data_dir, f);
        if exist(summary_path)
            optimization_summary = load(summary_path);
            bp = optimization_summary.output.best_params;
            row_vals = [{interval_start, interval_end, ...
                optimization_summary.output.fval_per_meas, ...
                sqrt(optimization_summary.output.mse)}, ...
                arrayfun(@(p) bp{char(p), 'Optimal'}, param_fields, 'UniformOutput', false).'];

            summ(i,:) = row_vals;
        else
            fprintf("\nOptimization summary not found at %s", summary_path)
        end
        fprintf("\r%s%%", num2str(round((i/length(data_names) )*100)))
    end

    [~,order] = sort(summ.interval_start);
    summ_ordered = summ(order,:);
    summ_ordered.window_name = repmat(window_set, height(summ_ordered), 1);
    summaries{w} = summ_ordered;

end % end outer loop

summaries_tab = vertcat(summaries{:});
summaries_tab.pat = repmat(pat, height(summaries_tab), 1);
summaries_tab.interval_size = summaries_tab.interval_end - summaries_tab.interval_start;


%% SAVE SUMMARY
out_path_dir = fullfile(out_dir,pat);
if ~exist(out_path_dir)
    mkdir(out_path_dir)
end
writetable(summaries_tab, fullfile(out_path_dir,"summary.csv")) % "/Users/canderson/odrive/home/melike-rotation/project001/outputs/param_summary/SM002"

