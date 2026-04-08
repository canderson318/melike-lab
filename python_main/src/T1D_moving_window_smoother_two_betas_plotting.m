% open repl
% matlab -nodesktop

%% setup
clear;clc;
set(groot, 'defaultFigureColor', 'white')
set(groot, 'defaultAxesColor', 'white')
set(groot, 'defaultTextColor', 'black')
set(groot, 'defaultAxesXColor', 'black')
set(groot, 'defaultAxesYColor', 'black')
set(groot, 'defaultColorbarColor', 'black') 

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
main_dir = '/Users/canderson/odrive/home/melike-rotation/project001/Tidepool_Exports';

pat = "SM002";

dirs = string({dir(main_dir).name});
window_sets = dirs(contains(dirs, "moving_window"));

% which file to read from each dir
f = "optimization_summary.mat";

summaries =  cell(length(window_sets), 1);

% w = 1;
for w=1:length(window_sets)
    window_set = window_sets(w);
    fprintf("%s\n", window_set)
    parent_data_dir = fullfile(main_dir,fullfile(window_set,pat));


    %% Load Data 

    % each windows result dir
    data_listing = dir(parent_data_dir);
    data_dirs = data_listing(~startsWith({data_listing.name}, '.'));
    data_names = {data_dirs.name};  % Just the names as cell array



    % what fields to grab
    fields = ["interval_start", "interval_end","Gb", "gamma", "sigma", "a", "b", "beta_day", "beta_night", "f_val","rmse"];
    types = repelem("double",length(fields));
    summ = table('Size', [length(data_names), length(fields)], 'VariableTypes',types, 'VariableNames', fields);

    % i = 1;
    for i=1:length(data_names)

        a_data_dir = char(data_names(i));

        splote = strsplit(a_data_dir, "_");
        [interval_start interval_end ] = deal(str2num(splote{2}),str2num(splote{3}));
        
        optimization_summary = load(fullfile(parent_data_dir, a_data_dir, f));

        summ(i,:) = {interval_start...
                    interval_end...
                    optimization_summary.output.best_params.Optimal(1)...%Gb
                    optimization_summary.output.best_params.Optimal(2)... %gamma
                    optimization_summary.output.best_params.Optimal(3)... %sigma
                    optimization_summary.output.best_params.Optimal(8)... %a
                    optimization_summary.output.best_params.Optimal(9)... %b
                    optimization_summary.output.best_params.Optimal(10)... %beta_day
                    optimization_summary.output.best_params.Optimal(11)... %beta_night
                    optimization_summary.output.fval_per_meas... 
                    sqrt(optimization_summary.output.mse)};
        
        fprintf("%s%%\n", num2str(round((i/length(data_names) )*100)) )
    end

    [~,order] = sort(summ.interval_start);
    summ_ordered = summ(order,:);
    summaries{w} = summ_ordered;

end % end outer loop


%% Plot for each window set
for i=1:length(window_sets)
    a_summ = summaries{i};
    window = window_sets(i);

    if ~ exist(fullfile(main_dir, "plots", window))
        mkdir(fullfile(main_dir, "plots", window))
    end

    figure; set(gcf, 'Color', 'white'); set(gca, 'Color', 'white'); 
    plot(a_summ.interval_start, a_summ.rmse)
    xlabel("Interval start")
    title(window, 'Interpreter', 'none')
    ylabel("RMSE")
    xt = get(gca, 'XTick');
    xticklabels(datestr(xt/60/24, 'HH:MM'))
    saveas(gcf, fullfile(main_dir, "plots", window, 'rmse_against_interval_start.png'))


    figure; set(gcf, 'Color', 'white'); set(gca, 'Color', 'white'); 
    plot(a_summ.interval_start, a_summ.gamma); hold on
    scatter(a_summ.interval_start, a_summ.gamma, 100, a_summ.rmse, 'filled')
    xt = get(gca, 'XTick');
    xticklabels(datestr(xt/60/24, 'HH:MM'))
    xlabel("Interval Start")
    ylabel("Gamma")
    title(window, 'Interpreter', 'none')
    colorbar
    hold off 
    saveas(gcf, fullfile(main_dir, "plots", window, 'gamma_against_interval_start_rmse_filled.png'))

end


%% PLot gammas at each time of day
cols = ["interval_start", "interval_end", "gamma"];
tmp = cellfun(@(s) s(:, cols), summaries, 'UniformOutput', false);
gammas = vertcat(tmp{:});
gammas.window = gammas.interval_end - gammas.interval_start;

x = [gammas.interval_start, gammas.interval_end, nan(height(gammas),1)]';
y = [gammas.gamma, gammas.gamma, nan(height(gammas),1)]';

figure; set(gcf, 'Color', 'white'); set(gca, 'Color', 'white'); 
plot(x(:),y(:)); hold on
scatter(x(1,:),y(1,:), 100, gammas.window, "filled"); hold on
scatter(x(2,:),y(2,:), 100, gammas.window); hold on
cb = colorbar;
cb.Label.String = "Interval Size";
ct = cb.Ticks;                                                                   
cb.TickLabels = datestr(ct/60/24, 'HH:MM');
hold off
xlabel("Time of Day")
ylabel("Gamma")
xt = get(gca, 'XTick');
xticklabels(datestr(xt/60/24, 'HH:MM'))
title("Gamma Estimate For Different Intervals")
saveas(gcf, fullfile(main_dir, "plots", 'gammas_at_diff_start_times.png'))


