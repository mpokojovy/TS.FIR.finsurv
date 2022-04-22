%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A Fast Initial Response Approach to Real-Time Financial Surveillance  %
%            (C) Michael Pokojovy and Andrews T. Anum (2022)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings
est_type = 'usual'; % Choose from 'mcd', 'mdpd' or 'usual'
alpha = 0.6; % Select TS parameter alpha from 3/4, 2/3, 5/8, 3/5, 1
all_ok_ARL = 370; % Select all-ok ARL from 100, 370 or 750
CUSUM_k = 0.5; % Select k from [0.1000  0.2000  0.2500  0.3000  0.4000  0.5000  0.6000  0.7000  0.7500  0.8000  0.9000  1.0000  1.1000  1.2000  1.2500  1.3000  1.4000  1.5000]
    
EDA_plot_flag = true;

%% Load training data
data2019 = readtable("data\FB_JAN1_2019_FEB7_2020.csv");

I = 1:size(data2019, 1);
prices2019 = data2019{flip(I), 5};
dates2019 = data2019{flip(I), 1};

T_obs = days(dates2019 - min(dates2019)) + 1;
T = (1:(days(max(dates2019) - min(dates2019)) + 1))';

dates2019  = (min(dates2019):days(1):max(dates2019))';
prices2019 = interp1(T_obs, prices2019, T);

%% Load monitoring data
data2020 = readtable("data\FB_FEB10_2020_MAR31_2020.csv");

I = 1:size(data2020, 1);
prices2020 = data2020{flip(I), 5};
dates2020 = data2020{flip(I), 1};

T_obs = days(dates2020 - min(dates2020)) + 1;
T = (1:(days(max(dates2020) - min(dates2020)) + 1))';

dates2020  = (min(dates2020):days(1):max(dates2020))';
prices2020 = interp1(T_obs, prices2020, T);

%% Plot raw prices
if (EDA_plot_flag)
    plot_raw_prices([dates2019; dates2020], [prices2019; prices2020], dates2020(1), 'Meta Platforms, Inc.~(FB)');
end
%% Plot log-differences
if (EDA_plot_flag)
    plot_raw_prices([dates2019; dates2020], [0; diff(log(prices2019)); diff(log([prices2019(end); prices2020]))], ...
                    dates2020(1), 'Meta Platforms, Inc.~(FB)');
end                

%% KS test to estimate the degrees of freedom in Student's t-distribution
x_train_raw = diff(log(prices2019));     

dfs = [1, linspace(2, 3, 11), 4 5];
[pvals, df_best] = student_KS_test(x_train_raw, dfs);

display(['df''s: ', num2str(dfs)]);
display(['pval''s: ', num2str(pvals)]);
display(['df best = ', num2str(df_best)]);

if (EDA_plot_flag)
    dfs = linspace(0.05, 5, 100); 
    [pvals, ~] = student_KS_test(x_train_raw, dfs);
    
    plot_pvals(dfs, pvals);
end
            
%% Student's t Q-Q-plot for log-differences
dfs = [1 2 2.7 3 4 5];

if (EDA_plot_flag)
    plot_student_QQ(x_train_raw, dfs, 2, 3);
end    

%% Perform Gaussianization of training data
loc_train = median(x_train_raw);
c_train = (tinv(0.75, df_best) - tinv(0.25, df_best))/(quantile(x_train_raw, 0.75) - quantile(x_train_raw, 0.25));
x_train = norminv(tcdf((x_train_raw - loc_train)*c_train, df_best));

%% Plot Gaussianized log-differences
log_diffs = [0; diff(log(prices2019)); diff(log([prices2019(end); prices2020]))];
log_diffs_gaussianized = norminv(tcdf((log_diffs - loc_train)*c_train, df_best));
if (EDA_plot_flag)
    plot_gaussianized_log_differences([dates2019; dates2020], log_diffs_gaussianized, ...
                                      dates2020(1), 'Meta Platforms, Inc.~(FB)');
end

%% Autocorrelation plot
if (EDA_plot_flag)
    plot_ACF(log_diffs_gaussianized);
end    

%% Estimate IC parameters
est_type = lower(est_type);

if (strcmp(est_type, 'mdpd'))
    par = 0.1182555;
elseif (strcmp(est_type, 'mcd'))
    par = 0.1;
else
    par = [];
end

[mu0, sigma0] = estimate_IC_standards(x_train, est_type, par);

display(['IC standards: ', num2str(mu0), ' and ', num2str(sigma0)]);

%% Plot standardized Gaussianized log-differences
log_diffs_gaussianized_standardized = (log_diffs_gaussianized - mu0)/sigma0;
if (EDA_plot_flag)
    plot_standardized_gaussianized_log_differences([dates2019; dates2020], log_diffs_gaussianized_standardized, ...
                                                   dates2020(1), 'Meta Platforms, Inc.~(FB)');  
end                                               
                                           
%% Prepare Process Monitoring
x_test = diff(log([prices2019(end); prices2020]));
x_test = norminv(tcdf((x_test - loc_train)*c_train, df_best));
x_test = (x_test - mu0)/sigma0;

run_control_charts(dates2020, x_test, CUSUM_k, alpha, all_ok_ARL);