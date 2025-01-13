clear; clc; close all

%% Define variables

data_classes = ["Linguistic", "SRT_BinShan", "SRT_Joanne"];
data_class = data_classes(1);

do_norm = 0; 
norm_range = [-1, 1];

if do_norm == 1, out_tags = " (normed)"; else, out_tags = ""; end

%% Setup file paths

if data_class == "Linguistic" % -------------------------------------------

    authors = ["Chang_et_al", "Tse_et_al"]; 
    author = authors(1); 
    
    tasks = ["Naming", "LD"]; 
    task = tasks(1); 

    data_folder = fullfile("Data_Linguistic", author, task); 

    data_vers = ["raw", "zscored"];
    data_ver = data_vers(2);

    if data_ver == "zscored"
        fn_regex = "zscored_sub_*.xlsx";
        sid_regex = "zscored_sub_%d.xlsx"; 
        out_tags = " (z-scored)" + out_tags; 
    else
        fn_regex = "sub_*.xlsx"; 
        sid_regex = "sub_%d.xlsx"; 
    end

    out_folder = fullfile("Output_Linguistic", author, task);    
    
elseif data_class == "SRT_BinShan" % --------------------------------------

    exp_names = [
        "Exp1_Spatial", ...
        "Exp2_Temporal", ...
        "Exp3_Intertwine(tDCS)", ...
        "Exp4_Intertwine"]; 
    exp_name = exp_names(1);

    data_vers = [
        ".", ... 
        "Basic_Shannon (2-back)", ...
        "Basic_Shannon (3-back)"];
    data_ver = data_vers(3);

    data_folder = fullfile("Data_SRT_BinShan", ...
        exp_name, "PreprocessedData", data_ver); 

    if data_ver == "."
        fn_regex = '*_Spatial_HL.csv';
        sid_regex = '%d_Spatial_HL.csv';

    else    
        fn_regex = 'sub_*.xlsx';
        sid_regex = 'sub_%d.xlsx';
    end

    out_folder = fullfile("Output_SRT_BinShan", ...
        exp_name, data_ver); 

elseif data_class == "SRT_Joanne" % ---------------------------------------
    
    inst_types = ["explicit", "implicit"]; 
    inst_type = inst_types(1);

    subj_groups = ["learned", "not_learned"];
    subj_group = subj_groups(1);

    data_vers = [
        "Basic_Shannon (2-back)", ...
        "Basic_Shannon (3-back)"];
    data_ver = data_vers(2);

    data_folder = fullfile("Data_SRT_Joanne", ...
        "PreprocessedData", data_ver, inst_type, subj_group); 

    fn_regex = 'sub_*.xlsx';
    sid_regex = 'sub_%d.xlsx';

    out_folder = fullfile("Output_SRT_Joanne", ...
        data_ver, inst_type, subj_group); 

end

% -------------------------------------------------------------------------

if ~ exist(out_folder, 'dir'), mkdir(out_folder), end

out_file_1 = join(['GLM_regression_results', out_tags, ".csv"], ''); 
out_file_2 = join(['GLM_residuals', out_tags, ".csv"], ''); 

%% Load data and preallocating memory

cd 'C:\Users\PinWei\my_CSR'; 

data_paths = dir(fullfile(data_folder, fn_regex)); 

% Find the maximum number of trials:
nT_max = 0; 
for i = 1:length(data_paths)
    data = readtable(fullfile(data_folder, data_paths(i).name));
    nT_max = max(nT_max, size(data, 1));  
end

nS = length(data_paths);    % number of subjects
nF = size(data, 2) - 1;     % number of factors
nP = nF + 1;                % number of parameters

model_measures = [
    "R_squared", "Adjusted_R2", "LogLik", "AIC", "AICc", "BIC", "NRMSE"];
nM = length(model_measures);

% Preallocating memory:
sid_list = zeros(nS, 1); 
results_list = zeros(nS, 1 + nP + nM);   
residuals_list = NaN(nT_max, nS);

%% Loop through each file

for i = 1:length(data_paths)

    fn = fullfile(data_folder, data_paths(i).name);
    data = readtable(fn); 

    sid = sscanf(data_paths(i).name, sid_regex); 
    sid_list(i) = sid;

    % The last column is the response variable Y:
    Y_name = data.Properties.VariableNames{end};    
    Y = data{:, end};

    % All other columns are factors:
    F_names = data.Properties.VariableNames(1:end-1); 
    F_vals = data{:, 1:end-1}; 

%     % Remove trials with extreme values (z-score > 3):
%     if exclude_outlier == 1
%         rows_to_remove = (F_vals(:, 1) > 3) | (F_vals(:, 2) > 3);
%         Y = Y(~ rows_to_remove, :);
%         F_vals = F_vals(~ rows_to_remove, :);
%     end

    % Normalization:
    if do_norm == 1, F_vals = normalize(F_vals, 'range', norm_range); end

    % Create the design matrix:
    nT = length(Y);       % number of trials
    nI = sum(nF-1:-1:1);  % number of interaction terms
    X = [  
        ones(nT, 1),  ... % constant term
        F_vals,       ... % linear terms
    ];  

    % Solve linear equation: 
    Coefs = pinv(X' * X) * X' * Y; 
    Y_pred = X * Coefs; 
    residuals = Y - Y_pred; 

    % Caculate R squared and Adjusted R squared:
    RSS = sum(residuals .^ 2);     % Residual Sum of Squares
    TSS = sum((Y - mean(Y)) .^ 2); % Total Sum of Squares
    Rsq = 1 - (RSS / TSS); 
    Rsq_adj = 1 - (RSS / (nT - nP - 1)) / (TSS / (nT - 1));  
        
    % Caculate various Information Criterion: 
    re_var = RSS / nT; % Residual variance
    log_lik = -0.5 * nT * (log(2 * pi * re_var) + 1); % Log Likelihood
    [~, ~, IC] = aicbic(log_lik, nP, nT); 
%     AIC = -2 * log_lik + 2 * nP;       % Akaike IC
%     AICc = AIC + (2 * nP ^ 2 + 2 * nP) / (nT - nP - 1); % Corrected AIC
%     BIC = -2 * log_lik + log(nT) * nP; % Bayesian IC
    
    % Caculate Normalized Root-Mean-Square Error:
    NRMSE = sqrt(mean((residuals ./ Y_pred) .^ 2)); 

    % Store the parameters and other results of each subject :
    param_info = [sid, Coefs', Rsq, Rsq_adj, ...
        log_lik, IC.aic, IC.aicc, IC.bic, NRMSE];
%     param_info = [sid, Coefs', Rsq, Rsq_adj, log_lik, AIC, AICc, BIC, NRMSE];
    results_list(i, :) = param_info;

    % Fill the residuals into column i (replace NaN)
    residuals_list(1:length(residuals), i) = residuals; 
end

fprintf('Average R-squared = %.3f', mean(results_list(:, end)));

%% Generate variable names and save to file

var_names = [
    'SID', 'X0', F_names, model_measures];
results_list_table = array2table( ...
    results_list, 'VariableNames', var_names);
writetable(results_list_table, fullfile(out_folder, out_file_1));

residuals_table = array2table( ...
    residuals_list, 'VariableNames', strcat("sub_", string(sid_list)));
writetable(residuals_table, fullfile(out_folder, out_file_2));
