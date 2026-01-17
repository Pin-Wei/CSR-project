clear; clc; close all

cfg = define_config();

data_paths = dir(fullfile(cfg.data_folder, cfg.fn_regex)); 

if ~ exist(cfg.out_folder, 'dir'), mkdir(cfg.out_folder), end

model_types = ["CSR", "GLM"];
model_type = model_types(2);

regression_results_out = fullfile(cfg.out_folder, join([model_type, '_', 'regression_results', cfg.out_tags, ".csv"], '')); 
residuals_table_out = fullfile(cfg.out_folder, join([model_type, '_', 'residuals', cfg.out_tags, ".csv"], ''));

%% Load data and preallocating memory

% Find the maximum number of trials:
nT_max = 0; 
for i = 1:length(data_paths)
    fp = fullfile(cfg.data_folder, data_paths(i).name);
    data = readtable(fp);
    nT_max = max(nT_max, size(data, 1));  
end

nS = length(data_paths); % number of subjects
nF = size(data, 2) - 1; % number of factors

switch model_type
    case "CSR"
        nP = (nF*3 + nF^2 + 2) / 2; % number of parameters

    case "GLM"
        nP = nF + 1;
end

model_measures = ["R_squared", "Adjusted_R2", "LogLik", "AIC", "AICc", "BIC", "NRMSE"];
nM = length(model_measures);

% Preallocating memory:
sid_list = zeros(nS, 1); 
results_list = zeros(nS, 1 + nP + nM);  
residuals_list = NaN(nT_max, nS);

%% Loop through each file

for i = 1:length(data_paths)

    fp = fullfile(cfg.data_folder, data_paths(i).name);
    data = readtable(fp); 

    sid = sscanf(data_paths(i).name, cfg.sid_regex); 
    sid_list(i) = sid;

    % The last column is the response variable Y:
    Y_name = data.Properties.VariableNames{end};    
    Y = data{:, end};

    % All other columns are factors:
    F_names = data.Properties.VariableNames(1:end-1); 
    F_vals = data{:, 1:end-1}; 

    % Remove trials with extreme values (z-score > 3):
    if cfg.exclude_outlier
        rows_to_remove = any(abs(F_vals) > 3, 2);
        Y = Y(~ rows_to_remove, :);
        F_vals = F_vals(~ rows_to_remove, :);
    end

    % Normalization:
    if cfg.do_norm
        F_vals = normalize(F_vals, 'range', cfg.norm_range);
    end

    % Create the design matrix:
    nT = length(Y);       % number of trials
    nI = sum(nF-1:-1:1);  % number of interaction terms

    switch model_type
        case "CSR"
            X = [  
                ones(nT, 1),  ... % constant term
                F_vals,       ... % linear terms
                F_vals .^ 2,  ... % quadratic terms
                zeros(nT, nI) ... % pre-allocate memory for interaction terms
            ];  
            col = 1 + nF + nF; 
            for j = 1:nF-1
                for k = j+1:nF
                    col = col + 1;
                    X(:, col) = F_vals(:, j) .* F_vals(:, k); 
                end
            end

        case "GLM"
            X = [  
                ones(nT, 1),  ... % constant term
                F_vals,       ... % linear terms
            ];  
    end

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
    param_info = [sid, Coefs', Rsq, Rsq_adj, log_lik, IC.aic, IC.aicc, IC.bic, NRMSE];
%     param_info = [sid, Coefs', Rsq, Rsq_adj, log_lik, AIC, AICc, BIC, NRMSE];
    results_list(i, :) = param_info;

    % Fill the residuals into column i (replace NaN):
    residuals_list(1:length(residuals), i) = residuals; 
end

idx_Rsq = 1 + nP + 1;
fprintf('Average R-squared = %.3f\n', mean(results_list(:, idx_Rsq)));

%% Generate variable names and save to file

switch model_type
    case "CSR"
        var_names = [
            'SID', 'X0', ...
            F_names, ... % used to be strcat('F', string(1:nF))
            strcat('F', string(1:nF), '^2'), ...
            strcat('F', string(repelem(1:nF-1, nF-1:-1:1)), ... 
                   'F', string(cell2mat(arrayfun(@(x) (x:nF), 2:nF, 'UniformOutput', false)))), ... 
            model_measures ...
        ];

    case "GLM"
        var_names = ['SID', 'X0', F_names, model_measures];
end

results_list_table = array2table(results_list, 'VariableNames', var_names);
writetable(results_list_table, regression_results_out);

residuals_table = array2table(residuals_list, 'VariableNames', strcat("sub_", string(sid_list)));
writetable(residuals_table, residuals_table_out);

fprintf('Finish generating output:\n%s\n', regression_results_out)