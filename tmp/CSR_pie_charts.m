clear; clc; close all

%% Load data table

authors = ["Chang_et_al", "Tse_et_al"]; 
tasks = ["Naming", "LD"]; 
data_folder = fullfile("Output_Linguistic", authors(1), tasks(1)); 
fig_folder = fullfile("Figs_Linguistic", authors(1), tasks(1), "Pies"); 

out_tags = ["", " (z-scored)"];
data_name = join(['CSR_regression_results', out_tags(2), ".csv"], '');

data = readtable(fullfile(data_folder, data_name)); 
F_names = ['LogCF' 'NS' 'REG' 'UNP' 'CON' 'PC' 'SC' 'SAR' 'IMG' 'AoA'];

%% Define the number and structure of parameters

nF = 8;                     % number of factors
nP = (nF^2 + 3*nF + 2) / 2; % number of parameters
nS = size(data, 1);         % number of subjects

%%
for i = 1:nS

    sid = data(i, 1); 
    const_value = data(i, 2);              % constant term
    linear_values = data(i, 2+1:2+nF);     % linear terms
    quad_values = data(i, 2+nF+1:2+2*nF);  % quadratic terms
    inter_values = data(i, 2+2*nF+1:nP+1); % interaction terms  

    % Draw pie chart
    figure('Name', sprintf('Subject %d', sid), 'Position', [100, 100, 1200, 600]);

    % 1. Pie chart showing the proportions of the four major categories
    category_values = [
        const_value, ...
        linear_values, ... % sum(linear_vals)
        quad_values, ...   % sum(quad_vals)
        sum(inter_values)];
    category_labels = [
        {'Constant'}, ...
        strcat(F_names, '-L'), ... % 'Linear'
        strcat(F_names, '-Q'), ... % 'Quadratic'
        {'Interaction'}];
    colors1 = [
        % light gray:
        0.8, 0.8, 0.8; 
        % intermediate colors:
        linspace(0.2, 0.6, nF)', linspace(0.5, 0.8, nF)', linspace(0.8, 1.0, nF)'; 
        % dark colors:
        linspace(0, 0.4, nF)', linspace(0.4, 0.7, nF)', linspace(0.6, 0.9, nF)';
         % light colors:
        repmat([0.9, 0.9, 0.6], length(inter_idx), 1)
    ];
    subplot(1, 2, 1);
    pie(category_values, category_labels);
    colormap(gca, colors);
    title(sprintf('Subject %d - Category Proportions', sid));

    % 2. Pie chart showing the proportions of each parameter
    param_values = [
        const_value, ...
        linear_values, ...
        quad_values, ...
        mean(coefs(inter_idx))];
    param_labels = [ ...
        {'Constant'}, ...
        arrayfun(@(x) sprintf('P%d (Linear)', x), 1:nF, 'UniformOutput', false), ...
        arrayfun(@(x) sprintf('P%d (Quadratic)', x), 1:nF, 'UniformOutput', false), ...
        'Interaction Mean'];
    subplot(1, 2, 2);
    pie(param_values, param_labels);
    colormap(gca, colors);
    title(sprintf('Subject %d - Parameter Proportions', sid));

    % Save image
    saveas(gcf, fullfile(fig_folder, sprintf('Subject_%d_PieCharts.png', sid)));
end
