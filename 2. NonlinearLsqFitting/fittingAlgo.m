% High-accuracy nonlinear least squares regression for the prime-counting
% function, pi(x).
% This script fits a model using a representative subset of the data to
% manage memory usage, then validates the RMSE on the complete dataset
% via chunked processing. It also visualizes the results and compares the
% model's accuracy against standard approximations.

%% Environment Setup
clear;
clc;
close all;
fprintf('Starting Prime-Counting Function Regression...\n');

% Configuration
data_filename = 'prime_counting_dataset_gpu.csv';

% Memory management parameters
num_fit_points = 1000000;  % Points for the fitting process (balances accuracy and RAM usage)
num_plot_points = 2000;    % Points for visual clarity in plots
rmse_chunk_size = 100000;  % Chunk size for final RMSE calculation on the full dataset

%% Loading the Prime-Counting Dataset
fprintf('Loading dataset: %s\n', data_filename);
if ~isfile(data_filename)
    error('Dataset file not found: %s', data_filename);
end

try
    T = readtable(data_filename);
    x_data_full = double(T.X(:));
    y_data_full = double(T.Pi_X(:));

    % Filtering data, as the model is undefined for x < 2.
    valid_indices = x_data_full >= 2;
    x_data_full = x_data_full(valid_indices);
    y_data_full = y_data_full(valid_indices);

    total_points = length(x_data_full);
    fprintf('Dataset loaded successfully with %d total data points.\n\n', total_points);
catch ME
    fprintf('Error loading or processing the CSV file.\n');
    rethrow(ME);
end

%% Defining the Regression Model
% The model is inspired by the asymptotic expansion of Li(x).
model_func = @(p, x) x ./ (log(x) - p(1) - p(2)./log(x) - p(3)./(log(x).^2));
p0 = [1.0, 1.0, 2.0]; % Initial guess for [p1, p2, p3]

fprintf('Model defined: pi(x) = x / (log(x) - p1 - p2/log(x) - p3/(log(x)^2))\n');
fprintf('Initial parameter guess p0 = [%.2f, %.2f, %.2f]\n\n', p0(1), p0(2), p0(3));

%% Performing Nonlinear Least Squares Regression
% Fitting on a subset of the data to prevent out-of-memory errors.
fprintf('Creating a fitting subset of %d points to prevent memory errors.\n', num_fit_points);
fit_indices = round(linspace(1, total_points, min(num_fit_points, total_points)));
x_fit_subset = x_data_full(fit_indices);
y_fit_subset = y_data_full(fit_indices);

fprintf('Setting up lsqcurvefit to run on the subset (%d points).\n', length(x_fit_subset));

% Setting optimization options
options = optimoptions('lsqcurvefit', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'iter-detailed', ...
    'MaxIterations', 100, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'UseParallel', true);

% Running the regression on the subset
tic;
[p_fit, resnorm, ~, exitflag] = my_lsqcurvefit(model_func, p0, x_fit_subset, y_fit_subset, [], [], options);
elapsed_time = toc;

fprintf('\n--- Regression Finished ---\n');
fprintf('Elapsed Time: %.2f seconds\n', elapsed_time);
fprintf('Exit Flag: %d (1 indicates success)\n', exitflag);
fprintf('Final Sum of Squared Residuals (on subset): %e\n', resnorm);
fprintf('\n--- Optimal Parameters Found ---\n');
fprintf('p1 = %.12f\n', p_fit(1));
fprintf('p2 = %.12f\n', p_fit(2));
fprintf('p3 = %.12f\n\n', p_fit(3));

% Printing the final fitted model equation
fprintf('--- Final Fitted Model Equation ---\n');
p1_str = sprintf('%.6f', p_fit(1));
% Handling signs for p2 and p3 for clean printing
if p_fit(2) >= 0
    p2_op = '-'; p2_val = p_fit(2);
else
    p2_op = '+'; p2_val = -p_fit(2);
end
if p_fit(3) >= 0
    p3_op = '-'; p3_val = p_fit(3);
else
    p3_op = '+'; p3_val = -p_fit(3);
end
fprintf('pi(x) approx x / (ln(x) - %s %s %.6f/ln(x) %s %.6f/ln(x)^2)\n\n', ...
    p1_str, p2_op, p2_val, p3_op, p3_val);

%% Saving Fitted Parameters
% Saving the p_fit vector to a .mat file for use in other functions.
fprintf('Saving fitted parameters to fitted_parameters.mat...\n');
save('fitted_parameters.mat', 'p_fit');
fprintf('Parameters saved successfully.\n\n');

%% Analyzing and Visualizing Results
fprintf('Analyzing the results...\n');

% Checking for benchmark function files.
if ~isfile('Li.m') || ~isfile('R.m')
    error('Benchmark functions Li.m and/or R.m not found.');
end

% Using a smaller, evenly spaced subset for clear and fast plotting.
plot_indices = round(linspace(1, total_points, num_plot_points));
x_plot = x_data_full(plot_indices);
y_plot = y_data_full(plot_indices);

% Calculating fitted model and benchmark approximations for the plot points.
fprintf('Calculating model values on %d points for plotting...\n', num_plot_points);
y_fit_plot = model_func(p_fit, x_plot);
y_li_plot = Li(x_plot);
y_r_plot = R(x_plot);

% Plot 1: Overall Fit (Log-Log Scale)
figure('Name', 'Model Fit vs. Data', 'NumberTitle', 'off');
loglog(x_data_full, y_data_full, 'k.', 'DisplayName', 'True pi(x) Data', 'MarkerSize', 4);
hold on;
% Plotting the fitted line using the smaller 'x_plot' subset to save memory
loglog(x_plot, y_fit_plot, 'r-', 'DisplayName', 'Fitted Model', 'LineWidth', 2);
grid on;
xlabel('x (log scale)');
ylabel('pi(x) (log scale)');
title('Fitted Model vs. True pi(x)');
legend('show', 'Location', 'SouthEast');

% Plot 2: Residual Error Comparison
figure('Name', 'Residual Error Analysis', 'NumberTitle', 'off');
error_fit_plot = y_plot - y_fit_plot;
error_li_plot = y_plot - y_li_plot;
error_r_plot = y_plot - y_r_plot;
semilogx(x_plot, error_fit_plot, 'r-', 'DisplayName', 'Error of Our Fitted Model', 'LineWidth', 2.0);
hold on;
semilogx(x_plot, error_r_plot, 'g-.', 'DisplayName', 'Error of R(x)', 'LineWidth', 1.5);
semilogx(x_plot, error_li_plot, 'b--', 'DisplayName', 'Error of Li(x)', 'LineWidth', 1.5);
grid on;
ax = gca;
ax.YAxis.Exponent = 0;
xlabel('x (log scale)');
ylabel('Error (True pi(x) - Approximation)');
title('Comparison of Residual Errors');
legend('show', 'Location', 'NorthWest');

% Plot 3: Intervals of Model Superiority
figure('Name', 'Model Superiority Intervals', 'NumberTitle', 'off');

% Calculating absolute errors for each model
abs_error_fit = abs(error_fit_plot);
abs_error_li = abs(error_li_plot);
abs_error_r = abs(error_r_plot);

% Finding which model has the minimum error for each point
% 1: Fitted Model, 2: R(x), 3: Li(x)
[~, best_model_idx] = min([abs_error_fit, abs_error_r, abs_error_li], [], 2);

% Defining colors and markers for clarity
colors = {'r', 'g', 'b'};
markers = {'s', 'o', '^'};
model_names = {'Fitted Model is Best', 'R(x) is Best', 'Li(x) is Best'};

hold on;
for i = 1:3
    % Finding all x-points where this model was the best
    wins = (best_model_idx == i);
    % Plotting these points as a distinct series
    semilogx(x_plot(wins), i * ones(sum(wins), 1), ...
        'Marker', markers{i}, ...
        'Color', colors{i}, ...
        'MarkerFaceColor', colors{i}, ...
        'LineStyle', 'none', ...
        'DisplayName', model_names{i});
end
hold off;

% Formatting the plot
set(gca, ...
    'YTick', 1:3, ...
    'YTickLabel', {'Fitted Model', 'R(x)', 'Li(x)'}, ...
    'YLim', [0.5, 3.5]);
grid on;
xlabel('x (log scale)');
ylabel('Best Performing Model');
title('Intervals of Model Superiority (Lowest Absolute Error)');
legend('show', 'Location', 'best');


fprintf('Plotting complete.\n');

%% Final Quantitative Evaluation
% Calculating RMSE on the FULL dataset using chunks to avoid memory errors.
fprintf('\nCalculating final RMSE over the entire dataset (%d points) in chunks...\n', total_points);

sum_sq_err_fit = 0;
sum_sq_err_li = 0;
sum_sq_err_r = 0;

num_chunks = ceil(total_points / rmse_chunk_size);
for i = 1:num_chunks
    start_idx = (i-1) * rmse_chunk_size + 1;
    end_idx = min(i * rmse_chunk_size, total_points);

    x_chunk = x_data_full(start_idx:end_idx);
    y_chunk = y_data_full(start_idx:end_idx);

    % Calculating errors for the current chunk
    err_fit_chunk = y_chunk - model_func(p_fit, x_chunk);
    err_li_chunk  = y_chunk - Li(x_chunk);
    err_r_chunk   = y_chunk - R(x_chunk);

    % Accumulating sum of squares
    sum_sq_err_fit = sum_sq_err_fit + sum(err_fit_chunk.^2);
    sum_sq_err_li  = sum_sq_err_li  + sum(err_li_chunk.^2);
    sum_sq_err_r   = sum_sq_err_r   + sum(err_r_chunk.^2);

    fprintf('Processed chunk %d/%d...\n', i, num_chunks);
end

% Calculating final RMSE
rmse_fit = sqrt(sum_sq_err_fit / total_points);
rmse_li  = sqrt(sum_sq_err_li / total_points);
rmse_r   = sqrt(sum_sq_err_r / total_points);

fprintf('\n--- Root Mean Squared Error (RMSE) Comparison (Full Dataset) ---\n');
fprintf('RMSE of Li(x):      %.4f\n', rmse_li);
fprintf('RMSE of R(x):       %.4f\n', rmse_r);
fprintf('RMSE of Fitted Model: %.4f\n', rmse_fit);