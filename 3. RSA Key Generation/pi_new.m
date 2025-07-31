function y = pi_new(x)
% Calculates the prime-counting function approximation using my new fitted model.
% This function uses a pre-calculated set of parameters that I've stored
% in 'fitted_parameters.mat'. It expects that file to be in the path.
% The model formula I'm using is:
% pi(x) approx x / (ln(x) - p1 - p2/ln(x) - p3/ln(x)^2)

% Using a persistent variable for the parameters. This way, I only have to
% load the .mat file once per session, which is much faster for repeated calls.
persistent fitted_params;

% Checking if the parameters are loaded yet. If not, I'm loading them now.
if isempty(fitted_params)
    param_file = 'fitted_parameters.mat';
    if ~isfile(param_file)
        error('Parameter file "%s" not found. Please run the fitting script first to generate it.', param_file);
    end
    loaded_data = load(param_file, 'p_fit');
    fitted_params = loaded_data.p_fit;
    fprintf('Loaded fitted parameters for pi_new() function.\n');
end

% Ensuring the input x is valid, since the model is undefined for x < 2.
if any(x < 2)
    warning('Input x contains values less than 2. The model is undefined for this range. Results for these values will be NaN.');
    x(x < 2) = NaN; % Setting invalid inputs to NaN to prevent errors in the formula.
end

% Extracting the parameters to make the formula below more readable.
p1 = fitted_params(1);
p2 = fitted_params(2);
p3 = fitted_params(3);

% Calculating the approximation using the vectorized formula.
log_x = log(x);
denominator = log_x - p1 - p2./log_x - p3./(log_x.^2);
y = x ./ denominator;

end