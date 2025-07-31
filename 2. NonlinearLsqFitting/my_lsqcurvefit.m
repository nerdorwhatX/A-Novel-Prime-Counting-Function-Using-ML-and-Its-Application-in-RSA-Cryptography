function [x,resnorm,residual,exitflag,output,lambda,jacobian] = my_lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
% My from-scratch implementation of the Levenberg-Marquardt algorithm.
% This function is designed to solve non-linear least squares problems 
% It finds the parameters 'x' that minimize the sum of squared differences 
% between a function 'fun(x, xdata)' and observed 'ydata'.

%% Handling arguments and setting default options.
if nargin < 7
    options = [];
end
if nargin < 6
    ub = [];
end
if nargin < 5
    lb = [];
end

%% Ensuring x0, lb, and ub are column vectors for consistency.
x0 = x0(:);
if ~isempty(lb)
    lb = lb(:);
else
    lb = -inf(size(x0));
end
if ~isempty(ub)
    ub = ub(:);
else
    ub = inf(size(x0));
end

%% Checking for valid inputs.
if any(lb > ub)
    error('Lower bounds must be less than or equal to upper bounds.');
end
if any(x0 < lb) || any(x0 > ub)
    error('Initial point x0 is not within the specified bounds.');
end

%% Getting options from the options structure, or using my defaults.
tolX = optimget(options, 'TolX', 1e-6);
tolFun = optimget(options, 'TolFun', 1e-6);
maxIter = optimget(options, 'MaxIterations', 400);
display = optimget(options, 'Display', 'final');

%% Initializing algorithm parameters.
x = x0; % Starting with the initial guess.
n = length(x);

%% Setting up Levenberg-Marquardt specific parameters.
lambda_val = 1e-2; % This is the initial damping parameter.
lambda_update_factor = 10; % Factor for increasing/decreasing damping.

%% Initializing counters and flags for the loop.
iter = 0;
exitflag = 0;
stop = false;

%% Setting up the header for iterative display.
if strcmp(display, 'iter')
    fprintf('\n                                         Norm of         First-order \n');
    fprintf(' Iteration  Func-count     f(x)          step          optimality\n');
end

%% Main Levenberg-Marquardt iteration loop.
while ~stop
    iter = iter + 1;

    % Evaluating the function and residual at the current point.
    F = fun(x, xdata);
    F = F(:); % Ensuring F is a column vector.
    residual = F - ydata(:);
    resnorm = sum(residual.^2);

    % Calculating the Jacobian matrix (J) using finite differences.
    J = zeros(length(ydata), n);
    h = 1e-6; % Step size for the finite difference approximation.
    for j = 1:n
        x_plus_h = x;
        x_plus_h(j) = x_plus_h(j) + h;
        F_plus_h = fun(x_plus_h, xdata);
        J(:, j) = (F_plus_h(:) - F) / h;
    end
    
    % Calculating the gradient (g) and approximate Hessian (H).
    g = J' * residual;
    H = J' * J;

    funcCount = 1 + n; % Tracking function evaluations (1 for F, n for Jacobian).

    % Solving for the next step.
    step_found = false;
    while ~step_found && iter <= maxIter
        % Applying the Levenberg-Marquardt damping to the Hessian.
        H_lm = H + lambda_val * eye(n);

        % Solving the linear system to find the step direction.
        try
            step = -H_lm \ g;
        catch
            % If H_lm is singular, I'm increasing damping and trying again.
            lambda_val = lambda_val * lambda_update_factor;
            continue;
        end

        % Checking if the proposed step is acceptable.
        x_new = x + step;

        % Projecting the new point back onto the bounds.
        x_new = max(lb, min(ub, x_new));

        % Evaluating the function at the potential new point.
        F_new = fun(x_new, xdata);
        residual_new = F_new(:) - ydata(:);
        resnorm_new = sum(residual_new.^2);
        funcCount = funcCount + 1;

        % Checking for improvement in the residual.
        if resnorm_new < resnorm
            % Success! Accepting the new point and decreasing damping.
            x = x_new;
            lambda_val = lambda_val / lambda_update_factor;
            step_found = true;
        else
            % Failure. Rejecting the step and increasing damping.
            lambda_val = lambda_val * lambda_update_factor;
        end

        % Checking for convergence.
        if norm(step) < tolX
            exitflag = 2;
            stop = true;
            break;
        end
        if abs(resnorm_new - resnorm) < tolFun
            exitflag = 3;
            stop = true;
            break;
        end
        if lambda_val > 1e16 % Heuristic for a "stuck" condition.
             exitflag = -3;
             stop = true;
             break;
        end
    end

    % Displaying iteration information if requested.
    if strcmp(display, 'iter')
        first_order_optimality = norm(g, inf);
        fprintf('%5.0f      %5.0f   %13.6g   %13.6g   %12.3g\n', iter, funcCount, resnorm, norm(step), first_order_optimality);
    end
    
    % Checking if we've hit the maximum number of iterations.
    if iter >= maxIter
        exitflag = 0;
        stop = true;
    end

end

%% Finalizing the outputs.
[residual, jacobian] = deal(residual, J); % Using last computed values.
resnorm = sum(residual.^2);

%% Building the output structure.
output.iterations = iter;
output.funcCount = funcCount;
output.algorithm = 'Levenberg-Marquardt';
output.firstorderopt = norm(g, inf);
output.message = getExitMessage(exitflag, maxIter, tolFun, tolX);

%% Building the lambda structure for Lagrange multipliers.
% For box constraints, these relate to the gradient at active bounds.
lambda.lower = zeros(n,1);
lambda.upper = zeros(n,1);
active_lower = (x <= lb + 1e-6);
active_upper = (x >= ub - 1e-6);
lambda.lower(active_lower) = g(active_lower);
lambda.upper(active_upper) = -g(active_upper);


if strcmp(display, 'final') || strcmp(display, 'iter')
    fprintf('\n%s\n', output.message);
end

end

%% Helper Functions

function val = optimget(options, name, default)
% A simplified helper function to get values from the options struct.
    if ~isempty(options) && isfield(options, name)
        val = options.(name);
    else
        val = default;
    end
end

function msg = getExitMessage(flag, maxIter, tolFun, tolX)
% A helper to create the exit message based on the final exitflag.
    switch flag
        case 1
            msg = 'Function converged to a solution.';
        case 2
            msg = sprintf('Change in x was less than the specified tolerance TolX = %g.', tolX);
        case 3
            msg = sprintf('Change in the residual was less than the specified tolerance TolFun = %g.', tolFun);
        case 0
            msg = sprintf('Number of iterations exceeded options.MaxIterations = %d.', maxIter);
        case -1
            msg = 'Algorithm was terminated by the user.';
        case -2
            msg = 'The problem is infeasible.';
        case -3
            msg = 'The trust region radius became too small.';
        otherwise
            msg = 'Optimization ended unexpectedly.';
    end
end