function y = R(x)
% Computes the Riemann R-function, a highly accurate benchmark for pi(x).
% I'm using the Gram Series formulation for the calculation:
% R(x) = 1 + sum_{k=1 to inf} [ (ln(x)^k) / (k * k! * zeta(k+1)) ]
% The series is truncated since it converges very quickly.

    % Checking for valid input.
    if any(x < 1)
        error('Input x must be greater than or equal to 1.');
    end
    
    % Setting the number of terms. The series converges so fast that 50
    % terms is more than enough for double precision.
    N_TERMS = 50;
    
    % Getting the natural log of x now to avoid recalculating it in the loop.
    log_x = log(x);
    
    % Initializing the output vector. The series starts with 1.
    y = ones(size(x));
    
    % Setting up variables for the loop.
    % I'm updating these iteratively to avoid recomputing large numbers.
    logx_k = ones(size(x)); % Stores (ln(x))^k
    k_factorial = 1;       % Stores k!
    
    % Looping to compute the sum of the series terms.
    for k = 1:N_TERMS
        % Updating the power of log(x) for the current k.
        logx_k = logx_k .* log_x;
        
        % Updating the factorial for the current k.
        k_factorial = k_factorial * k;
        
        % Calculating the k-th term of the Gram series.
        term = logx_k / (k * k_factorial * zeta(k + 1));
        
        % Adding the current term to the running sum.
        y = y + term;
    end

end