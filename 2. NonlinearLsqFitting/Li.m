function y = Li(x)
% Computes the offset logarithmic integral, Li(x), a standard benchmark
% for the prime-counting function pi(x).
% Li(x) is the integral of 1/ln(t) from 2 to x. My approach here is to
% calculate it efficiently using the relation Li(x) = li(x) - li(2),
% where li(x) is MATLAB's built-in logarithmic integral.

    % Checking for valid input.
    if any(x < 0)
        error('Input x must be non-negative.');
    end
    
    % Defining the constant value for li(2).
    li2 = 1.045163780117492; 
    
    % Calculating Li(x) using the fast, built-in logint function.
    y = logint(x) - li2;

end