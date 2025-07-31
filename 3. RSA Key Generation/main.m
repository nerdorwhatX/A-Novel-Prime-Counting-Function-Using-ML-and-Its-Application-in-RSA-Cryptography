% This is my main script for comparing the performance of two RSA key
% generation methodologies: the standard one and my enhanced version that
% uses the pi_new function to guide the prime search.

% Setting up a clean environment.
clear all;
clc;
close all;

% This block starts the parallel pool once for the entire script.
% Also ensuring the pool is automatically closed when the script is done.
fprintf('Starting the parallel pool for the session...\n');
if isempty(gcp('nocreate'))
    parpool;
end
cleanupObj = onCleanup(@() delete(gcp('nocreate')));
fprintf('====================================================\n');

fprintf('  RSA Key Generation Methodology Comparison  \n');
fprintf('====================================================\n\n');

%% Getting User Input for Key Size
keySizeBits = input('Enter the desired RSA key size in bits (e.g., 1024, 2048, 4096): ');
if ~isscalar(keySizeBits) || ~isnumeric(keySizeBits) || keySizeBits <= 0 || mod(keySizeBits, 1) ~= 0
    error('Invalid input. Please enter a positive integer for key size.');
end
fprintf('\nStarting comparison for %d-bit RSA keys...\n\n', keySizeBits);

%% Running the Standard RSA Key Generation
fprintf('----------------------------------------------------\n');
fprintf('  Running Standard RSA Key Generation Method...     \n');
fprintf('----------------------------------------------------\n');
[pubKey_primary, privKey_primary, runtime_primary, numPrimalityTests_primary] = generateRSAOld(keySizeBits);

% Printing the keys for the Standard Method
fprintf('\n--- Standard Method Keys ---\n');
fprintf('Public Key (n, e):\n');
fprintf('  n = %s\n', char(pubKey_primary.n));
fprintf('  e = %s\n', char(pubKey_primary.e));
fprintf('Private Key (n, d):\n');
fprintf('  n = %s\n', char(privKey_primary.n));
fprintf('  d = %s\n', char(privKey_primary.d));
fprintf('----------------------------\n');

%% Running the Enhanced RSA Key Generation
fprintf('\n----------------------------------------------------\n');
fprintf('  Running Enhanced RSA Key Generation Method...     \n');
fprintf('----------------------------------------------------\n');
[pubKey_enhanced, privKey_enhanced, runtime_enhanced, numPrimalityTests_enhanced, numIntervalsChecked_enhanced] = generateRSANew(keySizeBits);

% Printing the keys for the Enhanced Method
fprintf('\n--- Enhanced Method Keys ---\n');
fprintf('Public Key (n, e):\n');
fprintf('  n = %s\n', char(pubKey_enhanced.n));
fprintf('  e = %s\n', char(pubKey_enhanced.e));
fprintf('Private Key (n, d):\n');
fprintf('  n = %s\n', char(privKey_enhanced.n));
fprintf('  d = %s\n', char(privKey_enhanced.d));
fprintf('----------------------------\n');

%% Comparing Results and Giving a Verdict
fprintf('\n====================================================\n');
fprintf('             Performance Comparison                 \n');
fprintf('====================================================\n');
fprintf('\nMetric                  | Standard Method | Enhanced Method\n');
fprintf('------------------------|-----------------|-----------------\n');
fprintf('Runtime (seconds)       | %15.4f | %15.4f\n', runtime_primary, runtime_enhanced);
fprintf('Primality Tests         | %15d | %15d\n', numPrimalityTests_primary, numPrimalityTests_enhanced);
% The 'Intervals Checked' metric is specific to the new method.
fprintf('Intervals Checked       | %15s | %15d\n', 'N/A', numIntervalsChecked_enhanced);
fprintf('\n');

%% Final Verdict
fprintf('====================================================\n');
fprintf('                  Final Verdict                     \n');
fprintf('====================================================\n');

% Calculating performance differences
runtime_diff = runtime_primary - runtime_enhanced;
tests_diff = numPrimalityTests_primary - numPrimalityTests_enhanced;

% Providing a reasoned verdict based on performance metrics
if runtime_enhanced < runtime_primary
    % Enhanced method is faster
    verdict_winner = 'ENHANCED';
    runtime_improvement_percent = (runtime_diff / runtime_primary) * 100;
    
    fprintf('Verdict: The ENHANCED method is demonstrably SUPERIOR for %d-bit keys.\n\n', keySizeBits);
    fprintf('Reasoning:\n');
    fprintf(' - It was %.4f seconds faster, a significant %.2f%% performance improvement.\n', runtime_diff, runtime_improvement_percent);
    
    if tests_diff > 0
        fprintf(' - It was also more efficient, requiring %d fewer primality tests.\n', tests_diff);
    elseif tests_diff == 0
        fprintf(' - It required the same number of primality tests.\n');
    else
        fprintf(' - It surprisingly required %d more primality tests, but the faster runtime shows\n', -tests_diff);
        fprintf('   its overall strategy was more effective.\n');
    end
    fprintf('This indicates a more efficient prime searching strategy, leading to faster key generation.\n');

elseif runtime_primary < runtime_enhanced
    % Standard method is faster
    verdict_winner = 'STANDARD';
    
    fprintf('Verdict: The STANDARD method performed better in this instance for %d-bit keys.\n\n', keySizeBits);
    fprintf('Reasoning:\n');
    fprintf(' - It was %.4f seconds faster.\n', -runtime_diff);
    
    if tests_diff < 0
        fprintf(' - It also required %d fewer primality tests.\n', -tests_diff);
    else
        fprintf(' - Although it required %d more primality tests, its overall runtime was better.\n', tests_diff);
    end
    fprintf('This is an unexpected outcome, suggesting the overhead of the enhanced method''s logic\n');
    fprintf('did not pay off for this specific key size or run.\n');
    
else
    % Runtimes are virtually identical
    verdict_winner = 'INCONCLUSIVE';
    
    fprintf('Verdict: The performance was nearly IDENTICAL.\n\n');
    fprintf('Reasoning:\n');
    fprintf(' - The runtime difference was negligible (less than a millisecond).\n');
    if tests_diff > 0
        fprintf(' - The ENHANCED method has a slight edge in efficiency, requiring %d fewer primality tests.\n', tests_diff);
    elseif tests_diff < 0
        fprintf(' - The STANDARD method has a slight edge in efficiency, requiring %d fewer primality tests.\n', -tests_diff);
    else
        fprintf(' - All key metrics (runtime and tests) were identical. Both methods are equally performant.\n');
    end
end

fprintf('\nComparison complete.\n');
fprintf('====================================================\n');