function [publicKey, privateKey, runtime, numPrimalityTests, numIntervalsChecked] = generateRSANew(keySizeBits)

    fprintf('--- Starting Enhanced RSA Key Generation (using pi_new) ---\n');
    fprintf('Desired Key Size: %d bits\n', keySizeBits);

    tic;

    % Initializing counters and RSA parameters
    numPrimalityTests = 0;
    numIntervalsChecked = 0;
    p = sym(0);
    q = sym(0);
    e = sym(65537);
    primeBitLength = floor(keySizeBits / 2);

    %% Analyzing search space once using pi_new
    lowerBound = sym(2)^(primeBitLength - 1);
    upperBound = sym(2)^primeBitLength - 1;
    numSubIntervals = 100;
    intervalSize = (upperBound - lowerBound) / numSubIntervals;
    intervals = cell(numSubIntervals, 1);
    primeCounts = zeros(numSubIntervals, 1);
    
    fprintf('Analyzing search space in parallel with pi_new...\n');
    parfor i = 1:numSubIntervals
        start_i = lowerBound + (i-1) * intervalSize;
        end_i = lowerBound + i * intervalSize -1;
        if i == numSubIntervals, end_i = upperBound; end
        
        estimated_primes = floor(pi_new(double(end_i))) - floor(pi_new(double(start_i)));
        intervals{i} = [start_i, end_i];
        primeCounts(i) = max(0, estimated_primes);
    end

    [~, sortedIndices] = sort(primeCounts, 'descend');
    
    pool = gcp('nocreate');
    if isempty(pool), error('A parallel pool is required.'); end
    numCandidatesPerInterval = pool.NumWorkers * 128; % Increasing candidates

    % Generating small primes for trial division pre-screening
    smallPrimes = primes(3500); % Screening with primes up to 3500

    %% Searching for p and q in the rich intervals
    fprintf('Searching for two valid primes with trial division pre-screening...\n');
    
    for i = 1:numel(sortedIndices)
        currentIdx = sortedIndices(i);
        currentInterval = intervals{currentIdx};
        numIntervalsChecked = numIntervalsChecked + 1;
        
        candidates = rand_big_odd_int_in_range(currentInterval(1), currentInterval(2), numCandidatesPerInterval);
        
        validPrimes = [];
        testsInBatch = 0; % Counter for primality tests
        parfor k = 1:numel(candidates)
            c = candidates(k);

            % Trial Division
            isLikelyPrime = true;
            for j = 1:numel(smallPrimes)
                if mod(c, smallPrimes(j)) == 0
                    isLikelyPrime = false;
                    break;
                end
            end
            
            if isLikelyPrime
                testsInBatch = testsInBatch + 1;
                if isprime(c) && mod(c - 1, e) ~= 0
                    validPrimes = [validPrimes, c];
                end
            end
        end
        numPrimalityTests = numPrimalityTests + testsInBatch;
        
        for k = 1:numel(validPrimes)
            prime_k = validPrimes(k);
            if p == 0
                p = prime_k;
            elseif q == 0 && prime_k ~= p
                q = prime_k;
                break;
            end
        end
        
        if p > 0 && q > 0
            fprintf('Found valid primes in interval %d.\n', currentIdx);
            break;
        end
    end

    if p == 0 || q == 0
        error('Failed to find two valid primes. The search space may be too sparse or requires more candidates.');
    end

    fprintf('Found valid prime p: %s\n', char(p));
    fprintf('Found valid prime q: %s\n', char(q));

    %% Final RSA steps
    n = p * q;
    phi_n = (p - 1) * (q - 1);
    [~, d_inv, ~] = gcd(e, phi_n);
    d = mod(d_inv, phi_n);

    publicKey.n = n;
    publicKey.e = e;
    privateKey.n = n;
    privateKey.d = d;
    runtime = toc;
end

% Helper Functions
function r = rand_big_odd_int_in_range(lowerBound, upperBound, count)
    r = lowerBound + floor(rand(count, 1) .* (upperBound - lowerBound + 1));
    r(mod(r,2)==0) = r(mod(r,2)==0) + 1;
    r(r > upperBound) = [];
end