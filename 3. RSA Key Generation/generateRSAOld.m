function [publicKey, privateKey, runtime, numPrimalityTests] = generateRSAOld(keySizeBits)

    fprintf('--- Starting Standard RSA Key Generation ---\n');
    fprintf('Desired Key Size: %d bits\n', keySizeBits);
    
    tic;
    
    % Initializing counters and RSA parameters
    numPrimalityTests = 0;
    e = sym(65537);
    primeBitLength = floor(keySizeBits / 2);
    
    pool = gcp('nocreate');
    if isempty(pool), error('A parallel pool is required.'); end
    numCandidatesToTest = pool.NumWorkers * 128; % Increasing candidates per batch

    % Generating a list of small primes for trial division pre-screening.
    % Doing this once and broadcasting to the parallel workers.
    smallPrimes = primes(3500); % Screening with primes up to 3500

    %% Finding two distinct large primes, p and q, in a single parallel search.
    fprintf('Searching for two valid primes (p and q) with trial division pre-screening...\n');
    
    foundPrimes = sym([]);
    while numel(foundPrimes) < 2
        candidates = rand_big_odd_int(primeBitLength, numCandidatesToTest);
        
        % Using a temporary variable for reduction in the parfor loop
        validPrimesBatch = [];
        testsInBatch = 0; % Counter for primality tests in this batch
        
        parfor i = 1:numCandidatesToTest
            c = candidates(i);
            
            % Trial Division 
            % Checking for divisibility by small primes first. This is very fast.
            isLikelyPrime = true;
            for k = 1:numel(smallPrimes)
                if mod(c, smallPrimes(k)) == 0
                    isLikelyPrime = false;
                    break;
                end
            end
            
            % Only calling expensive isprime() if it passes the initial screen.
            if isLikelyPrime
                testsInBatch = testsInBatch + 1; % Incrementing test counter
                if isprime(c) && mod(c - 1, e) ~= 0
                    validPrimesBatch = [validPrimesBatch, c];
                end
            end
        end
        
        numPrimalityTests = numPrimalityTests + testsInBatch;
        
        if ~isempty(validPrimesBatch)
            % Adding unique new primes to our list of found primes
            foundPrimes = unique([foundPrimes, validPrimesBatch]);
        end
    end
    
    p = foundPrimes(1);
    q = foundPrimes(2);
    
    fprintf('Found prime p: %s (bits: %d)\n', char(p), bitlength(p));
    fprintf('Found prime q: %s (bits: %d)\n', char(q), bitlength(q));

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
function r = rand_big_odd_int(numBits, count)
    lowerBound = sym(2)^(numBits - 1);
    upperBound = sym(2)^numBits - 1;
    r = lowerBound + floor(rand(count, 1) .* (upperBound - lowerBound + 1));
    r(mod(r,2)==0) = r(mod(r,2)==0) + 1;
end

function bits = bitlength(num)
    if num == 0, bits = 0; return; end
    if isa(num, 'sym'), bits = floor(log2(num)) + 1;
    else, bits = floor(log2(double(num))) + 1;
    end
end