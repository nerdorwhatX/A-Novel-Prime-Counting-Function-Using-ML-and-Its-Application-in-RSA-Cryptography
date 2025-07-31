% MATLAB Code for Generating Pi(x) Dataset using Segmented Sieve
% Optimized for large N by collecting a sparse set of (x, pi(x)) pairs.

% Main Script for Demonstration

%% Clearing workspace and command window
clear;
clc;

%% Checking for GPU and Parallel Computing Toolbox
fprintf('Checking for GPU and Parallel Computing Toolbox...\n');
if ~license('test', 'Distrib_Computing_Toolbox')
    error('Parallel Computing Toolbox is not available. This script requires it to run on a GPU.');
end
if gpuDeviceCount == 0
    error('No compatible GPU found. Please ensure you have a supported NVIDIA GPU and updated drivers.');
else
    gpu = gpuDevice(1); % Getting the first GPU device
    fprintf('GPU Found: %s (Compute Capability: %s)\n', gpu.Name, gpu.ComputeCapability);
    fprintf('Starting Segmented Sieve for Pi(x) calculation on the GPU...\n\n');
end


%% Defining the minimum value of x (min_N) for which to collect pi(x).
min_N = 1; % Starting collecting from 1

%% User Input for Maximum N
max_N_input = input('Enter the desired maximum X value (e.g., 1e8, 1e10, 1e12): ');
if ~isscalar(max_N_input) || max_N_input < 1
    error('Maximum X value must be a positive number.');
end
max_N = max_N_input;

%% Ensuring min_N is not greater than max_N or negative
if min_N > max_N
    error('min_N cannot be greater than max_N.');
end
if min_N < 0
    error('min_N cannot be negative.');
end

%% User Input for Number of Data Points
num_points = input('Enter the desired number of data points (e.g., 100, 500, 2000): ');
if ~isscalar(num_points) || num_points < 1 || floor(num_points) ~= num_points
    error('Number of data points must be a positive integer.');
end

%% Generating logarithmically spaced x values to collect.
if max_N > 1
    x_values_to_collect_raw = unique(round(logspace(0, log10(max_N), num_points)));
else
    x_values_to_collect_raw = [0, 1];
end

%% Filtering x_values_to_collect to be within the [min_N, max_N] range.
x_values_to_collect = x_values_to_collect_raw(x_values_to_collect_raw >= min_N & x_values_to_collect_raw <= max_N);
if min_N == 0 && ~ismember(0, x_values_to_collect)
    x_values_to_collect = [0, x_values_to_collect];
end
x_values_to_collect = unique(sort(x_values_to_collect));

%% Defining the size of each segment
segment_size = 2^24; % 16.7 million

fprintf('Target range: x from %e to %e\n', min_N, max_N);
fprintf('Number of data points to collect: %d\n', length(x_values_to_collect));
fprintf('Segment Size: %e\n', segment_size);

tic; % Starting timer

%% Calling the segmented sieve function to get pi(x) for the specified x values.
[x_data_final, pi_x_for_dataset] = segmentedSievePi_GPU(max_N, segment_size, x_values_to_collect);

toc; % Ending timer

% Creating a table to store the results
ground_truth_dataset = table(x_data_final', pi_x_for_dataset', 'VariableNames', {'X', 'Pi_X'});

% Saving this dataset to a CSV file.
csv_filename = 'prime_counting_dataset_gpu.csv';
writetable(ground_truth_dataset, csv_filename, 'Delimiter', ',', 'WriteRowNames', false);
fprintf('\nDataset saved to %s\n', csv_filename);

%% Function Definitions

function primes = sievePrimes_GPU(limit)
% SIEVEPRIMES_GPU Implements Sieve of Eratosthenes on the GPU.
% primes = SIEVEPRIMES_GPU(limit) returns a row vector of primes up to 'limit'.
% The result is returned to the CPU memory.

    if limit < 2
        primes = [];
        return;
    end

    % Creating a logical array directly on the GPU
    is_prime = gpuArray(true(1, limit));
    is_prime(1) = false; % 1 is not a prime number

    % The outer loop running on the CPU, but operations inside affecting the GPU array.
    for p = 2:sqrt(limit)
        % This check reading a single value from the GPU, which is fast.
        if is_prime(p)
            % The large-scale marking operation happening on the GPU
            is_prime(p*p:p:limit) = false;
        end
    end

    % 'find' executing on the GPU, and 'gather' moving the result back to CPU memory
    primes = gather(find(is_prime));
end


function [x_data_collected, pi_x_data_collected] = segmentedSievePi_GPU(N, segmentSize, x_values_to_collect)
% SEGMENTEDSIEVEPI_GPU Calculates pi(x) for specified x values using a GPU-accelerated Segmented Sieve.
% N: The upper limit for prime counting.
% segmentSize: The size of each segment to process.
% x_values_to_collect: A sorted array of x values for which to store pi(x).

    if N < 1
        x_data_collected = [];
        pi_x_data_collected = [];
        return;
    end

    x_values_to_collect = unique(sort(x_values_to_collect(x_values_to_collect >= 0 & x_values_to_collect <= N)));

    if isempty(x_values_to_collect)
        warning('No valid x values to collect within the range [0, N].');
        x_data_collected = [];
        pi_x_data_collected = [];
        return;
    end

    % Pre-allocating arrays for results on the CPU
    x_data_collected = zeros(1, length(x_values_to_collect));
    pi_x_data_collected = zeros(1, length(x_values_to_collect));
    collection_idx = 1;
    target_x_ptr = 1;

    % Handling pi(0) and pi(1) explicitly
    while target_x_ptr <= length(x_values_to_collect) && x_values_to_collect(target_x_ptr) <= 1
        current_x = x_values_to_collect(target_x_ptr);
        if current_x == 0
            x_data_collected(collection_idx) = 0;
            pi_x_data_collected(collection_idx) = 0;
        elseif current_x == 1
            x_data_collected(collection_idx) = 1;
            pi_x_data_collected(collection_idx) = 0;
        end
        collection_idx = collection_idx + 1;
        target_x_ptr = target_x_ptr + 1;
    end

    % Generating primes up to sqrt(N) using the GPU-accelerated sieve
    limit_sqrt_N = floor(sqrt(N));
    fprintf('Generating base primes up to %d...\n', limit_sqrt_N);
    small_primes = sievePrimes_GPU(limit_sqrt_N);
    fprintf('Base primes generated. Starting segmented sieve...\n');

    % Moving the small primes list to the GPU for faster access in the loop
    small_primes_gpu = gpuArray(small_primes);

    current_total_prime_count = 0;

    % Iterating through segments
    for low = 2:segmentSize:N
        high = min(low + segmentSize - 1, N);

        % Progress Bar: Displaying the current progress
        % The '\r' at the end moving the cursor to the start of the line,
        % creating a dynamic progress indicator.
        percent_done = (high / N) * 100;
        fprintf('Processing... %.2f%% complete.\r', percent_done);

        % Creating the segment's boolean array on the GPU
        is_segment_prime_gpu = gpuArray(true(1, high - low + 1));

        % Using the small_primes on the GPU to mark composites within the current segment
        for p_val = small_primes_gpu
            p = double(p_val); % Casting to double for calculations
            start_val = max(p*p, ceil(low/p)*p);
            
            indices_to_mark = start_val:p:high;
            if ~isempty(indices_to_mark)
                is_segment_prime_gpu(indices_to_mark - low + 1) = false;
            end
        end
        
        % Gathering the boolean array from GPU to CPU for data collection
        segment_primes_cpu = gather(is_segment_prime_gpu);
        
        % Iterating through numbers in the now-local segment to update the total prime count
        % and collecting data points.
        for i = 1:length(segment_primes_cpu)
            current_num = low + i - 1;
            if current_num > N, break; end

            if segment_primes_cpu(i)
                current_total_prime_count = current_total_prime_count + 1;
            end

            % Checking if current_num is one of the x_values to collect
            while target_x_ptr <= length(x_values_to_collect) && x_values_to_collect(target_x_ptr) == current_num
                x_data_collected(collection_idx) = current_num;
                pi_x_data_collected(collection_idx) = current_total_prime_count;
                collection_idx = collection_idx + 1;
                target_x_ptr = target_x_ptr + 1;
            end
        end
    end

    % Printing a newline after the loop finishes
    fprintf('\nSieving complete.\n');

    % Trimming any unused pre-allocated space
    if collection_idx <= length(x_data_collected)
        x_data_collected = x_data_collected(1:collection_idx-1);
        pi_x_data_collected = pi_x_data_collected(1:collection_idx-1);
    end
end