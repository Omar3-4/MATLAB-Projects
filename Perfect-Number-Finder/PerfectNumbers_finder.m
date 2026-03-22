% =========================================================
% Perfect Number Finder using GPU Acceleration
% Searches for perfect numbers in ranges from 10^1 to 10^12
% A perfect number equals the sum of all its proper divisors
% =========================================================

tic;
total_start_time = tic;   % Timer for total elapsed time since program start
last_perfect_time = tic;  % Timer to measure gap between found perfect numbers
Perfect_N = [];            % Array to store all perfect numbers found
max_chunk = uint64(1e8);   % Max batch size to avoid overloading GPU memory

% --- Loop over 12 powers of 10 (ranges: 1-10, 10-100, ..., 10^11 - 10^12) ---
for j = 1:12
    total_start = uint64(10^(j-1));  % Start of current range
    total_end   = uint64(10^j);      % End of current range

    % Make sure we start from an even number (we only check even numbers)
    if mod(total_start, 2) ~= 0, total_start = total_start + 1
    end

    % Calculate chunk size to split the range into smaller GPU-friendly batches
    current_chunk = min(max_chunk, (total_end - total_start)/2);

    % --- Loop over chunks within the current range ---
    for current_start = total_start : uint64(2 * current_chunk) : total_end-1
        current_end = min(current_start + uint64(2 * current_chunk) - 2, total_end);

        % Load even numbers in this chunk onto the GPU
        n = gpuArray(double(current_start : 2 : current_end));

        % --- Filter 1: Keep only numbers ending in 6 or 8 ---
        % (All even perfect numbers end in 6 or 8)
        n = n(mod(n, 10) == 6 | mod(n, 10) == 8);
        if isempty(n)
            wait(gpuDevice); clear n; continue
        end

        % --- Filter 2: Keep numbers where digit sum is divisible by 9, or n=6 ---
        % (Based on digit-sum property of known perfect numbers)
        n = n(mod(n-1, 9) + 1 == 1 | n == 6);
        if isempty(n)
            wait(gpuDevice); clear n; continue
        end

        % --- Filter 3: Check if (8n + 1) is a perfect odd square ---
        % This comes from the formula: if n = k*(k+1)/2, then 8n+1 = (2k+1)^2
        roots = sqrt(8 * double(n) + 1);
        mask = (mod(roots, 1) == 0) & (mod(roots, 2) ~= 0); % Must be integer and odd
        n = n(mask);
        roots = roots(mask);

        % --- Process remaining candidates on the CPU ---
        if ~isempty(n)
            candidates = gather(uint64(n));          % Move data back from GPU to CPU
            candidate_roots = gather(roots);

            for i = 1:length(candidates)
                % Recover k from the triangular number formula: n = k*(k+1)/2
                k = uint64(0.5 * (candidate_roots(i) - 1));

                % Euler's theorem: even perfect numbers require k to be prime
                if isprime(k)

                    % --- Compute the sum of proper divisors of candidates(i) ---
                    f = floor(sqrt(double(candidates(i))));  % Only check up to sqrt(n)
                    F = uint64(f);
                    S = uint64(1);  % Start with 1 (always a divisor)

                    for T = 2:F
                        if mod(candidates(i), T) == 0       % T is a divisor
                            p = candidates(i) / T;           % Paired divisor
                            if T == p
                                S = S + p;                   % Perfect square case: add once
                            else
                                S = S + p + T;               % Add both divisors
                            end
                            if S > candidates(i), break; end % Early exit if sum exceeds n
                        end
                    end

                    % --- Check if sum of divisors equals the number itself ---
                    if S == candidates(i)
                        Perfect_N(end+1) = candidates(i);
                        fprintf('Perfect Number Found: %d\n', candidates(i));

                        % Calculate timing info
                        time_from_start = toc(total_start_time);  % Time since program started
                        time_since_last = toc(last_perfect_time);  % Time since last perfect number
                        last_perfect_time = tic;                   % Reset gap timer

                        % Print match details to console
                        fprintf('\n[MATCH FOUND: %lu]', candidates(i));
                        fprintf('\n- Total: %.2f s', time_from_start);
                        fprintf('\n- Gap:   %.2f s\n', time_since_last);

                        % Append result to log file immediately
                        fileID_instant = fopen('Instant_Perfect_Log.txt', 'a');
                        fprintf(fileID_instant, '%lu | Start: %.2f | Gap: %.2f | Time: %s\n', ...
                                candidates(i), time_from_start, time_since_last, datestr(now));
                        fclose(fileID_instant);

                        % Save backup of all found perfect numbers to .mat file
                        save('Perfect_Backup.mat', 'Perfect_N');
                    end
                end
            end
        end

        % Free memory after each chunk
        clear n roots mask candidates candidate_roots;
    end
end

% --- Final Summary ---
final_total_time = toc(total_start_time);
fprintf('\n-----------------------------------------');
fprintf('\nSearch Complete: %.2f mins\n', final_total_time/60);
disp('Final List:'); disp(Perfect_N);
