function [f_samples, F_Sigma, acc_rate, T_burnSteps] = MCMC_MH(log_f, x0, N, Sigma, K)
% MCMC_MH - Adaptive Metropolis-Hastings MCMC Sampler
%
% Features:
%   - Robbins-Monro stochastic adaptation for optimal proposal step size (Sigma)
%   - Geweke convergence diagnostic corrected for autocorrelation using FFT
%
% Inputs:
%   log_f       - Function handle to the LOG target probability density: log_f(x)
%   x0          - Initial state vector (1 x D)
%   N           - Number of thinned samples to extract
%   Sigma       - Initial proposal standard deviation (scalar or 1xD vector)
%   K           - Thinning factor (retains 1 out of every K samples)
%
% Outputs:
%   f_samples   - Generated thinned samples matrix (N x D)
%   F_Sigma     - Final adapted proposal standard deviation (Sigma)
%   acc_rate    - Acceptance rate during the production sampling phase
%   T_burnSteps - Total iterations spent in the burn-in phase

    % --- 0. Input Parameter Validation ---
    validateattributes(x0, {'numeric'}, {'2d', 'row'}, mfilename, 'x0', 2);
    validateattributes(N,  {'numeric'}, {'scalar', 'positive', 'integer'}, mfilename, 'N', 3);
    validateattributes(K,  {'numeric'}, {'scalar', 'positive', 'integer'}, mfilename, 'K', 5);

    D = size(x0, 2);               % Space dimension
    B_size = 1000;                 % Batch size for diagnostic and tuning
    Temp_batch = zeros(B_size, D); % Preallocate memory for batch samples
    T_burnSteps = 0;               % Burn-in steps counter
    x_c = x0;                      % Current state
    log_f_c = log_f(x_c);          % Evaluate current log-density
    Burn = false;                  % Convergence flag
    
    % --- Robbins-Monro Initialization ---
    batch_count = 0;               % Counter for adaptive batches
    if D == 1
        target_rate = 0.440;       % Optimal target acceptance rate for D = 1
    else
        target_rate = 0.234;       % Optimal target acceptance rate for D > 1
    end

    % =========================================================================
    % PHASE 1: Adaptive Burn-In Phase (Robbins-Monro + FFT Geweke)
    % =========================================================================
    while ~Burn
        acc_count = 0;

        % --- Run 1000 iterations for current batch ---
        for i = 1:B_size
            % Propose new candidate using Gaussian random walk
            x_n = x_c + (randn(1, D) .* Sigma);
            log_f_n = log_f(x_n);

            % Metropolis log-acceptance criterion
            log_a = log_f_n - log_f_c;

            if log(rand()) < log_a
                x_c = x_n;
                log_f_c = log_f_n;
                acc_count = acc_count + 1;
            end

            Temp_batch(i, :) = x_c;
        end

        T_burnSteps = T_burnSteps + B_size;
        batch_count = batch_count + 1;
        
        % Calculate actual batch acceptance rate
        actual_rate = acc_count / B_size;

        % ---------------------------------------------------------------------
        % Robbins-Monro Self-Tuning for Sigma
        % ---------------------------------------------------------------------
        % Learning rate decays with batch count to guarantee ergodic convergence
        gamma = 1 / sqrt(batch_count);
        
        % Log-space update to keep Sigma strictly positive:
        % log(Sigma_new) = log(Sigma_old) + gamma * (Actual_Rate - Target_Rate)
        Sigma = exp(log(Sigma) + gamma * (actual_rate - target_rate));

        % ---------------------------------------------------------------------
        % Geweke Diagnostic with FFT Autocorrelation Correction
        % ---------------------------------------------------------------------
        % Split batch into early 10% (P_A) and late 50% (P_B) segments
        N_A = round(0.1 * B_size);
        N_B = round(0.5 * B_size);

        P_A = Temp_batch(1:N_A, :);
        P_B = Temp_batch((B_size - N_B + 1):B_size, :);

        % Calculate segment means
        mean_A = mean(P_A, 1);
        mean_B = mean(P_B, 1);

        % Compute True Variances (Standard Error of the Mean) corrected for 
        % autocorrelation time (tau) using FFT routines
        sem_A = compute_fft_true_variance(P_A);
        sem_B = compute_fft_true_variance(P_B);

        % Standard Z-score ratio with corrected variances
        Z_scores = abs(mean_A - mean_B) ./ sqrt(sem_A + sem_B + eps);

        % Convergence check (|Z| < 1.96 corresponds to ~95% confidence level)
        if max(Z_scores) < 1.96
            Burn = true;
        end
    end

    % =========================================================================
    % PHASE 2: Main Production Sampling Phase (Frozen Proposal)
    % =========================================================================
    T_P_burnSteps = N * K;
    f_samples = zeros(N, D);
    p_acc_count = 0;
    sample_idx = 1;

    for i = 1:T_P_burnSteps
        x_n = x_c + (randn(1, D) .* Sigma);
        log_f_n = log_f(x_n);

        log_a = log_f_n - log_f_c;

        if log(rand()) < log_a
            x_c = x_n;
            log_f_c = log_f_n;
            p_acc_count = p_acc_count + 1;
        end

        % Save thinned sample every K steps
        if mod(i, K) == 0
            f_samples(sample_idx, :) = x_c;
            sample_idx = sample_idx + 1;
        end
    end

    % Final metrics
    F_Sigma = Sigma;
    acc_rate = p_acc_count / T_P_burnSteps;
end

% =========================================================================
% LOCAL HELPER: FFT-Based True Variance Computation
% =========================================================================
function sem = compute_fft_true_variance(Data)
% Calculates the True Variance (SEM) accounting for autocorrelation time (tau)
% Steps:
%   1. Subtract mean (center data)
%   2. Fast Fourier Transform (FFT) -> frequency domain
%   3. Compute Power Spectral Density (PSD = abs(fft)^2)
%   4. Inverse FFT (IFFT) -> returns Autocorrelation sequence
%   5. Integrate autocorrelation to find correlation time (tau)
%   6. True Variance = (var(Data) / N) * tau

    [N, D] = size(Data);
    sem = zeros(1, D);

    for d = 1:D
        x = Data(:, d);
        
        % 1. Mean Subtraction (Center sequence around 0)
        x_centered = x - mean(x);

        % 2 & 3. FFT & Power Spectral Density
        % Zero-padding to twice length to prevent circular autocorrelation aliasing
        nfft = 2^nextpow2(2 * N - 1);
        X_fft = fft(x_centered, nfft);
        psd = abs(X_fft).^2;

        % 4. Inverse FFT to compute Autocorrelation function
        autocorr_raw = ifft(psd);
        autocorr_raw = real(autocorr_raw(1:N)); % Retain positive lags

        % Normalize autocorrelation (rho_0 = 1)
        if autocorr_raw(1) == 0
            rho = zeros(N, 1);
        else
            rho = autocorr_raw / autocorr_raw(1);
        end

        % 5. Compute Correlation Time (tau)
        % Sum autocorrelations up to the first negative crossing to reduce tail noise
        first_neg = find(rho < 0, 1, 'first');
        if isempty(first_neg)
            max_lag = N;
        else
            max_lag = first_neg - 1;
        end

        % Integrated Autocorrelation Time tau = 1 + 2 * sum(rho(2:max_lag))
        tau = 1 + 2 * sum(rho(2:max_lag));
        tau = max(1, tau); % Ensure tau is strictly >= 1

        % 6. True Variance calculation: (var / N) * tau
        sem(d) = (var(x) / N) * tau;
    end
end
