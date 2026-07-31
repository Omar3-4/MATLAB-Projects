function [path_samples, eps_opt, acc_rate] = HMC(V_func, grad_V_func, D, N_samples, L_steps, m, dtau)
% HMC Performs Hamiltonian Monte Carlo using a spectral method.
%
% Inputs:
%   V_func      - Function handle for the potential V(x)
%   grad_V_func - Function handle for the gradient of V(x)
%   D           - Number of dimensions (default: 5000)
%   N_samples   - Number of MCMC samples to generate (default: 2000)
%   L_steps     - Number of integration steps per trajectory (default: 15)
%   m           - Mass parameter (default: 1.0)
%   dtau        - Time step for the path integral (default: 0.1)
%
% Outputs:
%   path_samples - (N_samples x D) matrix of generated samples
%   eps_opt      - Optimized step size from dual averaging
%   acc_rate     - Overall acceptance rate

    % PROTECTION: Validate inputs to prevent downstream dimension/type errors
    if nargin < 2
        error('At least V_func and grad_V_func must be provided.');
    end
    if ~isa(V_func, 'function_handle') || ~isa(grad_V_func, 'function_handle')
        error('V_func and grad_V_func must be valid function handles.');
    end
    if nargin < 3 || isempty(D) || ~isnumeric(D) || D <= 0, D = 5000; end
    if nargin < 4 || isempty(N_samples) || ~isnumeric(N_samples) || N_samples <= 0, N_samples = 2000; end
    if nargin < 5 || isempty(L_steps) || ~isnumeric(L_steps) || L_steps <= 0, L_steps = 15; end                 
    if nargin < 6 || isempty(m) || ~isnumeric(m) || m <= 0, m = 1.0; end
    if nargin < 7 || isempty(dtau) || ~isnumeric(dtau) || dtau <= 0, dtau = 0.1; end

    try
        % PROTECTION: Catch invalid matrix/math operations during parameter initialization
        k_idx   = (1:D)';
        lambda_k = 4 * (sin(pi * k_idx / (2 * (D + 1)))).^2;
        
        K_k     = (m / dtau) * lambda_k;
        M_k     = K_k + (m * dtau); 
        M_inv_k = 1.0 ./ M_k;
        sqrt_M  = sqrt(M_k);
        eps       = 0.25;               
        target_acc= 0.75;
        gamma_da  = 0.05;
        t0_da     = 10;
        kappa_da  = 0.75;
        mu_da     = log(10 * eps);
        log_eps_hat = log(eps);
        H_bar     = 0;
        burn_max  = 300;
        a_c = zeros(D, 1);
    catch ME
        % FALLBACK: Abort execution immediately if initialization math fails
        error('Initialization failed due to invalid dimensions or parameters.');
    end
    
    try
        % PROTECTION: Catch errors in initial action/gradient computation
        U_c    = compute_action_spectral(a_c, K_k, V_func, dtau);
        grad_c = compute_grad_spectral(a_c, K_k, grad_V_func, dtau);
    catch
        % FALLBACK: Return empty outputs to prevent hard crash
        path_samples = NaN(N_samples, D); eps_opt = NaN; acc_rate = 0; return;
    end

    for burn = 1:burn_max
        try
            % PROTECTION: Catch runtime crashes during burn-in Omelyan integration
            p_c = randn(D, 1) .* sqrt_M;
            H_c = U_c + 0.5 * sum((p_c.^2) .* M_inv_k);
            [a_prop, p_prop, U_prop, grad_prop, valid] = ...
                omelyan_spectral_step(a_c, p_c, grad_c, eps, L_steps, K_k, M_inv_k, V_func, grad_V_func, dtau);
        catch
            % FALLBACK: Flag step as invalid to force rejection and continue burn-in
            valid = false; U_prop = Inf; 
        end

        if valid && isfinite(U_prop)
            H_prop = U_prop + 0.5 * sum((p_prop.^2) .* M_inv_k);
            dH     = H_prop - H_c;
            alpha  = min(1.0, exp(-dH));
            if isnan(alpha), alpha = 0.0; end
        else
            alpha  = 0.0;
        end

        if rand() < alpha
            a_c    = a_prop;
            U_c    = U_prop;
            grad_c = grad_prop;
        end

        eta_da      = 1 / (burn + t0_da);
        H_bar       = (1 - eta_da) * H_bar + eta_da * (target_acc - alpha);
        log_eps     = mu_da - (sqrt(burn) / gamma_da) * H_bar;
        m_eta       = burn^(-kappa_da);
        log_eps_hat = m_eta * log_eps + (1 - m_eta) * log_eps_hat;
        eps         = exp(log_eps);
    end

    eps_opt = exp(log_eps_hat);
    spectral_samples = zeros(N_samples, D);
    acc_count = 0;

    for n = 1:N_samples
        try
            % PROTECTION: Catch runtime crashes during sampling Omelyan integration
            p_c = 0.85 * p_c + sqrt(1 - 0.85^2) * (randn(D, 1) .* sqrt_M);
            H_c = U_c + 0.5 * sum((p_c.^2) .* M_inv_k);
            [a_prop, p_prop, U_prop, grad_prop, valid] = ...
                omelyan_spectral_step(a_c, p_c, grad_c, eps_opt, L_steps, K_k, M_inv_k, V_func, grad_V_func, dtau);
        catch
            % FALLBACK: Flag step as invalid to force rejection
            valid = false; U_prop = Inf;
        end

        if valid && isfinite(U_prop)
            H_prop = U_prop + 0.5 * sum((p_prop.^2) .* M_inv_k);
            dH     = H_prop - H_c;
            if log(rand()) < -dH
                a_c       = a_prop;
                p_c       = p_prop;
                U_c       = U_prop;
                grad_c    = grad_prop;
                acc_count = acc_count + 1;
            else
                p_c       = -p_c; 
            end
        else
            p_c       = -p_c; 
        end
        spectral_samples(n, :) = a_c';
    end

    acc_rate = acc_count / N_samples;
    path_samples = zeros(N_samples, D);

    for n = 1:N_samples
        try
            % PROTECTION: Catch DST transform failures for individual samples
            path_samples(n, :) = dst_type1(spectral_samples(n, :)');
        catch
            % FALLBACK: Output NaNs for the failed sample row instead of crashing
            path_samples(n, :) = NaN(1, D);
        end
    end
    
    fprintf('\n===================================================\n');
    fprintf('   SPECTRAL HMC PATH INTEGRAL RESULTS (D = %d)\n', D);
    fprintf('===================================================\n');
    fprintf(' Optimal Step Size (eps) : %.5f\n', eps_opt);
    fprintf(' Final Acceptance Rate   : %.2f%%\n', acc_rate * 100);
end

function U = compute_action_spectral(a, K_k, V_func, dtau)
    try
        % PROTECTION: Catch evaluation errors in user-provided V_func
        x = dst_type1(a);
        U = 0.5 * sum((a.^2) .* K_k) + dtau * sum(V_func(x));
    catch
        % FALLBACK: Return Inf to force HMC step rejection
        U = Inf;
    end
end

function grad = compute_grad_spectral(a, K_k, grad_V_func, dtau)
    try
        % PROTECTION: Catch evaluation errors in user-provided grad_V_func
        x = dst_type1(a);
        grad_x = dtau * grad_V_func(x);
        grad = K_k .* a + dst_type1(grad_x);
    catch
        % FALLBACK: Return NaN array to invalidate gradient and step
        grad = NaN(size(a));
    end
end

function [a, p, U_val, grad_val, valid] = omelyan_spectral_step(a, p, grad, eps, L, K_k, M_inv_k, V_func, grad_V_func, dtau)
    valid = true;
    xi    = 0.1931833275037836; 
    U_val = Inf;
    grad_val = zeros(size(a));
    
    try
        % PROTECTION: Catch dimension mismatches or math crashes in leapfrog integration
        for l = 1:L
            p = p - xi * eps * grad;
            a = a + 0.5 * eps * (p .* M_inv_k);
            
            grad = compute_grad_spectral(a, K_k, grad_V_func, dtau);
            if any(~isfinite(a)) || any(~isfinite(grad)), valid = false; return; end
            
            p = p - (1 - 2 * xi) * eps * grad;
            a = a + 0.5 * eps * (p .* M_inv_k);
            
            grad = compute_grad_spectral(a, K_k, grad_V_func, dtau);
            if any(~isfinite(a)) || any(~isfinite(grad)), valid = false; return; end
            
            p = p - xi * eps * grad;
        end
        U_val    = compute_action_spectral(a, K_k, V_func, dtau);
        grad_val = grad;
        if ~isfinite(U_val), valid = false; end
    catch
        % FALLBACK: Flag step as invalid on critical math failure
        valid = false;
    end
end

function x = dst_type1(a)
    D = length(a);
    try
        % PROTECTION: Catch FFT errors or matrix dimension mismatches
        y = zeros(2*(D+1), 1);
        y(2:D+1) = a;
        y(D+3:end) = -flipud(a);
        
        y_fft = fft(y);
        x = -0.5 * imag(y_fft(2:D+1)) / sqrt((D + 1) / 2);
    catch
        % FALLBACK: Return NaN array to propagate invalidation safely
        x = NaN(D, 1);
    end
end