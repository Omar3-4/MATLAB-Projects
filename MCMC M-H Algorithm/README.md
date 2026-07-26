# Adaptive Metropolis-Hastings MCMC with FFT-Corrected Geweke Diagnostic

`MCMC_MH.m`, a highly optimized, self-tuning Markov Chain Monte Carlo (MCMC) sampler implemented in MATLAB. It utilizes the Metropolis-Hastings algorithm enhanced with stochastic optimization and spectral variance analysis to automate the most notoriously difficult aspects of MCMC sampling: proposal tuning and burn-in convergence detection.

## Core Algorithmic Features

### 1. Robbins-Monro Stochastic Adaptation

Choosing the correct proposal standard deviation ($\Sigma$) is critical for MCMC efficiency. If $\Sigma$ is too small, the chain explores too slowly (high acceptance, high autocorrelation). If $\Sigma$ is too large, proposals are constantly rejected.

This algorithm automates the tuning of $\Sigma$ during the burn-in phase using a **Robbins-Monro stochastic approximation**. The step size is updated in logarithmic space to guarantee positivity, driving the empirical acceptance rate ($\alpha$) toward the theoretical optimum ($\alpha^* = 0.44$ for 1D, $\alpha^* \approx 0.234$ for multi-dimensional spaces). The adaptation sequence follows a decaying learning rate ($\gamma_n = 1/\sqrt{n}$) to ensure ergodic convergence.

### 2. Automated FFT-Corrected Geweke Diagnostic

Unlike standard samplers that require the user to guess the burn-in length, this algorithm automatically detects stationarity. It utilizes the **Geweke diagnostic**, which compares the means of the early (10%) and late (50%) segments of the sample batches.

Crucially, standard variance calculations fail for MCMC due to sample correlation. This script employs a **Fast Fourier Transform (FFT)** to map the data to the frequency domain, calculate the Power Spectral Density (PSD), and derive the autocorrelation function via Inverse FFT. This allows for the exact calculation of the Integrated Autocorrelation Time ($\tau$), yielding the _True Variance_ (Standard Error of the Mean) to prevent premature convergence declarations.

## Function Signature

```
[f_samples, F_Sigma, acc_rate, T_burnSteps] = MCMC_MH(log_f, x0, N, Sigma, K)
```

### Inputs

|Parameter|Type|Description|
|---|---|---|
|`log_f`|Function Handle|The target probability density function in logarithmic space (to prevent numerical underflow). Format: `log_f(x)`.|
|`x0`|Vector (1 x D)|The initial starting state/coordinates for the Markov chain.|
|`N`|Integer|The desired number of final, _thinned_ samples to extract.|
|`Sigma`|Scalar / Vector|The initial guess for the proposal standard deviation.|
|`K`|Integer|The thinning factor (retains 1 out of every `K` samples) to minimize residual autocorrelation in the final dataset.|

### Outputs

|Parameter|Type|Description|
|---|---|---|
|`f_samples`|Matrix (N x D)|The final generated thinned samples from the target distribution.|
|`F_Sigma`|Scalar / Vector|The optimal proposal standard deviation found by the Robbins-Monro algorithm.|
|`acc_rate`|Float|The acceptance rate maintained during the Phase 2 production sampling.|
|`T_burnSteps`|Integer|The total number of iterations automatically spent in Phase 1 before convergence was achieved.|

## Internal Architecture

The script strictly divides execution into two phases to maintain detailed balance and Markovian properties in the final samples.

### Phase 1: Adaptive Burn-In

The algorithm processes samples in batches of 1,000. For each batch:

1. Samples are drawn using the current $\Sigma$.
    
2. The acceptance rate is measured, and $\Sigma$ is updated via the Robbins-Monro rule.
    
3. The FFT-based Geweke diagnostic is calculated.
    
4. If the absolute Z-score of the segment differences falls below 1.96 (95% confidence interval), the chain is declared stationary. Otherwise, a new batch begins.
    

### Phase 2: Production Sampling

Once stationarity is achieved, the adaptation immediately ceases. $\Sigma$ is frozen at its optimal value. The algorithm then generates `N * K` samples, saving every `K`-th sample to the output matrix to provide a statistically independent dataset representing the target distribution.

## Usage Example: 2D Correlated Gaussian

This example demonstrates how to sample from a highly correlated 2D Gaussian distribution, a common test for MCMC efficiency.

```
% 1. Define the target distribution parameters
mu = [0, 0];
Covariance = [1.0, 0.8; 
              0.8, 1.0];
invCov = inv(Covariance);

% 2. Define the LOG target density function
% Note: Normalization constants can be omitted in MCMC
log_pdf = @(x) -0.5 * (x - mu) * invCov * (x - mu)';

% 3. Set MCMC parameters
x0 = [10, -10];     % Intentional bad starting point to test burn-in
N = 5000;           % Request 5000 samples
Sigma_init = 1.0;   % Initial proposal step size
K = 10;             % Thinning factor

% 4. Run the Adaptive Sampler
tic;
[samples, optimal_sigma, acc_rate, burn_steps] = MCMC_MH(log_pdf, x0, N, Sigma_init, K);
execution_time = toc;

% 5. Display Results
fprintf('MCMC Completed in %.2f seconds.\n', execution_time);
fprintf('Total Burn-in Steps required: %d\n', burn_steps);
fprintf('Optimal Proposal Sigma: %.4f\n', optimal_sigma);
fprintf('Production Acceptance Rate: %.2f%%\n', acc_rate * 100);

% 6. Visualization
scatter(samples(:,1), samples(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.5);
title('MCMC Samples: 2D Correlated Gaussian');
xlabel('X_1'); ylabel('X_2');
grid on;
```

## Mathematical Background

**Robbins-Monro Log-Update Rule:**

$$\log(\Sigma_{n+1}) = \log(\Sigma_n) + \frac{1}{\sqrt{n}} (\alpha_{\text{actual}} - \alpha_{\text{target}})$$

**Geweke Z-Score:**

$$Z = \frac{\bar{\theta}_A - \bar{\theta}_B}{\sqrt{\text{var}(\bar{\theta}_A) + \text{var}(\bar{\theta}_B)}}$$

**FFT-Derived True Variance:** Given a sequence $X$, the power spectral density $S(f)$ is computed via FFT. The autocorrelation function $\rho(k)$ is extracted via IFFT. The integrated autocorrelation time $\tau$ and the resulting true standard error of the mean (SEM) are:

$$\tau = 1 + 2 \sum_{k=1}^{\infty} \rho(k)$$$$\text{True Variance} = \frac{\sigma^2}{N} \times \tau$$

## Requirements

- MATLAB (R2016a or later recommended).
    
- No external toolboxes are required. Core MATLAB functions (`fft`, `ifft`, `mean`, `var`) handle all advanced spectral computations.
