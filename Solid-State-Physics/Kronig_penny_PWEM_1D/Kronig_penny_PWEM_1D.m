% =========================================================================
% Method: 1D Plane Wave Expansion Method (PWEM)
% Model: Kronig-Penney Periodic Rectangular Potential (Well-Barrier-Well)
% Outputs Generated:
%   1. Reduced Zone Scheme Band Structure with Shaded Forbidden Gaps
%   2. Real-Space Potential V(x) & Bloch Probability Densities at Gamma (k=0)
%   3. Parametric Evolution of Band Gaps (Gap 1 & Gap 2) vs Barrier Height V0
% =========================================================================

clear, clc
close all

% 1. Physical Parameters & Lattice Inputs
w1 = 0.3; % well-1 width
b  = 0.3; % Barrier width
w2 = 0.3; % well-2 width
d  = w1 + b + w2;  % Lattice constant (nm)
V0 = 4;      % Potential barrier height (eV)
hm = 0.0381; % hbar^2 / (2 * m_e) in (eV * nm^2)

% 2. Plane Wave Basis Setup
N = 10;                     % Maximum plane wave index
m = -N:N;                   
M = length(m);              
Mk = 10 * M;                % Number of k-points for plotting
k = linspace(-pi/d, pi/d, Mk); % Wave vector array across the 1st Brillouin Zone

% 3. Reciprocal Lattice Vectors & Potential Matrix (V)
gm = m * (2*pi/d);          % Reciprocal lattice vectors (G_m)
V = zeros(M);               

% The Fourier components of the periodic potential V(x)
for i_m = 1:M
    for i_n = 1:M
        p = m(i_m) - m(i_n); % Difference between plane wave indices
        if p == 0
            V(i_m, i_n) = V0 * b / d; % diagonal elements
        else
            % Off-diagonal Fourier components of the barrier
            term1 = exp(-1i * 2*pi * p * (w1 + b) / d);
            term2 = exp(-1i * 2*pi * p * w1 / d);
            V(i_m, i_n) = (V0 / (-1i * 2*pi * p)) * (term1 - term2);
        end
    end
end

% 4. Hamiltonian Construction & Eigenvalue 
E_all = zeros(M, Mk); 
psi_all = cell(Mk, 1);

for j = 1:length(k)
    % Diagonal kinetic energy matrix 
    KE = diag(hm * (k(j) + gm).^2);
    
    % Total Hamiltonian matrix
    H = KE + V;
    H = (H + H') / 2; % Ensure Hermiticity
    
    % Solve for eigenvalues and eigenvectors
    [vecs, vals] = eig(H);
    [E_sorted, sort_idx] = sort(real(diag(vals)));
    E_all(:, j) = E_sorted;
    psi_all{j} = vecs(:, sort_idx);
end

% 5. Physical Outputs Calculation Setup
num_bands_plot = 4;

% Identify Band Gaps between consecutive bands
gaps = [];
for n = 1:(num_bands_plot - 1)
    top_lower = max(E_all(n, :));
    bot_upper = min(E_all(n + 1, :));
    if bot_upper > top_lower
        gaps = [gaps; top_lower, bot_upper]; % [gap_start, gap_end]
    end
end

% Real-space Wavefunction reconstruction across 3 unit cells
x = linspace(-d, 2*d, 600);
Vx = zeros(size(x));
x_mod = mod(x, d);
Vx(x_mod >= w1 & x_mod <= (w1 + b)) = V0;

% Wavefunctions for Band 1 & Band 2 at Gamma-point (k = 0)
[~, gamma_idx] = min(abs(k));
psi1 = zeros(size(x));
psi2 = zeros(size(x));
c1 = psi_all{gamma_idx}(:, 1);
c2 = psi_all{gamma_idx}(:, 2);

for i_m = 1:M
    phase_factor = exp(1i * (k(gamma_idx) + gm(i_m)) * x);
    psi1 = psi1 + c1(i_m) * phase_factor;
    psi2 = psi2 + c2(i_m) * phase_factor;
end

prob1 = abs(psi1).^2;
prob2 = abs(psi2).^2;
scale_h = 1.2;
prob1_plot = (prob1 / max(prob1)) * scale_h + E_all(1, gamma_idx);
prob2_plot = (prob2 / max(prob2)) * scale_h + E_all(2, gamma_idx);

% Periodic Bloch functions u_n(x) at k = 0
u1 = zeros(size(x));
for i_m = 1:M
    u1 = u1 + c1(i_m) * exp(1i * gm(i_m) * x);
end

% Derivatives: Group Velocity and Effective Mass
dk = k(2) - k(1);
vg    = zeros(num_bands_plot, Mk);
m_eff = zeros(num_bands_plot, Mk);

for n = 1:num_bands_plot
    dE_dk   = gradient(E_all(n, :), dk);
    d2E_dk2 = gradient(dE_dk, dk);
    vg(n, :)    = dE_dk / (2 * hm);
    m_eff(n, :) = (2 * hm) ./ (d2E_dk2 + eps);
end

% 6. Parametric Studies (Gap vs V0, Bandwidth vs b, and Convergence vs N)

% Study 1: Band Gap vs Potential Height V0
V0_vals = linspace(0, 20, 40);
gap1_vs_V0 = zeros(size(V0_vals));
gap2_vs_V0 = zeros(size(V0_vals));

for idx_v = 1:length(V0_vals)
    V0_temp = V0_vals(idx_v);
    V_temp = zeros(M);
    for i_m = 1:M
        for i_n = 1:M
            p = m(i_m) - m(i_n);
            if p == 0
                V_temp(i_m, i_n) = V0_temp * b / d;
            else
                term1 = exp(-1i * 2*pi * p * (w1 + b) / d);
                term2 = exp(-1i * 2*pi * p * w1 / d);
                V_temp(i_m, i_n) = (V0_temp / (-1i * 2*pi * p)) * (term1 - term2);
            end
        end
    end
    
    % Check at Gamma (k=0) and Zone Boundary (k=pi/d) for band extrema
    H_gamma = diag(hm * (0 + gm).^2) + V_temp;
    H_edge  = diag(hm * (pi/d + gm).^2) + V_temp;
    H_gamma = (H_gamma + H_gamma') / 2;
    H_edge  = (H_edge + H_edge') / 2;
    
    E_g = sort(real(eig(H_gamma)));
    E_e = sort(real(eig(H_edge)));
    
    % Gap 1 occurs at k = pi/d (between band 1 and band 2)
    gap1_vs_V0(idx_v) = max(0, E_e(2) - E_e(1));
    % Gap 2 occurs at k = 0 (between band 2 and band 3)
    gap2_vs_V0(idx_v) = max(0, E_g(3) - E_g(2));
end

% Study 2: Bandwidth vs Barrier Width b (at fixed well width w1 = w2 = 0.3 nm)
b_vals = linspace(0.05, 0.6, 30);
bw_band1 = zeros(size(b_vals));
bw_band2 = zeros(size(b_vals));

for idx_b = 1:length(b_vals)
    b_curr = b_vals(idx_b);
    d_curr = w1 + b_curr + w2;
    gm_curr = m * (2*pi/d_curr);
    
    V_b = zeros(M);
    for i_m = 1:M
        for i_n = 1:M
            p = m(i_m) - m(i_n);
            if p == 0
                V_b(i_m, i_n) = V0 * b_curr / d_curr;
            else
                term1 = exp(-1i * 2*pi * p * (w1 + b_curr) / d_curr);
                term2 = exp(-1i * 2*pi * p * w1 / d_curr);
                V_b(i_m, i_n) = (V0 / (-1i * 2*pi * p)) * (term1 - term2);
            end
        end
    end
    
    H_g_b = diag(hm * (0 + gm_curr).^2) + V_b;
    H_e_b = diag(hm * (pi/d_curr + gm_curr).^2) + V_b;
    H_g_b = (H_g_b + H_g_b') / 2;
    H_e_b = (H_e_b + H_e_b') / 2;
    
    E_g_b = sort(real(eig(H_g_b)));
    E_e_b = sort(real(eig(H_e_b)));
    
    bw_band1(idx_b) = abs(E_e_b(1) - E_g_b(1));
    bw_band2(idx_b) = abs(E_e_b(2) - E_g_b(2));
end

% Study 3: PWEM Numerical Convergence vs Plane Wave Cutoff N
N_range = 1:30;
E1_conv = zeros(size(N_range));
E2_conv = zeros(size(N_range));

for idx_N = 1:length(N_range)
    N_temp = N_range(idx_N);
    m_temp = -N_temp:N_temp;
    M_temp = length(m_temp);
    gm_temp = m_temp * (2*pi/d);
    
    V_temp = zeros(M_temp);
    for i_m = 1:M_temp
        for i_n = 1:M_temp
            p = m_temp(i_m) - m_temp(i_n);
            if p == 0
                V_temp(i_m, i_n) = V0 * b / d;
            else
                term1 = exp(-1i * 2*pi * p * (w1 + b) / d);
                term2 = exp(-1i * 2*pi * p * w1 / d);
                V_temp(i_m, i_n) = (V0 / (-1i * 2*pi * p)) * (term1 - term2);
            end
        end
    end
    
    H_conv = diag(hm * (0 + gm_temp).^2) + V_temp;
    H_conv = (H_conv + H_conv') / 2;
    evals_conv = sort(real(eig(H_conv)));
    E1_conv(idx_N) = evals_conv(1);
    E2_conv(idx_N) = evals_conv(2);
end

% 7. Plotting Each Relationship in a Separate Figure

% Figure 1: Reduced Zone Scheme Band Structure with Shaded Gaps
figure('Color', 'w', 'Name', 'Reduced Zone Scheme');
hold on;
for g = 1:size(gaps, 1)
    fill([-pi/d, pi/d, pi/d, -pi/d], ...
         [gaps(g, 1), gaps(g, 1), gaps(g, 2), gaps(g, 2)], ...
         [0.9 0.9 0.9], 'EdgeColor', 'none', 'DisplayName', 'Forbidden Gap');
end
plot(k, E_all(1:num_bands_plot, :)', 'b', 'LineWidth', 1.8);
xline(-pi/d, '--r', '-\pi/d', 'LineWidth', 1.2);
xline(pi/d, '--r', '+\pi/d', 'LineWidth', 1.2);
xlim([-pi/d, pi/d]);
ylim([0, 25]);
xlabel('Wavevector k (nm^{-1})', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('1. Band Structure with Shaded Gaps (Reduced Zone)', 'FontSize', 13);
grid on;


% Figure 2: Real-Space Potential and Wavefunctions
figure('Color', 'w', 'Name', 'Real Space Wavefunctions');
plot(x, Vx, 'k--', 'LineWidth', 1.3, 'DisplayName', 'V(x)');
hold on;
plot(x, prob1_plot, 'b-', 'LineWidth', 1.8, 'DisplayName', '|\psi_1(x)|^2');
plot(x, prob2_plot, 'r-', 'LineWidth', 1.8, 'DisplayName', '|\psi_2(x)|^2');
yline(E_all(1, gamma_idx), 'b:');
yline(E_all(2, gamma_idx), 'r:');
xlim([-d, 2*d]);
ylim([0, V0 + 5]);
xlabel('Position x (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy / Potential (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('2. Real-Space Potential & Probability Densities at \Gamma (k = 0)', 'FontSize', 13);
legend('Location', 'northeast');
grid on;

% Figure 3: Band Gaps vs Barrier Height V0
figure('Color', 'w', 'Name', 'Band Gaps vs V0');
plot(V0_vals, gap1_vs_V0, 'b-o', 'LineWidth', 1.6, 'MarkerSize', 4, 'DisplayName', '1st Forbidden Gap');
hold on;
plot(V0_vals, gap2_vs_V0, 'r-s', 'LineWidth', 1.6, 'MarkerSize', 4, 'DisplayName', '2nd Forbidden Gap');
xlabel('Barrier Height V_0 (eV)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Band Gap Width (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('3. Evolution of Band Gaps vs. Barrier Height V_0', 'FontSize', 13);
legend('Location', 'northwest');
grid on;

% 8. Save Figures 
[script_folder, script_name] = fileparts(mfilename('fullpath'));
if isempty(script_name)
    script_name = 'Kronig_penny_PWEM';
    script_folder = pwd;
end
output_dir = fullfile(script_folder, script_name);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fig_handles = findobj('Type', 'figure');
for i_fig = 1:length(fig_handles)
    fig = fig_handles(i_fig);
    if isvalid(fig)
        fig_name = get(fig, 'Name');
        if isempty(fig_name)
            fig_name = sprintf('figure_%d', get(fig, 'Number'));
        else
            fig_name = regexprep(fig_name, '[^\w\-_]', '_');
        end
        file_path = fullfile(output_dir, sprintf('%s.png', fig_name));
        exportgraphics(fig, file_path, 'Resolution', 300);
    end
end


