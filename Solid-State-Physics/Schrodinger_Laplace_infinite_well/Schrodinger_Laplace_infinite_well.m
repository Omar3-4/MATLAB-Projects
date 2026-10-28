% =========================================================================
% Method: Analytical Solution of 1D TISE via Laplace Transformation
% Model: 1D Infinite Potential Well (Particle in a Box)
% Outputs:
%   1. Wavefunctions psi_n(x) Shifted on Energy Levels with Infinite Barriers
%   2. Probability Densities |psi_n(x)|^2 on Quantized Energy Levels
%   3. Quantized Energy Spectrum (E_n proportional to n^2)
%   4. Quantum Confinement Effect (E_n vs Well Width L)
%   5. Momentum Space Wavefunctions |psi(p)|^2 (Fourier Domain)
% =========================================================================

clear, clc, close all

%% 1. Setting Differential Equation & Laplace Transform Solution
syms k s x c1 m E hbar Phi
syms L positive

% Using Laplace Transformation to solve the differential equation: psi''(x) + k^2*psi(x) = 0
laplaceEqn = (s^2*Phi - s*0 - c1) + k^2*Phi;

% Separating the sides of the equation
Phi_sol = solve(laplaceEqn, Phi);

% Use inverse Laplace transform to get real-space function
Laplaceinv_sol = ilaplace(Phi_sol, s, x);

% Apply boundary condition: psi(L) = 0
boundary_eq = subs(Laplaceinv_sol, x, L);

% Quantization from sin(kL) = 0
syms n positive
k_quant = n*pi/L;

% Update solution by substituting quantized k
Laplaceinv_sol = subs(Laplaceinv_sol, k, k_quant);

% Using normalization to get coefficient c1: integral(psi^2, 0, L) = 1
normalization = int(Laplaceinv_sol^2, x, 0, L) - 1;

% Solve for c1
C1 = solve(normalization, c1);

% Substitute sin(2*n*pi) = 0 and take positive amplitude
C1_val = subs(abs(C1(1)), sin(2*n*pi), 0);

% Substitute c1 into solution
Psi_Solution = subs(Laplaceinv_sol, c1, C1_val);

% Simplification of final analytic wavefunction
Psi_Final = simplify(Psi_Solution)
E_analytic = solve(k_quant^2 - 2*m*E/hbar^2, E)

%% 2. Numerical Substitution & Evaluation Directly from Analytic Results
L_val = 1e-9;
hp = 6.626e-34;
hbar_val = hp / (2*pi);
me = 9.109e-31;
q = 1.602e-19;
c_light = 3e8;

x_v = linspace(0, L_val, 200);
x_nm = x_v * 1e9;
L_nm = L_val * 1e9;

% Evaluate ground-state energy E1 directly from analytic expression E_analytic (at n=1)
E1 = double(subs(E_analytic, {hbar, m, L, n}, {hp/(2*pi), me, L_val, 1})) / q;

% Scaling factor exactly as in original code
sc = E1 / sqrt(2 / L_val);

num_states = 5;
En_vals = zeros(1, num_states);
psi_matrix = zeros(num_states, length(x_v));
prob_matrix = zeros(num_states, length(x_v));

% Evaluate wavefunctions and energy levels directly from analytic symbolic formulas: Psi_Final and E_analytic
for n_val = 1:num_states
    En_vals(n_val) = double(subs(E_analytic, {hbar, m, L, n}, {hp/(2*pi), me, L_val, n_val})) / q;
    psi_sym_n = subs(Psi_Final, {L, n}, {L_val, n_val});
    psi_matrix(n_val, :) = double(subs(psi_sym_n, x, x_v));
    prob_matrix(n_val, :) = psi_matrix(n_val, :).^2;
end

%% 3. Independent Figures for Each Relationship

% Figure 1: Wavefunctions psi_n(x) Shifted on Energy Levels with Infinite Walls
figure('Color', 'w', 'Name', 'Wavefunctions_on_Energy_Levels');
hold on;
for n_idx = 1:num_states
    En = En_vals(n_idx);
    % Baseline energy level
    plot(x_nm, En * ones(size(x_nm)), '--k', 'LineWidth', 1);
    % Shifted wavefunction using sc directly from original code
    plot(x_nm, En + sc * psi_matrix(n_idx, :), 'LineWidth', 2);
    text(L_nm * 1.02, En, sprintf('n = %d  (%.2f eV)', n_idx, En), 'FontSize', 10, 'FontWeight', 'bold');
end
% Draw infinite potential barriers
xline(0, '-k', 'LineWidth', 2.5);
xline(L_nm, '-k', 'LineWidth', 2.5);
xlim([-0.15 * L_nm, 1.25 * L_nm]);
ylim([0, En_vals(end) + 2*E1]);
xlabel('Position x (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy (eV) / Wavefunction \psi_n(x)', 'FontSize', 12, 'FontWeight', 'bold');
title('1. Wavefunctions \psi_n(x) on Quantized Energy Levels', 'FontSize', 13);
grid on;

% Figure 2: Probability Densities |psi_n(x)|^2
figure('Color', 'w', 'Name', 'Probability_Densities');
hold on;
sc_prob = E1 / (max(prob_matrix(:)) * 1.3);
for n_idx = 1:num_states
    En = En_vals(n_idx);
    plot(x_nm, En * ones(size(x_nm)), '--k', 'LineWidth', 0.8);
    plot(x_nm, En + sc_prob * prob_matrix(n_idx, :), 'LineWidth', 2);
    text(L_nm * 1.02, En, sprintf('n = %d', n_idx), 'FontSize', 10, 'FontWeight', 'bold');
end
xline(0, '-k', 'LineWidth', 2.5);
xline(L_nm, '-k', 'LineWidth', 2.5);
xlim([-0.15 * L_nm, 1.25 * L_nm]);
ylim([0, En_vals(end) + 2*E1]);
xlabel('Position x (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy (eV) / Probability Density |\psi_n(x)|^2', 'FontSize', 12, 'FontWeight', 'bold');
title('2. Probability Densities |\psi_n(x)|^2 on Energy Levels', 'FontSize', 13);
grid on;

% Figure 3: Quantized Energy Spectrum (E_n proportional to n^2)
figure('Color', 'w', 'Name', 'Energy_Levels_Spectrum');
hold on;
for n_idx = 1:num_states
    En = En_vals(n_idx);
    plot([0.5, 2.5], [En, En], 'b-', 'LineWidth', 2.5);
    text(2.6, En, sprintf('E_%d = %d^2 E_1 = %.3f eV', n_idx, n_idx, En), 'FontSize', 11, 'FontWeight', 'bold');
end
xlim([0, 4.2]);
ylim([0, En_vals(end) + 2*E1]);
set(gca, 'XTick', []);
ylabel('Energy (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('3. Quantized Energy Spectrum (E_n \propto n^2)', 'FontSize', 13);
grid on;

% Figure 4: Quantum Confinement: Energy vs. Well Width L
figure('Color', 'w', 'Name', 'Energy_vs_Well_Width');
L_sweep = linspace(0.5e-9, 3.0e-9, 100);
L_sweep_nm = L_sweep * 1e9;
hold on;
colors = lines(num_states);
for n_idx = 1:num_states
    En_sweep = (n_idx^2 * hp^2 ./ (8 * me * L_sweep.^2)) / q;
    plot(L_sweep_nm, En_sweep, 'LineWidth', 2, 'Color', colors(n_idx, :), ...
         'DisplayName', sprintf('n = %d', n_idx));
end
xlim([min(L_sweep_nm), max(L_sweep_nm)]);
ylim([0, 15]);
xlabel('Well Width L (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy E_n (eV)', 'FontSize', 12, 'FontWeight', 'bold');
title('4. Quantum Confinement Effect (E_n \propto 1/L^2)', 'FontSize', 13);
legend('Location', 'northeast');
grid on;

% Figure 5: Momentum Space Representation |phi_n(p)|^2 (Heisenberg Uncertainty)
figure('Color', 'w', 'Name', 'Momentum_Space_Distribution');
p_max = 5 * pi * hbar_val / L_val;
p_grid = linspace(-p_max, p_max, 500);
hold on;
for n_idx = 1:3
    kn = n_idx * pi / L_val;
    % Analytical Fourier transform of sin(kn * x) box state
    term_plus  = sinc((p_grid/hbar_val - kn) * L_val / (2*pi));
    term_minus = sinc((p_grid/hbar_val + kn) * L_val / (2*pi));
    phi_p = abs(term_plus - term_minus).^2;
    phi_p = phi_p / max(phi_p);
    plot(p_grid / (hbar_val / 1e-9), phi_p, 'LineWidth', 2, ...
         'DisplayName', sprintf('n = %d', n_idx));
end
xlabel('Momentum p (h_bar / nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Probability Density |\psi(p)|^2', 'FontSize', 12, 'FontWeight', 'bold');
title('5. Wavefunctions in Momentum Space (Fourier Domain)', 'FontSize', 13);
legend('Location', 'northeast');
grid on;

%% 4. Save Figures 
[script_folder, script_name] = fileparts(mfilename('fullpath'));
if isempty(script_name)
    script_name = 'Schrodinger_Laplace_infinite_well';
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

