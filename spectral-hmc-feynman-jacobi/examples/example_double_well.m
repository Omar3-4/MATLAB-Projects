% EXAMPLE_DOUBLE_WELL
% Run the Feynman-Jacobi framework on a symmetric double-well potential.
% Demonstrates quantum mechanical tunneling (instanton paths).
% Runtime on a standard laptop: ~60-120 seconds for D=2000, N=1500.

clear; clc;

% Add src/ and utils/ to MATLAB path if needed
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'utils'));

disp('========================================');
disp('Example 2: Symmetric Double-Well');
disp('========================================');

% Quartic Double-Well Potential
m = 1.0;
Vx  = @(x) 0.25 * (x.^2 - 1.0).^2;
dVx = @(x) x .* (x.^2 - 1.0);

% Boundary conditions forcing a transition
xa = -1.0; xb = 1.0; T = 5.0;
D = 2000; N = 1500; L_steps = 15; NO_Test = 300;

% Run solver
[Classical_Path, Quantum_Fluctuations, L_class, L_quantum_min] = ...
    Feynman_Jacobi(xa, xb, T, m, Vx, dVx, D, N, L_steps, NO_Test, true, 'Symmetric Double Well');

disp(' ');
disp('Results:');
fprintf('Classical Geodesic Length : %.4f\n', L_class);
fprintf('Min Quantum Path Length   : %.4f\n', L_quantum_min);
