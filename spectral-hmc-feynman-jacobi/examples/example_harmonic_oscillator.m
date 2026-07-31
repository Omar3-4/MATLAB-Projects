% EXAMPLE_HARMONIC_OSCILLATOR
% Run the Feynman-Jacobi framework on a quantum harmonic oscillator.
% This is the foundational test case with exact Gaussian fluctuations.
% Runtime on a standard laptop: ~30-60 seconds for D=1000, N=1000.

clear; clc;

% Add src/ and utils/ to MATLAB path if needed
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'utils'));

disp('========================================');
disp('Example 1: Harmonic Oscillator');
disp('========================================');

% Physical parameters
m = 1.0; omega = 1.0;
Vx  = @(x) 0.5 * m * (omega^2) * (x.^2);
dVx = @(x) m * (omega^2) * x;

% Boundary conditions and grid parameters
xa = -1.0; xb = 1.0; T = 2.0;
D = 1000; N = 1000; L_steps = 10; NO_Test = 200;

% Run solver with full diagnostics
[Classical_Path, Quantum_Fluctuations, L_class, L_quantum_min] = ...
    Feynman_Jacobi(xa, xb, T, m, Vx, dVx, D, N, L_steps, NO_Test, true, 'Harmonic Oscillator');

disp(' ');
disp('Results:');
fprintf('Classical Geodesic Length : %.4f\n', L_class);
fprintf('Min Quantum Path Length   : %.4f\n', L_quantum_min);
