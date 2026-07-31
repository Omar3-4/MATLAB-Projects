% EXAMPLE_ASYMMETRIC_WELL
% Run the Feynman-Jacobi framework on an asymmetric anharmonic well.
% Broken reflection symmetry creates competing local and global minima.
% Runtime on a standard laptop: ~90-180 seconds for D=3000, N=2000.

clear; clc;

% Add src/ and utils/ to MATLAB path if needed
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src', 'utils'));

disp('========================================');
disp('Example 3: Asymmetric Anharmonic Well');
disp('========================================');

% Asymmetric Anharmonic Potential
m = 1.0;
Vx  = @(x) 0.25*(x.^4) - 0.5*(x.^2) + 0.2*x;
dVx = @(x) (x.^3) - x + 0.2;

% Boundary conditions across asymmetric barrier
xa = -1.2; xb = 0.8; T = 6.0;
D = 3000; N = 2000; L_steps = 20; NO_Test = 500;

% Run solver
[Classical_Path, Quantum_Fluctuations, L_class, L_quantum_min] = ...
    Feynman_Jacobi(xa, xb, T, m, Vx, dVx, D, N, L_steps, NO_Test, true, 'Asymmetric Anharmonic Well');

disp(' ');
disp('Results:');
fprintf('Classical Geodesic Length : %.4f\n', L_class);
fprintf('Min Quantum Path Length   : %.4f\n', L_quantum_min);
