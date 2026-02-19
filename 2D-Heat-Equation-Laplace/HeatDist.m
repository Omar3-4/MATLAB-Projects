function T = HeatDist(N, T_top, T_bottom, T_left, T_right)
% HeatDist Solves the 2D Steady-State Heat Equation (Laplace Equation).
%   T = HeatDist(N, T_top, T_bottom, T_left, T_right) calculates the 
%   temperature distribution on a square plate of size N-by-N using the 
%   Finite Difference Method (Gauss-Seidel iteration).
%
%   INPUTS:
%       N        - Number of grid points along one side (Resolution)&
%       Dimensions of the plate which N is the square length.
%       T_top    - Constant temperature boundary at the top edge.
%       T_bottom - Constant temperature boundary at the bottom edge.
%       T_left   - Constant temperature boundary at the left edge.
%       T_right  - Constant temperature boundary at the right edge.
%
%   OUTPUT:
%       T        - A 2D matrix representing the calculated temperature at 
%                  each grid point.

% Boundary Conditions
T = zeros(N,N);
T(1,:) =  T_top;
T(N,:)= T_bottom;
T(2:N-1,1) = T_right;
T(2:N-1,N) = T_left; 

% Convergence Settings
error = 1;
tol = 1e-6;

% Iteration Loop
while error>tol
    T_old = T;
    for i = 2:length(T)-1
        for j=2:length(T)-1
            T(i,j) = 1/4 *(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1));
        end
    end
    error = max(abs(T(:)-T_old(:)));
end

% ---- plotting ----
subplot(1,2,1);
imagesc(T);
title('2D Temperature Distribution');
colorbar; xlabel('x'); ylabel('y');
subplot(1,2,2);
meshc(T);
shading interp;
title('3D Surface Plot');
colorbar; xlabel('x'); ylabel('y'); zlabel('Temp');
end
