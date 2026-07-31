function L = Jacobi_metric(x, E, m, Vx)
%JACOBI_METRIC Computes the Jacobi-Maupertuis arc length of a path.
%
%   L = JACOBI_METRIC(x, E, m, Vx) computes the geometric length of a path
%   x in configuration space under the Jacobi conformal metric:
%       ds = sqrt(2*m*(E - V(x))) * |dx|
%
%   This metric emerges from Maupertuis-Jacobi principle and maps classical
%   mechanics onto Riemannian geometry.
%
%   Inputs:
%       x   - Path vector (1 x D or D x 1)
%       E   - Total conserved energy (scalar)
%       m   - Particle mass (scalar)
%       Vx  - Function handle for potential V(x)
%
%   Output:
%       L   - Scalar Jacobi arc length
%
%   References:
%       Jacobi, C.G.J., "Vorlesungen uber Dynamik" (1866), p. 44.
%       Feynman, R.P., Rev. Mod. Phys. 20, 371 (1948).

    dx = diff(x);
    x_mid = (x(1:end-1) + x(2:end)) / 2;
    V = Vx(x_mid);
    
    % Numerical safeguard: ensure non-negative argument under sqrt
    K_Term = max(0, E - V + eps);
    
    ds = sqrt(2 * m * K_Term) .* abs(dx);
    L = sum(ds);
end
