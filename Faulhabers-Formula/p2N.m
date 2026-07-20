function [f, S_sum] = p2N(p, N_value)
     % p2N is a Matlab function used to find the polynomial 
     % equivalent to the "Sum of Powers" p as a function of N
     % which is called "Faulhaber's Formula"
     
     % p: the power of the integers
     % S: row vector of the cumulative sums
     % M: n*n matrix (Vandermonde-like matrix)
     % V: symbolic row vector containing powers of N

% Input Validation: Check if p is provided, non-empty, and valid
if nargin < 1 || isempty(p) || ~isnumeric(p) || ~isscalar(p) || p < 0 || isnan(p) || isinf(p)
    error('Input Validation: p must be a valid non-negative numeric finite scalar.');
end

% Input Validation: Check if N_value is provided and valid
if nargin < 2 || isempty(N_value) || ~isnumeric(N_value)
    % Fallback: Assign default 0 to prevent downstream substitution crash
    N_value = 0; 
end

% Startup
n = p+1;

% Preallocate memory for S, M, and V arrays to prevent dynamic resizing overhead
S = zeros(1, n + 1);
M = zeros(n, n);
V = sym(zeros(1, n));

S(1) = 0;

% 1st loop: Calculate sums for first n points 
for i = 1:n
    S(i+1) = S(i) + i^p;     % vec (1*n) [row]
end

% 2nd loop: Build the system of equations matrix
for i = 1:n
    for j = 1:n
        M(i,j) = i^j;     %Mat(n*n)
    end
end

syms N
V = N.^(1:n);

% Defensive Check: Ensure matrix M is well-conditioned/invertible before solving
if rcond(M) < 1e-15
    f = []; S_sum = []; return;
end

% Solve using symbolic conversion to get exact fractions
K = sym(M) \ sym(S(2:end)).';     % Mat(n*n)*vec(n*1) = vec(n*1) % K [column]
f = V*K;

% Defensive Check: Ensure substitution components match structural types
if ~isempty(f) && isa(f, 'sym')
    S_sum = double(subs(f, N, N_value));
else
    S_sum = NaN;
end
end
