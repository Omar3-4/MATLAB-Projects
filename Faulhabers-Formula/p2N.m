function [f S_sum] = p2N(p,N_value)
     % p2N is a Matlab function used to find the polynomial 
     % equivalent to the "Sum of Powers" p as a function of N
     % which is called "Faulhaber's Formula"
     
     % p: the power of the integers
     % S: row vector of the cumulative sums
     % M: n*n matrix (Vandermonde-like matrix)
     % V: symbolic row vector containing powers of N

% Startup
n = p+1;
S = [];
M = [];
V = [];
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
% Solve using symbolic conversion to get exact fractions
K = sym(M) \ sym(S(2:end)).';     % Mat(n*n)*vec(n*1) = vec(n*1) % K [column]
f = V*K;

S_sum =double(subs(f,N,N_value));
end

