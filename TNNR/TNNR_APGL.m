function [X1] = TNNR_APGL(Omega, M, M_Omega, lambda, tol, max_iter)
% Inputs:
%   Omega - set of observed entries.  Should be linearly indexed.
%
%	M - the original matrix
%
%	M_Omega - the sampled matrix
%
%	lambda - defult 0.06
%
%   tol - stopping criteria
%
%   maxiter - maximum number of iterations
%% TNNR-APGL
[m,n] = size(M);
r = rank(M);
t = 1;
X0 = zeros(m, n);
X0(Omega) = M_Omega;
X = X0;
Y = X;
for k = 1: max_iter
    [u, s, v] = svd(X);
    u = u(:, 1: r);
    v = v(:, 1: r);
    X1 = D_shrinke(Y+ t*(u*v'- lambda*(Proj(Y, Omega)- X0)), t);
    t1 = (1+ sqrt(1+4*t^2))/2;
    Y = X1+ (t-1)*(X1- X)/t1;
    t = t1;
    X = X1;
    RMSE(1,k) = norm(X - M, 'fro')/norm(M, 'fro');
    if (RMSE(1,k) < tol)
        break;
    end
end
end

function D = D_shrinke(X, beta)
[u, s, v] = svd(X, 'econ') ;
s = diag(max(diag(s) - beta, 0));
D = u*s*v';
end


function P = Proj(X, index)
P = zeros(size(X));
P(index) = X(index);
end
