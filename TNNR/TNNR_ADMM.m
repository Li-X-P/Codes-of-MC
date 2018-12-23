function [X] = TNNR_ADMM(Omega,M,M_Omega, beta, tol, max_iter)
% Inputs:
%   Omega - set of observed entries.  Should be linearly indexed.
%
%	M - the original matrix
%
%	M_Omega - the sampled matrix
%
%	beta - defult 0.1
%
%   tol - stopping criteria
%
%   maxiter - maximum number of iterations
% TNNR-ADMM
[m,n] = size(M);
r = rank(M);
X0 = zeros(m, n);
X0(Omega) = M_Omega;
X = X0;
W = X;
Y = X;
for k = 1: max_iter
    [u, s, v] = svd(X);
    u = u(:, 1: r);
    v = v(:, 1: r);
    X = D_shrinke(W- Y/beta, 1/beta);
    W = X+ (u* v'+ Y)/beta;
    W(Omega) = M_Omega;
    Y = Y+ beta*(X- W);
    X0 = X;
    RMSE(1, k) = norm(X - M, 'fro')/norm(M, 'fro');
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
