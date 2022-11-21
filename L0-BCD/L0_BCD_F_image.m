function [Out_X, MSE ] = L0_BCD_F_image(M, M_Omega, rak, maxiter)
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.
%%
% Fast version for image inpainting 

%%
zero_entires = find(M_Omega == 0);
Omega_array = ones(size(M_Omega));
Omega_array(zero_entires) = 0;
one_entries = find(M_Omega == 1);
Omega_array(one_entries) = 0;

[r,c] = size(M_Omega);
U = randn(r,rak);
V = randn(rak,c);
MSE = zeros(1, maxiter);

%%
for iter = 1 : maxiter
%% Update V       
    for j = 1:c
        row = find(Omega_array(:,j) == 1);
        U_I =  U(row,:);
        b_I = M_Omega(row,j);
        V(:,j) = pinv(U_I)* b_I;
    end
%% Update U       
    for i = 1 : r
        col = find(Omega_array(i,:) == 1);
        V_I = V(:,col);
        b_I = M_Omega(i,col);
        U(i,:) = b_I * pinv(V_I);
    end
%% Update    
    X = U*V;
    MSE(1,iter) = norm(M - X,'fro').^2/(r*c);

end
Out_X = X;
end