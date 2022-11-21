function [Out_X, MSE  ] = L0_BCD_image(M, M_Omega, Omega_array, rak, maxiter)
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.
%%
% for image inpainting

%% determine outlier_set
miss_outlier_entires = find(M_Omega == 0);
miss_entires = find(Omega_array == 0);
outlier_set_zero = setdiff(miss_outlier_entires,miss_entires);
outlier_set_one = find(M_Omega == 1);
outlier_set = union(outlier_set_zero,outlier_set_one);

%% initialization
[r,c] = size(M_Omega);
U = randn(r,rak);
V = randn(rak,c);
S = zeros(size(M_Omega));
S(outlier_set) = M_Omega(outlier_set);
MSE = zeros(1, maxiter);


for iter = 1 : maxiter
    M_Omega_update = M_Omega - S;
    
%% Update V   
    for j = 1:c
        row = find(Omega_array(:,j) == 1);
        U_I =  U(row,:);
        b_I = M_Omega_update(row,j);
        V(:,j) = pinv(U_I)* b_I;
    end
    
%% Update U       
    for i = 1 : r
        col = find(Omega_array(i,:) == 1);
        V_I = V(:,col);
        b_I = M_Omega_update(i,col);
        U(i,:) = b_I * pinv(V_I);
    end
    
%% Update S    
    E = (M_Omega - U*V).*Omega_array;
    S = zeros(size(M_Omega));
    S(outlier_set) = E(outlier_set);
    
%%    
    X = U*V;
    MSE(1,iter) = norm(M - X,'fro').^2/(r*c);


end
Out_X = X;
end