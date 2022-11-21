function [Out_X, MSE,loss] = L0_BCD_F(M, M_Omega, Omega_array, rak, maxiter)
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.
%%
% Fast version - update Omega_array 

%%
[r,c] = size(M_Omega);
U = randn(r,rak);
V = randn(rak,c);
MSE = zeros(1, maxiter);
loss = zeros(1,maxiter);
g = @(U,V,S,thres) norm((M_Omega-U*V-S).*Omega_array,'fro')^2+ thres*nnz(S.*Omega_array);
thres = 10000;


%% initialization
Omega = find(Omega_array==1);
posi = initializationS(M_Omega(Omega));
outliers_posi = Omega(posi);
outliers = M_Omega(outliers_posi);
S = zeros(size(M_Omega));
S(outliers_posi) = outliers;

Omega_array_update = Omega_array;
Omega_array_update(outliers_posi) = 0;


%%
for iter = 1 : maxiter

    M_Omega_update = (M_Omega - S).*Omega_array_update;
%% Update V   
    for j = 1:c
        row = find(Omega_array_update(:,j) == 1);
        U_I =  U(row,:);
        b_I = M_Omega_update(row,j);
        V(:,j) = (U_I'*U_I)\(U_I'* b_I);
    end
%% Update U       
    for i = 1 : r
        col = find(Omega_array_update(i,:) == 1);
        V_I = V(:,col);
        b_I = M_Omega_update(i,col);
        U(i,:) = b_I * V_I' /(V_I*V_I');
    end
    
%% Update S    
    Y_Omega = (M_Omega - U*V).*Omega_array;
    Omega = find(Omega_array==1);
    posi = Laplace_kernel(abs(Y_Omega(Omega)));
    outliers_posi = Omega(posi);
    outliers = abs(Y_Omega(outliers_posi));
    temp_thres = thres;
    thres = min(min(outliers(:))^2,temp_thres);
    s = reshape(Y_Omega(Omega),[],1).* (reshape(abs(Y_Omega(Omega)),[],1) >= sqrt(thres));
    S = zeros(size(Y_Omega));
    S(Omega) = s;
    
%% update Omega_array    
    Omega_array_update = Omega_array;
    Omega_array_update(outliers_posi) = 0;
    
%% 
    X = U*V;
    MSE(1,iter) = norm(M - X,'fro').^2/(r*c);
    loss(iter) = g(U,V,S,thres);
%         if norm(M - X,'fro')^2/(r*c) < 0.1
%             break
%         end

end
Out_X = X;
end


function y = initializationS(input_data)
%%
% To find the outliers in input.
%%
Input_data = input_data;
Input_data = sort(Input_data);
A_E = std(input_data);%标准差%A_E表示σE
R = (prctile(Input_data-mean(input_data),75) - prctile(Input_data-mean(input_data),25))/1;
A = 1.06*min(A_E,R/1.34)*length(input_data)^-0.2;
w = exp(-(abs(input_data-mean(input_data))./(A)));%P表示权重
y = find(w<1e-20);  
end