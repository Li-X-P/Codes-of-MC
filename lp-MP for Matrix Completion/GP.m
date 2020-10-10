function [Y, MSE] = GP( M_Omega, maxiteration, array_Omega, p)
 %%
 % M_Omega:         The input incomplete image
 % maxiteration:    Maximum iteration
 % array_Omega:     The subset of Omega
 % p:               The p of lp-norm
 %%
[m ,n] = size(M_Omega);
X_j = zeros(m, n);
u_j = randn(m, 1);
v_j = randn(n, 1);
MSE = [];
in = 0;

for iter = 1 : maxiteration
    R_j = M_Omega - X_j.*array_Omega; 
    %% calculate u_j
    in2 = 0;
    while (in2 < 2)
%         a = norm(R_j - u_j * v_j.','fro')^2/norm(R_j,'fro')^2;
        a = lp_norm(R_j - u_j * v_j.',p)^(1/p)/lp_norm(R_j,p)^(1/p);
        for i = 1 : m
            g_l = R_j(i,:).';
            index = array_Omega(i,:).';
            u_j(i,1) = IRLS(v_j, g_l,p,20,index);
        end
        
        %% calculate v_j
        for i = 1 : n
            h_q = R_j(:,i);
            index = array_Omega(:,i);
            v_j(i,1) = IRLS(u_j, h_q, p,20,index);
        end
        b = lp_norm(R_j - u_j * v_j.',p)^(1/p)/lp_norm(R_j,p)^(1/p);
        if a-b < 0.00001
            in2 = in2 + 1;
        end
    end
    X_j = X_j + u_j * v_j.';
    %% judgement
    MSE = [MSE, lp_norm(M_Omega-X_j.*array_Omega,p)^(1/p)/lp_norm(M_Omega,p)^(1/p)];
    Y = X_j;
    if iter~=1
        step_MSE = MSE(iter-1) - MSE(iter);
        if step_MSE < 0.0005
            in = in  + 1;
        end
        if in > 1
            break;
        end
    end

end

end