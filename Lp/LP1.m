function [Out_X] = LP1(M_Omega,Omega,rak,maxiter)

% This matlab code implements the l1-reg 
% method for Matrix Completion.
%
% M - m x n matrix of observations/data (required input)
%
% Omega is the subset
%
% maxIter - maximum number of iterations
%         - DEFAULT 500, if omitted or -1.

iter = 0;
[n1,n2] = size(M_Omega);
U = randn(n1,rak);
t = 2;

while 1
    iter = iter + 1;
    clear row col;
for j = 1:n2
    clear row;
    row_i = find(Omega(2,:) == j);
    [~,col_n] = size(row_i);
    for n = 1:col_n
        row(n,1) = Omega(1,row_i(1,n));
    end
    for i = 1:length(row)
        U_I(i,:) =  U(row(i,1),:);
        b_I(i) = M_Omega(row(i,1),j);
    end
    
    [~,W_col] = size(b_I);
    if iter == 1
        W = eye(W_col);
    else
        ksi = U_I * V(:,j) - b_I';
        clear W;
        for W_index = 1 : W_col
            W(W_index,W_index) = 1/((abs(ksi(W_index,1)))^2 + 0.0001)^(1/4);
        end
    end
    for inner_iter = 1 : t
        V(:,j) = pinv(U_I' * W' * W * U_I)* U_I' * W' * W * b_I';
        ksi = U_I * V(:,j) - b_I';
        clear W;
        for W_index = 1 : W_col
            W(W_index,W_index) = 1/((abs(ksi(W_index,1)))^2 + 0.0001)^(1/4);
        end
    end
    clear U_I b_I;
end
clear row col;
for i = 1:n1
    clear col;
    col_i_new = find(Omega(1,:) == i);
    [~,col_n] = size(col_i_new);
    for n = 1:col_n
        col (1,n) = Omega(2,col_i_new(1,n));
    end    
    for j =1:length(col)
        V_I(:,j) = V(:,col(1,j));
        b_I(j) = M_Omega(i,col(1,j));
    end
    [~,W_col] = size(b_I);
    ksi = U(i,:) * V_I - b_I;
    clear W;
    for W_index = 1 : W_col
        W(W_index,W_index) = 1/((abs(ksi(1,W_index)))^2 + 0.0001)^(1/4);
    end
    
    for inner_iter = 1 : t
        U(i,:) = b_I * W * W'* V_I' * pinv(V_I * W * W'* V_I');
        ksi = U(i,:) * V_I - b_I;
        clear W;
        for W_index = 1 : W_col
            W(W_index,W_index) = 1/((abs(ksi(1,W_index)))^2 + 0.0001)^(1/4);
        end
    end
    clear V_I b_I;  
end
    %**************
    %Judging
    %**************
    X = U*V;
    if (iter>maxiter)
        break;
    end
end
    Out_X = X;
end