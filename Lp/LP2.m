function [Out_X] = LP2(M_Omega,Omega,rak,maxiter)
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
while 1
    iter = iter + 1;
    clear row col V;
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
        V(:,j) = pinv(U_I)* b_I';
        clear U_I b_I;
    end
    clear row col U;
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
        U(i,:) = b_I * pinv(V_I);
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