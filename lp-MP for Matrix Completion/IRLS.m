function x = IRLS(a0,b0,p,max_iter,index)
% IRLS for min ||a*x-b||_p^p, where a is a column vector,x is scalar
% Author: Wen-Jun Zeng (wenjzeng@gmail.com)

a = [];
b = [];
for i = 1 : length(index)
    if index(i) == 1
        a = [a; a0(i)];
        b = [b; b0(i)];
    end
end

m = length(a);
% Initialize using LS solution
x = (a'*b)/(a'*a);

r = a*x-b;
f_old = norm(r,p)^p;
error = 1;
epsilon = 1e-8;
k = 0;
while(error>epsilon)
    k = k+1;
    w = abs(r).^((p-2)/2);
    wa = w.*a;
    wb = w.*b;
    x = (wa'*wb)/(wa'*wa);
    r = a*x-b;
    f_new = norm(r,p)^p;
    error = abs(f_new-f_old)/f_new;
    f_old = f_new;
    if(k>max_iter)
        break;
    end
end