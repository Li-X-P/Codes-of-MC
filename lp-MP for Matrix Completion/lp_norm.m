function y = lp_norm(X,p)
[m, n] = size(X);
sum = 0;
for i = 1 : m
    for j = 1 : n
        sum = sum + abs(X(i,j))^p;
    end
end
y = sum^(1/p);
end