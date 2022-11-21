function position = Laplace_kernel(noise)
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.
%%
% To find the outliers in noise.
%%
Noise = noise;
Noise = sort(Noise);
A_E = std(noise);%标准�?%A_E表示σE
R = (prctile(Noise,75) - prctile(Noise,25));
sigma = 1.06*min(A_E,R/1.34)*length(noise)^-0.2;
w = exp(-(abs(noise)./(sigma)));%P表示权重 before 2*pi*sqrt(A)
position = find(w<1e-20);  

end
