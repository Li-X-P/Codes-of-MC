%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.

%%
clear variables
close all hidden
[r, c, rak] = deal(400,500,10);
M = randn(400, 10)* randn(10,500);

per = 0.5;  % Obervation ratio
dB =3;      % SNR
%%

array_Omega = binornd( 1, per, [ r, c ] );
M_Omega = M.*array_Omega;
omega = find(array_Omega(:)==1);
noise = Gaussian_noise(M_Omega(omega),'GM',dB);
Noise = zeros(size(M_Omega));
Noise(omega) = noise;
M_Omega = M_Omega + Noise;
maxiter = 50;



[X_F,MSE_F,loss_F] = L0_BCD_F(M,M_Omega,array_Omega,rak, maxiter);


[X, MSE, loss] = L0_BCD(M,M_Omega,array_Omega,rak, maxiter);




figure
semilogy(MSE_F,'r','LineWidth',1.2);
hold on
semilogy(MSE,'k--','LineWidth',1.2);

xlabel('Iteration number');
ylabel('MSE');
legend('$\ell_0$-BCD-F','$\ell_0$-BCD','Interpreter','LaTex');
figure_setting(1.5, 15, 600, 500)



