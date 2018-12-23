clear all;

load('E:\Code\Data\M.mat');
[m, n] = size(M);
rak = rank(M);
p= 0.45;
array_Omega = binornd( 1, p, [ m, n ] );
%sampling
M_Omega = M .*array_Omega;

%add noise
dB = 6;
s_n = 10^(dB/10);    
length_Omega = m*n*p;
sigma_v_2 = norm(M_Omega,'fro')^2/(length_Omega * s_n) ;
sigma_1 = sqrt(sigma_v_2 / 10.9);
model = 'GM';
switch model
    case 'GW'
        sigma_1 = 0.1;
    case 'SB'
        sigma_1 = 5;
    case 'GM'
        dB = 6;
        s_n = 10^(dB/10);
        sigma_v_2 = (norm(M_Omega,'fro')^2)/(length_Omega * s_n) ;
        sigma_1 = sqrt(sigma_v_2 / 10.9);     
end

Q = Gaussian_noise(sigma_1,model);
addnoise = M + Q;
% 
% inexact alm mc
tic;
[A, E] = inexact_alm_mc(M.*array_Omega, 1e-6, 500);
tElapsed = toc;
% semilogy(RMSE);
% xlabel('Iteration number');
% ylabel('Nomalized RMSE');
% legend('IALM');

result = A.U *A.V';
error = norm(result - M, 'fro')/norm(M,'fro');

disp('Iteration');
disp(iter);
disp('Time');
disp(tElapsed);
disp('rank of A');
disp(svp);
disp('|A-M|_F / |M|_F');
disp(error);

