
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.

%% Image
clear variables
close all hidden
image = imread('Windows.jpg');
[width,height,z]=size(image);
if(z>1)
    image=rgb2gray(image);
end
% unit8 to double
M = mat2gray(image);
[ r, c ] = size(M);


%% mask
% random
per = 0.5;
array_Omega = binornd( 1, per, [ r, c ] );
M_Omega = M.*array_Omega;

% fixed
% [M_Omega, array_Omega] = mask_image(M);


%% add noise
omega = find(array_Omega(:)==1);
dB = 5;
noise = imnoise(M_Omega(omega),'salt & pepper',1/dB) ;
Noise = zeros(size(M_Omega));
Noise(omega) = noise; 
M_Omega = Noise;


rak = 10;
maxiter = 50;

%%
figure
imshow(M_Omega)
tic
[X_F, MSE_F] = L0_BCD_F_image(M, M_Omega, rak, maxiter);
toc
fprintf('PSNR:%d; SSIM:%d\n',psnr(X_F,M),ssim(X_F,M))
tic
[X, MSE] = L0_BCD_image(M,M_Omega,array_Omega,rak, maxiter);
toc
fprintf('PSNR:%d; SSIM:%d\n',psnr(X,M),ssim(X,M))


