clear variables
close all hidden

image = imread('7.jpg');
[width,height,z]=size(image);
if(z>1)
    image=rgb2gray(image);
end
% unit8 to double
image = mat2gray(image);

[m,n] = size(image);
M = image;
per = 0.6;
array_Omega = binornd( 1, per, [ m, n ] );

%%
DB = 9;
M_Omega = M.* array_Omega;
M_noise = imnoise(M_Omega, 'salt & pepper', 1/DB);
[Y, MSE] = GP(M_noise, 100, array_Omega, 1);
[peaksnr, snr] = psnr(Y, M);
fprintf('\n The Peak-SNR value is %0.4f\n', peaksnr);


semilogy(MSE);
xlabel('Iteration number');
ylabel('MSE');

hs = self_subplot(1, 4, [0.005, 0.005], [0.01, 0.01], [0.01, 0.01]);
axes(hs(1));
imshow(M);
title('Original','fontname','Times New Roman');
axes(hs(2));
imshow(M_Omega);
title('With  missing data','fontname','Times New Roman');
axes(hs(3));
imshow(M_noise);
title('With  missing data + noise','fontname','Times New Roman');
axes(hs(4));
imshow(Y);
title('Recovery','fontname','Times New Roman');
