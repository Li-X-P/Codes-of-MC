function [M_Omega, array_Omega] = mask_image(M)
    array_Omega = load('mask_image.mat');
    array_Omega = array_Omega.mask_image;
    array_Omega = imresize(array_Omega,size(M));
    array_Omega = round(array_Omega);
    M_Omega = M.*array_Omega;
end