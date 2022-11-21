function [noise,sigma_1] = Gaussian_noise(Signal,model,SNR)
%%
% If you use this code, please cite the following paper in your corresponding work. Thanks!
% X. P. Li, Z.-L. Shi, Q. Liu and H. C. So, "Fast robust matrix completion
% via ?0-norm minimization" IEEE Transactions on Cybernetics, 2022.
%%
Signal_size = size(Signal);
num = 1;
for i = 1 : length(Signal_size)
    num = num * Signal_size(i);
end
Signal_power = sum(abs(Signal(:)).^2) / num;
s_n = 10^(SNR/10);


switch model
    case 'GW'
        Noise_power = Signal_power/s_n;
        noise = sqrt(Noise_power) * randn( Signal_size);
        
    case 'GM'
        c_1 = 0.9;
        c_2 = 0.1;
        sigma_v_2 = Signal_power/ s_n ; % sigma_v_2 = c_1 * sigma_1^2 + c_2 * sigma_2^2 = 0.9*sigma_1^2 + 0.1*sigma_2^2; sigma_2 = 10 * sigma_1
        sigma_1 = sqrt(sigma_v_2 / (c_1 + 10^2*c_2));
        p = c_1/1;        % the percentage of outliers c_1/c_2
        sigma = [ 10*sigma_1 sigma_1 ];
        flag = binornd( 1, p, Signal_size );
        noise = sigma(1) * randn( Signal_size).*(ones(Signal_size) - flag) + sigma(2)*randn( Signal_size).*flag;
    case 'Laplacian'
        u1 = rand(Signal_size);
        u2 = rand(Signal_size);
        noise = log(u1./u2);
        noise = bsxfun(@minus, noise, mean(noise));
        noise = bsxfun(@rdivide, noise, std(noise));
        Noise_power = Signal_power/s_n;
        noise = sqrt(Noise_power)*noise;
end
end